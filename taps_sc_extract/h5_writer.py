"""
Amethyst-compatible HDF5 writer with single-pass and sharded directory support.

Writes per-barcode base-resolution methylation datasets with the schema:
  - /CG/<barcode>/1
  - /CH/<barcode>/1
dtype: [('chr', 'S10'), ('pos', '<i8'), ('pct', '<f8'), ('t', '<i8'), ('c', '<i8')]
Note: 't' (unmethylated) comes before 'c' (methylated).

Supports:
1. Incremental append mode (append_data).
2. Single-pass direct creation mode (create_cell_dataset) for 50x faster writing.
3. Multi-file sharded directory writing (write_sharded_hdf5_dir) with portable relative ExternalLink master.h5.
"""

import concurrent.futures
import hashlib
import os
from collections import defaultdict
from typing import Any, Dict, List, Optional, Union

import h5py
import numpy as np

# Disable HDF5 file locking on network / mounted filesystems
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

SUPPORTED_COMPRESSION = ("gzip", "lzf", "none", "blosc", "blosc-zstd")


def h5_compression_kwargs(compression: str, compression_level: int = 1) -> Dict[str, Any]:
    """Return h5py ``create_dataset`` compression keyword arguments."""
    if compression in (None, "", "none"):
        return {}
    if compression == "gzip":
        return {"compression": "gzip", "compression_opts": int(compression_level)}
    if compression == "lzf":
        return {"compression": "lzf"}
    if compression in ("blosc", "blosc-zstd"):
        try:
            import hdf5plugin
        except ImportError as exc:
            raise ImportError(
                "compression=%r requires the hdf5plugin package (pip install hdf5plugin)"
                % compression
            ) from exc
        cname = "zstd" if compression == "blosc-zstd" else "lz4"
        clevel = int(compression_level)
        if not 0 <= clevel <= 9:
            clevel = 1
        return dict(
            hdf5plugin.Blosc(
                cname=cname,
                clevel=clevel,
                shuffle=hdf5plugin.Blosc.SHUFFLE,
            )
        )
    raise ValueError(
        f"Unsupported HDF5 compression {compression!r}. "
        f"Use one of: {', '.join(SUPPORTED_COMPRESSION)}."
    )


def configure_blosc_threads(
    n_writers: int,
    compression: str,
    compression_threads: Optional[int] = None,
) -> Optional[int]:
    """
    Set BLOSC_NTHREADS for intra-chunk parallel compression.

    Blosc's thread pool is process-global, so when several shard writers share
    a process the per-compress thread count is ``cpus // n_writers`` (unless
    ``compression_threads`` is set explicitly).
    """
    if compression not in ("blosc", "blosc-zstd"):
        return None
    if compression_threads is not None:
        n = max(1, int(compression_threads))
    else:
        n = max(1, (os.cpu_count() or 1) // max(1, n_writers))
    os.environ["BLOSC_NTHREADS"] = str(n)
    return n

# Amethyst / Facet structured array dtype
METH_DTYPE = np.dtype([
    ('chr', 'S10'),
    ('pos', '<i8'),
    ('pct', '<f8'),
    ('t', '<i8'),
    ('c', '<i8')
])


class AmethystH5Writer:
    """Writer for Amethyst-compatible HDF5 files."""

    def __init__(
        self,
        filepath: str,
        mode: str = 'w',
        compression: str = 'gzip',
        compression_level: int = 6,
    ):
        self.filepath = filepath
        self.compression = compression
        self.compression_level = compression_level
        self.h5 = h5py.File(filepath, mode)
        # Ensure CG and CH groups exist
        self.cg_group = self.h5.require_group('CG')
        self.ch_group = self.h5.require_group('CH')
        self.write_metadata()

    def write_metadata(self, version: str = "amethyst2.0.0"):
        """Write metadata version attribute/dataset."""
        meta_group = self.h5.require_group('metadata')
        if 'version' not in meta_group:
            meta_group.create_dataset('version', data=version.encode('utf-8'))

    def create_cell_dataset(self, context: str, barcode: str, records: np.ndarray):
        """
        Create a complete methylation dataset for a barcode and context in a single pass.

        Bypasses dataset resizing and B-tree allocation churn.
        `records` must be a coordinate-sorted NumPy structured array with dtype METH_DTYPE.
        """
        if len(records) == 0:
            return

        group = self.cg_group if context == 'CG' else self.ch_group
        bc_group = group.require_group(barcode)

        if '1' in bc_group:
            del bc_group['1']

        chunk_size = min(len(records), 65536)
        bc_group.create_dataset(
            '1',
            data=records,
            chunks=(chunk_size,),
            dtype=METH_DTYPE,
            **h5_compression_kwargs(self.compression, self.compression_level),
        )

    def append_data(self, context: str, barcode: str, records: np.ndarray):
        """
        Append structured methylation records for a given barcode and context ('CG' or 'CH').

        `records` must be a numpy array with dtype METH_DTYPE.
        """
        if len(records) == 0:
            return

        group = self.cg_group if context == 'CG' else self.ch_group
        bc_group = group.require_group(barcode)

        if '1' not in bc_group:
            chunk_size = max(min(len(records), 32768), 1024)
            bc_group.create_dataset(
                '1',
                data=records,
                maxshape=(None,),
                chunks=(chunk_size,),
                dtype=METH_DTYPE,
                **h5_compression_kwargs(self.compression, self.compression_level),
            )
        else:
            ds = bc_group['1']
            curr_len = ds.shape[0]
            new_len = curr_len + len(records)
            ds.resize((new_len,))
            ds[curr_len:new_len] = records

    def close(self):
        """Close the HDF5 file."""
        if self.h5:
            self.h5.close()
            self.h5 = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def _write_single_shard(
    shard_filepath: str,
    shard_barcodes: List[str],
    cg_dict: Dict[str, np.ndarray],
    ch_dict: Dict[str, np.ndarray],
    compression: str = "gzip",
    compression_level: int = 1,
    version: str = "amethyst2.0.0",
):
    """Write a single shard file containing a subset of barcodes."""
    with AmethystH5Writer(shard_filepath, mode="w", compression=compression, compression_level=compression_level) as writer:
        writer.write_metadata(version)
        for bc in shard_barcodes:
            if bc in cg_dict:
                writer.create_cell_dataset("CG", bc, cg_dict[bc])
            if bc in ch_dict:
                writer.create_cell_dataset("CH", bc, ch_dict[bc])


def write_sharded_hdf5_dir(
    output_dir: str,
    cg_data: Dict[str, np.ndarray],
    ch_data: Dict[str, np.ndarray],
    n_shards: int = 16,
    compression: str = "gzip",
    compression_level: int = 1,
    version: str = "amethyst2.0.0",
) -> str:
    """
    Write cell methylation datasets partitioned across multiple shard HDF5 files in a directory.

    Creates:
      - `output_dir/shard_00.h5`, `output_dir/shard_01.h5`, ...
      - `output_dir/master.h5` containing relative ExternalLinks to all shards.

    Returns:
      Path to `master.h5` (compatible with `amethyst::createObject()`).
    """
    os.makedirs(output_dir, exist_ok=True)
    all_barcodes = sorted(set(cg_data.keys()) | set(ch_data.keys()))

    # Partition barcodes deterministically across shards
    shard_assignments: Dict[int, List[str]] = defaultdict(list)
    for bc in all_barcodes:
        shard_idx = int(hashlib.md5(bc.encode("ascii")).hexdigest(), 16) % n_shards
        shard_assignments[shard_idx].append(bc)

    shard_filenames = [f"shard_{i:03d}.h5" for i in range(n_shards)]

    n_writers = min(n_shards, 16)
    configure_blosc_threads(n_writers, compression)
    # Write shard files in parallel using ThreadPool
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_writers) as executor:
        futures = []
        for i in range(n_shards):
            shard_bcs = shard_assignments[i]
            shard_path = os.path.join(output_dir, shard_filenames[i])
            futures.append(
                executor.submit(
                    _write_single_shard,
                    shard_path,
                    shard_bcs,
                    cg_dict=cg_data,
                    ch_dict=ch_data,
                    compression=compression,
                    compression_level=compression_level,
                    version=version,
                )
            )
        for f in concurrent.futures.as_completed(futures):
            f.result()

    # Create master.h5 with relative external links
    master_path = os.path.join(output_dir, "master.h5")
    with h5py.File(master_path, "w") as master_h5:
        meta_group = master_h5.require_group("metadata")
        meta_group.create_dataset("version", data=version.encode("utf-8"))
        cg_master = master_h5.require_group("CG")
        ch_master = master_h5.require_group("CH")

        for i in range(n_shards):
            shard_bcs = shard_assignments[i]
            if not shard_bcs:
                continue
            shard_filename = shard_filenames[i]
            for bc in shard_bcs:
                if bc in cg_data:
                    cg_master[bc] = h5py.ExternalLink(shard_filename, f"CG/{bc}")
                if bc in ch_data:
                    ch_master[bc] = h5py.ExternalLink(shard_filename, f"CH/{bc}")

    return master_path
