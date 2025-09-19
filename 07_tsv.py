import polars as pl
from pathlib import Path


parquet_files = (
    list(Path("celltypes").glob("*.parquet"))
    + list(Path("gene_locations").glob("*.parquet"))
    + list(Path("snp_locations").glob("*.parquet"))
)

for parquet_file in parquet_files:
    print(f"Converting {parquet_file}...")
    try:
        # df = pl.read_csv(tsv_file, separator="\t")
        df = pl.read_parquet(parquet_file)
        tsv_path = parquet_file.with_suffix(".tsv")
        df.write_csv(tsv_path, separator="\t")
        print(f"Saved to {tsv_path}")
        parquet_file.unlink()
        print(f"Deleted original TSV file: {parquet_file}")
    except Exception as e:
        print(f"Failed to convert {parquet_file}: {e}")
