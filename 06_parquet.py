import polars as pl
from pathlib import Path


tsv_files = (
    list(Path("celltypes").glob("*.tsv"))
    + list(Path("gene_locations").glob("*.tsv"))
    + list(Path("snp_locations").glob("*.tsv"))
)

for tsv_file in tsv_files:
    print(f"Converting {tsv_file}...")
    try:
        df = pl.read_csv(tsv_file, separator="\t")
        parquet_path = tsv_file.with_suffix(".parquet")
        df.write_parquet(parquet_path)
        print(f"Saved to {parquet_path}")
        tsv_file.unlink()
        print(f"Deleted original TSV file: {tsv_file}")
    except Exception as e:
        print(f"Failed to convert {tsv_file}: {e}")
