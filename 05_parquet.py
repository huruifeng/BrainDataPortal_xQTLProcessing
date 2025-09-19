import polars as pl
from pathlib import Path

dataset_name = "eQTLsummary_demo"

# List all TSV files in the specified directories
tsv_files = (
    list(Path(f"{dataset_name}/celltypes").glob("*.tsv"))
    + list(Path(f"{dataset_name}/gene_locations").glob("*.tsv"))
    + list(Path(f"{dataset_name}/snp_locations").glob("*.tsv"))
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
