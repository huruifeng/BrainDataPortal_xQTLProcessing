import polars as pl
from pathlib import Path
from glob import glob
import os
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import time

dataset_name = "eQTLsummary_demo"

print("Collecting significant pairs...")
sig_dfs = []
for filepath in tqdm(sorted(glob(dataset_name + "/filtered_celltypes/*.tsv")), desc="Reading filtered files"):
    df = pl.read_csv(filepath, separator="\t", columns=["snp_id", "gene_id"])
    sig_dfs.append(df)

significant_df = pl.concat(sig_dfs).unique()
print(f"  > Collected {significant_df.height} unique significant SNP–gene pairs")

os.makedirs(dataset_name + "/celltypes", exist_ok=True)
unfiltered_files = sorted(glob(dataset_name + "/unfiltered_celltypes/*.tsv"))

print("Filtering celltypes to match significant pairs...")

def process_file(filepath):
    region = Path(filepath).stem
    output_path = f"{dataset_name}/celltypes/{region}.tsv"
    print(f"  > Starting {region}.tsv...")

    start_time = time.time()
    chunk_size = 1_000_000
    header = pl.read_csv(filepath, separator="\t", n_rows=0)
    total_rows = sum(1 for _ in open(filepath)) - 1  # -1 for header

    chunks = []
    for offset in tqdm(range(1, total_rows + 1, chunk_size), desc=f"  Processing {region}", unit="row"):
        df_chunk = pl.read_csv(
            filepath,
            separator="\t",
            skip_rows=offset,
            n_rows=chunk_size,
            has_header=False,
            new_columns=header.columns
        )
        df_filtered = df_chunk.join(significant_df, on=["snp_id", "gene_id"], how="semi")
        chunks.append(df_filtered)

    df = pl.concat(chunks) if chunks else pl.DataFrame(schema=header.schema)
    df.write_csv(output_path, separator="\t", float_precision=6)

    elapsed = time.time() - start_time
    return (region, df.height, elapsed)

results = []
for filepath in tqdm(unfiltered_files, desc="Processing files", unit="file"):
    results.append(process_file(filepath))

for region, count, elapsed in results:
    print(f"  - {region}.tsv: {count} entries ({elapsed:.1f}s)")

## move the celltype_parquet.json to celltypes
os.rename(dataset_name + "/celltype_parquet.json", dataset_name + "/celltypes/celltype_parquet.json")


print("All done. Filtered files are in 'celltypes/'.")
