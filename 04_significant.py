import pandas as pd
from pathlib import Path
from glob import glob
import os
from tqdm import tqdm

print("Collecting significant pairs...")
sig_pairs = []

for filepath in sorted(glob("filtered_celltypes/*.tsv")):
    df = pd.read_csv(filepath, sep="\t", usecols=["snp_id", "gene_id"])
    sig_pairs.append(df)

significant_df = pd.concat(sig_pairs).drop_duplicates()
print(f"  > Collected {len(significant_df)} unique significant SNP–gene pairs")

os.makedirs("celltypes", exist_ok=True)

print("Filtering celltypes to match significant pairs...")

for filepath in sorted(glob("unfiltered_celltypes/*.tsv")):
    print(f"  * Processing {filepath}...")
    region = Path(filepath).stem
    output_path = Path("celltypes") / f"{region}.tsv"

    filtered_chunks = []

    with open(filepath) as f:
        total_lines = sum(1 for _ in f) - 1

    chunk_reader = pd.read_csv(filepath, sep="\t", chunksize=500_000)

    for chunk in tqdm(
        chunk_reader,
        total=(total_lines // 500_000) + 1,
        desc=f"  > Filtering {region}",
        unit="chunk",
    ):
        filtered = chunk.merge(significant_df, on=["snp_id", "gene_id"], how="inner")
        filtered_chunks.append(filtered)

    final_df = pd.concat(filtered_chunks, ignore_index=True)
    final_df.to_csv(output_path, sep="\t", index=False, float_format="%.6g")
    print(f"  - {region}.tsv created with {len(final_df)} entries")

print("All done. Filtered files are in 'celltypes/'.")
