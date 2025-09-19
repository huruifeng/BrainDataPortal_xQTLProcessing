## sampling the data records from the dataset
import os
import time
import pandas as pd
from glob import glob
from pathlib import Path

unfiltered_files = qtl_data_files = sorted(glob("data/eQTLsummary_unfiltered/*_combined_nominals.tsv"))
output_dir = "data/eQTLsummary_sampled"

print(f"Found {len(unfiltered_files)} unfiltered celltype files.")
os.makedirs(output_dir, exist_ok=True)

## for each cell type, sample 1M records, save to new file
print("Sampling 100k records from each unfiltered celltype file...")
for i, filepath in enumerate(unfiltered_files):
    path = Path(filepath)
    celltype = path.stem.split("_")[0]
    output_file = output_dir + f"/{celltype}_sampled.tsv"

    print(f"  [{i+1}/{len(unfiltered_files)}] Sampling from {celltype}...")
    file_start = time.time()

    # Read the data file
    df = pd.read_csv(filepath, sep="\t")

    # Sample 1m records (or all if less than 1M)
    sampled_df = df.sample(n=min(1000000, len(df)), random_state=42)

    # Save the sampled data to a new file
    sampled_df.to_csv(output_file, sep="\t", index=False)

    file_end = time.time()
    print(f"    Sampled {len(sampled_df)} records to {output_file} in {file_end - file_start:.2f} seconds.")