import pandas as pd
from pathlib import Path
from glob import glob
import os
from tqdm import tqdm

os.makedirs("filtered_celltypes", exist_ok=True)

files = sorted(glob("unfiltered_celltypes/*.tsv"))

print("Filtering for significant QTLs in each celltype...")

for filepath in tqdm(files, desc="Processing files", unit="file"):
    region = Path(filepath).stem
    output_path = Path("filtered_celltypes") / f"{region}.tsv"

    df = pd.read_csv(
        filepath,
        sep="\t",
        dtype={"snp_id": str, "gene_id": str, "p_value": float, "beta_value": float},
    )
    filtered_df = df[df["p_value"] < 0.01]

    filtered_df.to_csv(output_path, sep="\t", index=False, float_format="%.6g")

print("Done. Filtered files are in 'filtered_celltypes/'.")
