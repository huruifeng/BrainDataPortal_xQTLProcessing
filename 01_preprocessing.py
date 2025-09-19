import os
import pandas as pd
from pathlib import Path
from glob import glob
import time
import re
from tqdm import tqdm


def safe_filename(name):
    return re.sub(r"[^a-zA-Z0-9_\-]", "_", name)


def delete_old_data(dataset_name):
    directories = [
        "celltypes",
        "filtered_celltypes",
        "gene_jsons",
        "gene_locations",
        "snp_json",
        "snp_locations",
        "unfiltered_celltypes",
    ]
    for directory in directories:
        if os.path.exists(dataset_name + "/" + directory):
            print(f"Deleting old data in {dataset_name }/{directory}...")
            for file in glob(os.path.join(dataset_name , directory, "*")):
                try:
                    os.remove(file)
                except Exception as e:
                    print(f"  Error deleting {file}: {e}")
            os.rmdir( dataset_name + "/" + directory)
            print(f"  Cleared {directory}")
        else:
            print(f"{directory} does not exist, skipping deletion.")

############################
## define parameters, modify as needed
dataset_name = "eQTLsummary_demo"
qtl_data_files = sorted(glob("data/eQTLsummary_sampled/*_sampled.tsv"))
geneid_col = "gene_id"
snpid_col = "variant_id"
pvalue_col = "pval_nominal"
beta_col = "slope"

############################`
## main script
print("Deleting old data...")
delete_old_data(dataset_name)

print(f"{len(qtl_data_files)} QTL data files found\n")

# Organize QTL data into celltype-specific files
print("Processing QTL data files...")
start_time = time.time()

os.makedirs(dataset_name + "/" + "unfiltered_celltypes", exist_ok=True)

for i, filepath in enumerate(qtl_data_files):
    path = Path(filepath)
    region = path.stem.split("_")[0]
    output_file = Path(dataset_name + "/" + "unfiltered_celltypes") / f"{region}.tsv"

    print(f"  [{i+1}/{len(qtl_data_files)}] Processing {region}...")
    file_start = time.time()

    file_size = os.path.getsize(filepath)
    processed_bytes = 0
    first_chunk = not output_file.exists()

    with open(filepath, "rb") as f:
        reader = pd.read_csv(f,sep="\t", usecols=[snpid_col, geneid_col, pvalue_col, beta_col],chunksize=1000000,)

        for chunk in reader:
            chunk = chunk.rename(columns={
                    geneid_col: "gene_id",
                    snpid_col: "snp_id",
                    pvalue_col: "p_value",
                    beta_col: "beta_value",
                }
            )

            chunk.to_csv(output_file,mode="a",sep="\t",index=False,header=first_chunk,float_format="%.6g")
            first_chunk = False

            current_pos = f.tell()
            processed_bytes = min(current_pos, file_size)  # Ensure we don't exceed file size
            percent = (processed_bytes / file_size) * 100

            elapsed = time.time() - file_start
            if percent > 0:  # Avoid division by zero
                eta = (elapsed / percent) * (100 - percent)
                print(f"    Progress: {percent:.1f}% | Elapsed: {elapsed:.1f}s | ETA: {eta:.1f}s", end="\r")
            else:
                print("    Progress: 0% | Starting...")

    elapsed = time.time() - file_start
    print(f"    Progress: 100.0% | Completed in {elapsed:.1f}s")

print(f"Processed all QTL files in {time.time()-start_time:.2f} seconds\n")

#############################
## filter for significant QTLs
print("Filtering for significant QTLs in each celltype...")
os.makedirs(dataset_name + "/" + "filtered_celltypes", exist_ok=True)

files = sorted(glob(dataset_name + "/" + "unfiltered_celltypes/*.tsv"))

for filepath in tqdm(files, desc="Processing files", unit="file"):
    region = Path(filepath).stem
    output_path = Path(dataset_name + "/" + "filtered_celltypes") / f"{region}.tsv"

    df = pd.read_csv(
        filepath,
        sep="\t",
        dtype={"snp_id": str, "gene_id": str, "p_value": float, "beta_value": float},
    )
    filtered_df = df[df["p_value"] < 0.01]

    filtered_df.to_csv(output_path, sep="\t", index=False, float_format="%.6g")

print(f"All done. Filtered files are in '{dataset_name}/filtered_celltypes/'.")

