import os
import pandas as pd
from pathlib import Path
from glob import glob
import time
import re


def safe_filename(name):
    return re.sub(r"[^a-zA-Z0-9_\-]", "_", name)


def delete_old_data():
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
        if os.path.exists(directory):
            print(f"Deleting old data in {directory}...")
            for file in glob(os.path.join(directory, "*")):
                try:
                    os.remove(file)
                except Exception as e:
                    print(f"  Error deleting {file}: {e}")
            os.rmdir(directory)
            print(f"  Cleared {directory}")
        else:
            print(f"{directory} does not exist, skipping deletion.")


print("Deleting old data...")
delete_old_data()


# Organize QTL data into celltype-specific files
print("Processing QTL data files...")
start_time = time.time()
qtl_data_files = sorted(glob("eQTLsummary_unfiltered/*_combined_nominals.tsv"))
os.makedirs("unfiltered_celltypes", exist_ok=True)

for i, filepath in enumerate(qtl_data_files):
    path = Path(filepath)
    region = path.stem.split("_")[0]
    output_file = Path("unfiltered_celltypes") / f"{region}.tsv"

    print(f"  [{i+1}/{len(qtl_data_files)}] Processing {region}...")
    file_start = time.time()

    file_size = os.path.getsize(filepath)
    processed_bytes = 0
    first_chunk = not output_file.exists()

    with open(filepath, "rb") as f:
        reader = pd.read_csv(
            f,
            sep="\t",
            usecols=["variant_id", "gene_id", "pval_nominal", "slope"],
            chunksize=1000000,
        )

        for chunk in reader:
            chunk = chunk.rename(
                columns={
                    "gene_id": "gene_id",
                    "variant_id": "snp_id",
                    "pval_nominal": "p_value",
                    "slope": "beta_value",
                }
            )

            chunk.to_csv(
                output_file,
                mode="a",
                sep="\t",
                index=False,
                header=first_chunk,
                float_format="%.6g",
            )
            first_chunk = False

            current_pos = f.tell()
            processed_bytes = min(
                current_pos, file_size
            )  # Ensure we don't exceed file size
            percent = (processed_bytes / file_size) * 100

            elapsed = time.time() - file_start
            if percent > 0:  # Avoid division by zero
                eta = (elapsed / percent) * (100 - percent)
                print(
                    f"    Progress: {percent:.1f}% | "
                    f"Elapsed: {elapsed:.1f}s | "
                    f"ETA: {eta:.1f}s"
                )
            else:
                print("    Progress: 0% | Starting...")

    elapsed = time.time() - file_start
    print(f"    Progress: 100.0% | Completed in {elapsed:.1f}s")

print(f"Processed all QTL files in {time.time()-start_time:.2f} seconds\n")

# Process gene locations
print("Processing gene locations...")
start_time = time.time()
gene_location_files = sorted(glob("geneLoc/*.gene.hg38.bed"))
output_dir = Path("gene_locations")
output_dir.mkdir(exist_ok=True)

chrom_data = {}
total_files = len(gene_location_files)

for i, filepath in enumerate(gene_location_files):
    print(f"  [{i+1}/{total_files}] Reading {Path(filepath).name}...")
    file_start = time.time()

    chunks = []
    for chunk in pd.read_csv(
        filepath,
        sep="\t",
        usecols=["gene_name", "chr", "start", "end", "strand"],
        dtype={"gene_name": str, "chr": str, "start": int, "end": int, "strand": str},
        chunksize=100000,
    ):
        chunks.append(chunk)

    df = pd.concat(chunks)
    df = df.rename(
        columns={
            "gene_name": "gene_id",
            "chr": "chromosome",
            "start": "position_start",
            "end": "position_end",
            "strand": "strand",
        }
    )

    # Group by chromosome
    for chrom, group in df.groupby("chromosome"):
        if chrom not in chrom_data:
            chrom_data[chrom] = []
        chrom_data[chrom].append(group)

    print(f"  - Processed in {time.time()-file_start:.2f} seconds")

# Write chromosome files with deduplication
print("\nWriting gene location files...")
for chrom, groups in chrom_data.items():
    df = pd.concat(groups)
    df.drop_duplicates(
        subset=["gene_id", "position_start", "position_end", "strand"], inplace=True
    )
    df = df.drop(columns=["chromosome"])
    df.to_csv(output_dir / f"{chrom}.tsv", sep="\t", index=False)
    print(f"  - Chromosome {chrom}: {len(df):,} genes")

print(f"Processed gene locations in {time.time()-start_time:.2f} seconds\n")

# Process SNP locations with vectorized operations and chunking
print("Processing SNP locations...")
start_time = time.time()
snp_location_files = sorted(glob("snpLoc/*.txt"))
output_dir = Path("snp_locations")
output_dir.mkdir(exist_ok=True)
total_files = len(snp_location_files)

for i, filepath in enumerate(snp_location_files):
    chromosome = Path(filepath).stem.split("-")[0]
    print(
        f"  [{i+1}/{total_files}] Processing {Path(filepath).name} as {chromosome}..."
    )
    file_start = time.time()

    def process_snp_chunk(chunk):
        chunk["snp_id"] = chunk["snp"].str.split("_").str[-1]
        return chunk.drop(columns=["snp", "chr"]).rename(columns={"pos": "position"})

    chunks = []
    for chunk in pd.read_csv(
        filepath, sep="\t", usecols=["snp", "chr", "pos"], chunksize=500000
    ):
        chunks.append(process_snp_chunk(chunk))

    df = pd.concat(chunks)
    df.drop_duplicates(subset=["snp_id", "position"], inplace=True)
    df = df[["snp_id", "position"]]

    df.to_csv(output_dir / f"{chromosome}.tsv", sep="\t", index=False)
    print(f"  - Wrote {len(df):,} SNPs in {time.time()-file_start:.2f} seconds")

print(f"Processed all SNP locations in {time.time()-start_time:.2f} seconds")
