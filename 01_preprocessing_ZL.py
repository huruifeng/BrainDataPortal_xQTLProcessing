import os
import pandas as pd
from pathlib import Path
from glob import glob
import time
import re
from collections import defaultdict


def safe_filename(name):
    return re.sub(r"[^a-zA-Z0-9_\-]", "_", name)


# Organize QTL data into celltype-specific files
print("Processing QTL data files...")
start_time = time.time()
qtl_data_files = sorted(glob("eQTLsummary/eQTL.*.with_allele_info.tsv"))
os.makedirs("unfiltered_celltypes", exist_ok=True)

for i, filepath in enumerate(qtl_data_files):
    path = Path(filepath)
    region = path.stem.split(".")[1]
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
            usecols=["SNP", "gene", "p-value", "beta"],
            chunksize=1000000,
        )

        for chunk in reader:
            chunk = chunk.rename(
                columns={
                    "SNP": "snp_id",
                    "gene": "gene_id",
                    "p-value": "p_value",
                    "beta": "beta_value",
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
gene_location_files = ["annotations/gene_positions_GRC38.tsv"]
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
        usecols=["gene_id", "chromosome", "position_start", "position_end", "strand"],
        dtype={
            "chromosome": str,
            "position_start": int,
            "position_end": int,
            "strand": str,
        },
        chunksize=100000,
    ):
        chunks.append(chunk)

    df = pd.concat(chunks)
    df = df.rename(
        columns={
            "gene_id": "gene_id",
            "chromosome": "chromosome",
            "position_start": "position_start",
            "position_end": "position_end",
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
snp_location_files = [
    "annotations/allele_align_to_PD_GWAS_PD5D.Merged.bulk_tissue.PD5D_MTG.tsv"
]
output_dir = Path("snp_locations")
output_dir.mkdir(exist_ok=True)
total_files = len(snp_location_files)

for i, filepath in enumerate(snp_location_files):
    print(f"  [{i+1}/{total_files}] Processing {Path(filepath).name}...")
    file_start = time.time()
    file_snp_count = 0

    seen_per_chromosome = defaultdict(set)

    chunk_reader = pd.read_csv(
        filepath, sep="\t", usecols=["rsID", "chr", "pos"], chunksize=500000
    )

    for chunk in chunk_reader:
        chunk = chunk.rename(columns={"rsID": "snp_id", "pos": "position"})

        for chrom, group in chunk.groupby("chr"):
            if group.empty:
                continue

            snp_tuples = list(zip(group.snp_id, group.position))

            mask = []
            new_tuples = []
            for t in snp_tuples:
                if t not in seen_per_chromosome[chrom]:
                    mask.append(True)
                    new_tuples.append(t)
                else:
                    mask.append(False)

            seen_per_chromosome[chrom].update(new_tuples)

            new_snps = group.loc[mask, ["snp_id", "position"]]
            new_count = len(new_snps)
            file_snp_count += new_count

            if new_count > 0:
                chrom_file = output_dir / f"chr{chrom}.tsv"
                write_header = not chrom_file.exists()
                new_snps.to_csv(
                    chrom_file, sep="\t", mode="a", header=write_header, index=False
                )

    print(f"  - Wrote {file_snp_count:,} SNPs in {time.time()-file_start:.2f} seconds")

print(f"Processed all SNP locations in {time.time()-start_time:.2f} seconds")
