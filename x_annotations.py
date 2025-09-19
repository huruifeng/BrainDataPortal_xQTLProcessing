## sampling the data records from the dataset
import os
import time
import pandas as pd
from glob import glob
from pathlib import Path

# Process gene locations
print("Processing gene locations...")
start_time = time.time()
gene_location_files = sorted(glob("data/mingming/geneLoc/*.gene.hg38.bed"))
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

## merge all chromosome files into one file with chromosome column
all_genes = []
for filepath in sorted(output_dir.glob("*.tsv")):
    df = pd.read_csv(filepath, sep="\t")
    df["chromosome"] = Path(filepath).stem
    all_genes.append(df)

df = pd.concat(all_genes)
df.to_csv(output_dir / "gene_annotations.tsv", sep="\t", index=False)
print(f"  - Wrote {len(df):,} genes in {time.time()-start_time:.2f} seconds")

print(f"Processed gene locations in {time.time()-start_time:.2f} seconds\n")

# Process SNP locations with vectorized operations and chunking
print("Processing SNP locations...")
start_time = time.time()
snp_location_files = sorted(glob("data/mingming/snpLoc/*.txt"))
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


## merge all chromosome files into one file
all_snps = []
for filepath in sorted(output_dir.glob("*.tsv")):
    df = pd.read_csv(filepath, sep="\t")
    df["chr"] = Path(filepath).stem
    all_snps.append(df)

df = pd.concat(all_snps)
df.drop_duplicates(subset=["snp_id", "chr", "position"], inplace=True)
df.to_csv(output_dir / "snp_annotations.tsv", sep="\t", index=False)
print(f"  - Wrote {len(df):,} SNPs in {time.time()-start_time:.2f} seconds")
