import os
from pathlib import Path
import pandas as pd
import json
from glob import glob
import re
import time
from collections import defaultdict

############################################
## define parameters, modify as needed
dataset_name = "eQTLsummary_demo"
gene_annotation_file = "data/gene_annotations.tsv"
geneid_col = "gene_id"
gene_start_col = "position_start"
gene_end_col = "position_end"
gene_chrom_col = "chromosome" 
gene_strand_col = "strand"  

snp_annotation_file = "data/snp_annotations.tsv"
snpid_col = "snp_id"
snp_chrom_col = "chr"
snp_pos_col = "position"

#############################################
def safe_filename(name):
    return re.sub(r"[^a-zA-Z0-9_\-]", "_", name)


# Process gene locations
print("Processing gene locations...")
start_time = time.time()
output_dir = Path(dataset_name + "/gene_locations")
output_dir.mkdir(exist_ok=True)

chrom_data = {}

file_start = time.time()

chunks = []
df = pd.read_csv( gene_annotation_file,
    sep="\t",
    usecols=[geneid_col, gene_chrom_col, gene_start_col, gene_end_col, gene_strand_col],
    dtype={
        geneid_col: str,
        gene_chrom_col: str,
        gene_start_col: int,
        gene_end_col: int,
        gene_strand_col: str,
    },
)
df = df.rename(
    columns={
        geneid_col: "gene_id",
        gene_chrom_col: "chromosome",
        gene_start_col: "position_start",
        gene_end_col: "position_end",
        gene_strand_col: "strand",
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



#############################################
# Process SNP locations with vectorized operations and chunking
print("Processing SNP locations...")
start_time = time.time()
output_dir = Path(dataset_name + "/snp_locations")
output_dir.mkdir(exist_ok=True)

file_start = time.time()
file_snp_count = 0

seen_per_chromosome = defaultdict(set)

df = pd.read_csv(snp_annotation_file, sep="\t", usecols=[snpid_col, snp_chrom_col, snp_pos_col])
df.rename(columns={snpid_col: "snp_id", snp_chrom_col: "chr", snp_pos_col: "position"}, inplace=True)
 
for chrom, group in df.groupby("chr"):
    chrom_file = output_dir / f"{chrom}.tsv"
    group = group[["snp_id", "position"]]
    group = group[~group["snp_id"].isin(seen_per_chromosome[chrom])]
    seen_per_chromosome[chrom].update(group["snp_id"])
    if not group.empty:
        if chrom_file.exists():
            group.to_csv(chrom_file, sep="\t", index=False, mode="a", header=False)
        else:
            group.to_csv(chrom_file, sep="\t", index=False)
        file_snp_count += len(group)

print(f"  - Wrote {file_snp_count:,} SNPs in {time.time()-file_start:.2f} seconds")

print(f"Processed all SNP locations in {time.time()-start_time:.2f} seconds")


######################################### 
# Annotate genes and SNPs with celltype information

# Load celltype mapping
celltype_mapping_file = Path(dataset_name + "/celltypes/celltype_parquet.json")

celltype_mapping = json.loads(celltype_mapping_file.read_text())
filename_to_pretty = {v: k for k, v in celltype_mapping.items()}
celltype_files = sorted(glob(dataset_name + "/celltypes/*.tsv"))

# Step 1: Collect unique genes and SNPs from celltype files
print("Collecting unique genes and SNPs from celltype files...")
unique_genes = set()
unique_snps = set()

start_time = time.time()
for filepath in celltype_files:
    # Process in chunks to reduce memory usage
    for chunk in pd.read_csv(
        filepath, sep="\t", usecols=["snp_id", "gene_id"], chunksize=500000
    ):
        unique_genes.update(chunk["gene_id"].unique())
        unique_snps.update(chunk["snp_id"].unique())

print(
    f"Found {len(unique_genes):,} unique genes and {len(unique_snps):,} unique SNPs in {time.time()-start_time:.2f} seconds"
)

# Step 2: Create reverse mapping from chromosome to genes/SNPs
print("Creating chromosome index...")
start_time = time.time()

# Create mappings: chromosome -> set of genes/SNPs
chrom_to_genes = {}
chrom_to_snps = {}

# Process gene locations
for filepath in glob(dataset_name + "/gene_locations/*.tsv"):
    chrom = Path(filepath).stem
    # Only read gene_id column to save memory
    df = pd.read_csv(filepath, sep="\t", usecols=["gene_id"])
    # Filter to only genes we care about
    chrom_genes = set(df["gene_id"]).intersection(unique_genes)
    chrom_to_genes[chrom] = chrom_genes

# Process SNP locations
for filepath in glob(dataset_name + "/snp_locations/*.tsv"):
    chrom = Path(filepath).stem
    # Only read snp_id column to save memory
    df = pd.read_csv(filepath, sep="\t", usecols=["snp_id"])
    # Filter to only SNPs we care about
    chrom_snps = set(df["snp_id"]).intersection(unique_snps)
    chrom_to_snps[chrom] = chrom_snps

print(f"Created chromosome index in {time.time()-start_time:.2f} seconds")

# Step 3: Build location maps in bulk with progress tracking
print("Building location maps...")
start_time = time.time()
gene_info_map = {}
snp_info_map = {}

# Process genes in bulk with progress
total_genes = sum(len(genes) for genes in chrom_to_genes.values())
processed_genes = 0
print(f"Processing {total_genes:,} genes across {len(chrom_to_genes)} chromosomes...")

for chrom, genes in chrom_to_genes.items():
    if not genes:
        continue

    filepath = f"{dataset_name}/gene_locations/{chrom}.tsv"
    df = pd.read_csv(
        filepath,
        sep="\t",
        usecols=["gene_id", "position_start", "position_end", "strand"],
    )
    df = df[df["gene_id"].isin(genes)]

    for _, row in df.iterrows():
        gene_id = row["gene_id"]
        gene_info_map[gene_id] = {
            "chromosome": chrom,
            "position_start": int(row["position_start"]),
            "position_end": int(row["position_end"]),
            "strand": row["strand"],
        }
        processed_genes += 1
        if processed_genes % 10000 == 0:  # Print every 10k genes
            elapsed = time.time() - start_time
            print(
                f"  Genes: {processed_genes:,}/{total_genes:,} ({processed_genes/total_genes:.1%}) | "
                f"Elapsed: {elapsed:.1f}s | "
                f"ETA: {(elapsed/processed_genes)*(total_genes-processed_genes):.1f}s"
            )

# Process SNPs in bulk with progress
total_snps = sum(len(snps) for snps in chrom_to_snps.values())
processed_snps = 0
print(f"\nProcessing {total_snps:,} SNPs across {len(chrom_to_snps)} chromosomes...")

for chrom, snps in chrom_to_snps.items():
    if not snps:
        continue

    filepath = f"{dataset_name}/snp_locations/{chrom}.tsv"
    df = pd.read_csv(filepath, sep="\t", usecols=["snp_id", "position"])
    df = df[df["snp_id"].isin(snps)]

    for _, row in df.iterrows():
        snp_id = row["snp_id"]
        snp_info_map[snp_id] = {
            "chromosome": chrom,
            "position": int(row["position"]),
        }
        processed_snps += 1
        if processed_snps % 100000 == 0:  # Print every 100k SNPs
            elapsed = time.time() - start_time
            print(
                f"  SNPs: {processed_snps:,}/{total_snps:,} ({processed_snps/total_snps:.1%}) | "
                f"Elapsed: {elapsed:.1f}s | "
                f"ETA: {(elapsed/processed_snps)*(total_snps-processed_snps):.1f}s"
            )

print(
    f"\nMapped {len(gene_info_map):,} genes and {len(snp_info_map):,} SNPs in {time.time()-start_time:.2f} seconds"
)

# Step 4: Process celltype files with chunking
print("Processing celltype files...")
start_time = time.time()
gene_to_regions = {}
snp_to_regions = {}

for filepath in celltype_files:
    path = Path(filepath)
    region_filename = path.stem + ".tsv"
    region_pretty = filename_to_pretty.get(region_filename, region_filename)
    print(f"Processing {region_pretty}...")

    # Process in large chunks
    for chunk in pd.read_csv(
        filepath, sep="\t", usecols=["snp_id", "gene_id"], chunksize=1000000
    ):
        # Process genes
        for gene in chunk["gene_id"].unique():
            if gene in gene_info_map:
                gene_to_regions.setdefault(
                    gene, {**gene_info_map[gene], "celltypes": set()}
                )
                gene_to_regions[gene]["celltypes"].add(region_pretty)

        # Process SNPs
        for snp in chunk["snp_id"].unique():
            if snp in snp_info_map:
                snp_to_regions.setdefault(
                    snp, {**snp_info_map[snp], "celltypes": set()}
                )
                snp_to_regions[snp]["celltypes"].add(region_pretty)

print(f"Processed celltype files in {time.time()-start_time:.2f} seconds")

# Step 5: Write output files with progress tracking
print("Writing JSON files...")
os.makedirs(f"{dataset_name}/gene_jsons", exist_ok=True)
os.makedirs(f"{dataset_name}/snp_jsons", exist_ok=True)

# Per-gene JSONs with progress
gene_ids = list(gene_to_regions.keys())
total_genes = len(gene_ids)
start_time = time.time()
print(f"Writing {total_genes:,} gene JSON files...")

for i, (gene, data) in enumerate(gene_to_regions.items(), 1):
    gene_path = Path(f"{dataset_name}/gene_jsons") / f"{safe_filename(gene)}.json"
    with gene_path.open("w") as f:
        json.dump(
            {
                "chromosome": data["chromosome"],
                "position_start": data["position_start"],
                "position_end": data["position_end"],
                "strand": data["strand"],
                "celltypes": sorted(data["celltypes"]),
            },
            f,
        )

    if i % 1000 == 0 or i == total_genes:  # Print every 1k genes or at end
        elapsed = time.time() - start_time
        rate = i / elapsed if elapsed > 0 else 0
        print(
            f"  Genes: {i:,}/{total_genes:,} ({i/total_genes:.1%}) | "
            f"Elapsed: {elapsed:.1f}s | "
            f"Rate: {rate:.1f} files/sec | "
            f"ETA: {(total_genes - i)/rate:.1f}s"
            if rate > 0
            else "ETA: Calculating..."
        )

# Per-SNP JSONs with progress
snp_ids = list(snp_to_regions.keys())
total_snps = len(snp_ids)
start_time = time.time()
print(f"\nWriting {total_snps:,} SNP JSON files...")

for i, (snp, data) in enumerate(snp_to_regions.items(), 1):
    snp_path = Path(f"{dataset_name}/snp_jsons") / f"{safe_filename(snp)}.json"
    with snp_path.open("w") as f:
        json.dump(
            {
                "chromosome": data["chromosome"],
                "position": data["position"],
                "celltypes": sorted(data["celltypes"]),
            },
            f,
        )

    if i % 10000 == 0 or i == total_snps:  # Print every 10k SNPs or at end
        elapsed = time.time() - start_time
        rate = i / elapsed if elapsed > 0 else 0
        print(
            f"  SNPs: {i:,}/{total_snps:,} ({i/total_snps:.1%}) | "
            f"Elapsed: {elapsed:.1f}s | "
            f"Rate: {rate:.1f} files/sec | "
            f"ETA: {(total_snps - i)/rate:.1f}s"
            if rate > 0
            else "ETA: Calculating..."
        )

# Write lists
Path(f"{dataset_name}/gene_list.json").write_text(json.dumps(gene_ids, separators=(",", ":")))
Path(f"{dataset_name}/snp_list.json").write_text(json.dumps(snp_ids, separators=(",", ":")))

print("\nDone. Output files created.")
