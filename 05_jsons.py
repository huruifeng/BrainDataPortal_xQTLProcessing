import os
from pathlib import Path
import pandas as pd
import json
from glob import glob
import re
import time


def safe_filename(name):
    return re.sub(r"[^a-zA-Z0-9_\-]", "_", name)


celltype_mapping_file = Path("celltype_mapping.json")
celltype_mapping = json.loads(celltype_mapping_file.read_text())
filename_to_pretty = {
    v.replace(".parquet", ".tsv"): k for k, v in celltype_mapping.items()
}
files = sorted(glob("celltypes/*.tsv"))

print("Collecting unique genes and SNPs from celltype files...")
unique_genes = set()
unique_snps = set()

start_time = time.time()
for filepath in files:
    print(filepath)
    for chunk in pd.read_csv(
        filepath, sep="\t", usecols=["snp_id", "gene_id"], chunksize=500000
    ):
        unique_genes.update(chunk["gene_id"].unique())
        unique_snps.update(chunk["snp_id"].unique())

print(
    f"Found {len(unique_genes):,} unique genes and {len(unique_snps):,} unique SNPs in {time.time()-start_time:.2f} seconds"
)

print("Creating chromosome index...")
start_time = time.time()

chrom_to_genes = {}
chrom_to_snps = {}

for filepath in glob("gene_locations/*.tsv"):
    chrom = Path(filepath).stem
    df = pd.read_csv(filepath, sep="\t", usecols=["gene_id"])
    chrom_genes = set(df["gene_id"]).intersection(unique_genes)
    chrom_to_genes[chrom] = chrom_genes

for filepath in glob("snp_locations/*.tsv"):
    chrom = Path(filepath).stem
    df = pd.read_csv(filepath, sep="\t", usecols=["snp_id"])
    chrom_snps = set(df["snp_id"]).intersection(unique_snps)
    chrom_to_snps[chrom] = chrom_snps

print(f"Created chromosome index in {time.time()-start_time:.2f} seconds")

print("Building location maps...")
start_time = time.time()
gene_info_map = {}
snp_info_map = {}

total_genes = sum(len(genes) for genes in chrom_to_genes.values())
processed_genes = 0
print(f"Processing {total_genes:,} genes across {len(chrom_to_genes)} chromosomes...")

for chrom, genes in chrom_to_genes.items():
    if not genes:
        continue

    filepath = f"gene_locations/{chrom}.tsv"
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

total_snps = sum(len(snps) for snps in chrom_to_snps.values())
processed_snps = 0
print(f"\nProcessing {total_snps:,} SNPs across {len(chrom_to_snps)} chromosomes...")

for chrom, snps in chrom_to_snps.items():
    if not snps:
        continue

    filepath = f"snp_locations/{chrom}.tsv"
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

print("Processing celltype files...")
start_time = time.time()
gene_to_regions = {}
snp_to_regions = {}

for filepath in files:
    path = Path(filepath)
    region_filename = path.stem + ".tsv"
    region_pretty = filename_to_pretty.get(region_filename, region_filename)
    print(f"Processing {region_pretty}...")

    for chunk in pd.read_csv(
        filepath, sep="\t", usecols=["snp_id", "gene_id"], chunksize=1000000
    ):
        for gene in chunk["gene_id"].unique():
            if gene in gene_info_map:
                gene_to_regions.setdefault(
                    gene, {**gene_info_map[gene], "celltypes": set()}
                )
                gene_to_regions[gene]["celltypes"].add(region_pretty)

        for snp in chunk["snp_id"].unique():
            if snp in snp_info_map:
                snp_to_regions.setdefault(
                    snp, {**snp_info_map[snp], "celltypes": set()}
                )
                snp_to_regions[snp]["celltypes"].add(region_pretty)

print(f"Processed celltype files in {time.time()-start_time:.2f} seconds")

print("Writing JSON files...")
os.makedirs("gene_jsons", exist_ok=True)
os.makedirs("snp_jsons", exist_ok=True)

gene_ids = list(gene_to_regions.keys())
total_genes = len(gene_ids)
start_time = time.time()
print(f"Writing {total_genes:,} gene JSON files...")

for i, (gene, data) in enumerate(gene_to_regions.items(), 1):
    gene_path = Path("gene_jsons") / f"{safe_filename(gene)}.json"
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

snp_ids = list(snp_to_regions.keys())
total_snps = len(snp_ids)
start_time = time.time()
print(f"\nWriting {total_snps:,} SNP JSON files...")

for i, (snp, data) in enumerate(snp_to_regions.items(), 1):
    snp_path = Path("snp_jsons") / f"{safe_filename(snp)}.json"
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

Path("gene_list.json").write_text(json.dumps(gene_ids, separators=(",", ":")))
Path("snp_list.json").write_text(json.dumps(snp_ids, separators=(",", ":")))

print("\nDone. Output files created.")
