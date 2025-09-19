import pandas as pd
import os

dataset_name = "eQTLsummary_demo"
input_file = "PD_GWAS_with_rsIDs_hg38.tab"
gwas_chrom_col = "#chrom"
gwas_pos_col = "chromEnd"
gwas_snp_col = "rsID"
gwas_beta_col = "b"
gwas_pval_col = "p"


output_dir = f"{dataset_name}/gwas"
process_gwas_data(input_file, output_dir)

os.makedirs(output_dir, exist_ok=True)

df = pd.read_csv(input_file, sep="\t")

df = df.rename(
    columns={
        gwas_chrom_col: "chromosome",
        gwas_pos_col: "position",
        gwas_beta_col: "beta_value",
        gwas_pval_col: "p_value",
        gwas_snp_col: "snp_id",
    }
)

output_cols = ["chromosome", "snp_id", "position", "beta_value", "p_value"]
df = df[output_cols]

for chrom, group in df.groupby("chromosome"):
    formatted_group = group.copy()
    formatted_group["beta_value"] = formatted_group["beta_value"].map(
        "{0:.6g}".format
    )
    formatted_group["p_value"] = formatted_group["p_value"].map("{0:.6g}".format)
    formatted_group = formatted_group.drop(columns=["chromosome"])

    output_file = os.path.join(output_dir, f"{chrom}.tsv")

    formatted_group.to_csv(output_file, sep="\t", index=False, float_format="%0.6g")
    print(f"Saved {len(formatted_group)} variants to {output_file}")

