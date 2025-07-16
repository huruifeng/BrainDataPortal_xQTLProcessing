import json
from pathlib import Path
import sys


directory = Path("filtered_celltypes").resolve()
output_file = Path("celltype_mapping.json").resolve()
mapping = {}

if not directory.is_dir():
    print(f"Directory '{directory}' does not exist.")
    sys.exit(1)

files = [f for f in directory.iterdir() if f.is_file() and f.suffix == ".tsv"]

if not files:
    print(f"No files found in '{directory}'")
    sys.exit(1)

print(
    "Please provide display names (e.g., Astrocytes, GABAergic_Neurons, etc.) for the following files:"
)
for filename in files:
    stem = filename.stem
    print(f"File: {filename.name}")
    pretty_name = (
        input(f"Enter display name (Blank for `{stem}`): ").strip().replace(" ", "_")
    )  # Don't allow spaces
    if not pretty_name:
        pretty_name = stem.replace("_", " ").title().replace(" ", "_")
    mapping[pretty_name] = filename.name.replace(".tsv", ".parquet")

output_file.parent.mkdir(parents=True, exist_ok=True)
output_file.write_text(json.dumps(mapping, indent=2))

print(f"\nMapping saved to {output_file.name}")
