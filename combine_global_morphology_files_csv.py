# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 20:24:38 2025

@author: Emily Gemmill
"""
import os
import pandas as pd

# Set your path
csv_dir = "PATH/TO/DIR"  # rename to main dir
output_dir = os.path.join(csv_dir, "compiled_properties")
os.makedirs(output_dir, exist_ok=True)

property_tables = {}  # {property: {genotype: list of values}}

# Iterate over each genotype CSV
for fname in os.listdir(csv_dir):
    if not fname.endswith(".csv"):
        continue

    genotype = fname.replace("_morphology_summary.csv", "")
    fpath = os.path.join(csv_dir, fname)

    df = pd.read_csv(fpath, index_col=0)

    for prop in df.index:
        if prop not in property_tables:
            property_tables[prop] = {}

        values = df.loc[prop].dropna().values.tolist()
        property_tables[prop][genotype] = values

# Save one CSV per morphological property
for prop, geno_dict in property_tables.items():
    # Convert to DataFrame by padding shorter columns with NaNs
    max_len = max(len(v) for v in geno_dict.values())
    padded = {g: v + [None] * (max_len - len(v)) for g, v in geno_dict.items()}

    stacked_df = pd.DataFrame(padded)
    stacked_df = stacked_df[sorted(stacked_df.columns)]  # optional sort

    out_path = os.path.join(output_dir, f"{prop}.csv")
    stacked_df.to_csv(out_path, index=False)

    print(f"Saved: {out_path}")
