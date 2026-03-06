from __future__ import annotations
import os
import warnings
import csv
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from sklearn.cluster import AffinityPropagation


from .core import (
    compute_distance_matrix,
    group_fragments,
    apply_aggressive_grouping,
    replace_cut_sites_and_pad,
)

# ----------------------------------------------------------------------
# Defines the lengths of the fragments 
MIN_LENGTH = 301
MAX_LENGTH = 499
# ----------------------------------------------------------------------

def DNA_clustering(
    csv_path: os.PathLike,
    *,
    aggressive: bool,
    progress: Callable[[float], None] | None = None,
) -> Path:
    if progress:
        progress(0)
        
    # Determine if .csv uses , or ; as separator
    with open(csv_path, newline='') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read())
    

    if dialect.delimiter == ',':
        df = pd.read_csv(csv_path)            # Import the csv with a comma as the separator
    elif dialect.delimiter == ';':
        df = pd.read_csv(csv_path, sep=';')   # Import the csv with a semicolon as the separator

    if progress:
        progress(2)

    df = df.apply(lambda x: x.astype(str).str.upper())
    df_large = df[df["Sequence"].apply(len) >= 500]
    # filter short fragments and reset index so numpy arrays align
    df = df[df["Sequence"].apply(len) < MAX_LENGTH].reset_index(drop=True)

    # distance matrix + clustering
    fragments = df["Sequence"].to_numpy()
    sim = compute_distance_matrix(fragments, progress)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        affprop = AffinityPropagation(affinity="precomputed", random_state=0)
        affprop.fit(sim)

    clusters = np.empty(len(df), dtype=int)
    for cid in np.unique(affprop.labels_):
        clusters[affprop.labels_ == cid] = cid
    df["Cluster"] = clusters
    if progress:
        progress(90)

    df1 = group_fragments(df, progress=progress)
    if aggressive:
        df1 = apply_aggressive_grouping(df1)

    grouped = (
        df1.groupby("Group").agg({"Sequence": "".join, "Name": list}).reset_index()
    )
    grouped["Length"] = grouped["Sequence"].str.len()

    # add back large fragments
    max_group = grouped["Group"].max() if not grouped.empty else 0
    for _, row in df_large.iterrows():
        max_group += 1
        grouped = pd.concat(
            [
                grouped,
                pd.DataFrame(
                    {
                        "Group": [max_group],
                        "Sequence": row["Sequence"],
                        "Name": row["Name"],
                        "Length": [len(row["Sequence"])],
                    }
                ),
            ],
            ignore_index=True,
        )

    grouped = replace_cut_sites_and_pad(grouped)
    if progress:
        progress(98)

    if aggressive == True:
        out_path = Path(csv_path).with_name(Path(csv_path).stem + "_aggressive_grouped.csv")
    else:
        out_path = Path(csv_path).with_name(Path(csv_path).stem + "_grouped.csv")
    
    # Output file will have either ; or , depending on the input
    if dialect.delimiter == ',':
        grouped.to_csv(out_path, index=False, sep=",")
    elif dialect.delimiter == ';':
        grouped.to_csv(out_path, index=False, sep=";")

    if progress:
        progress(100)
    return out_path



# Step 0: Defining the overhangs and types

start_constant = "gcatcgtctcatcggtctca"
end_constant = "tgagacctgagacggcat"

overhangs = pd.DataFrame(
    np.array([
    ["1", "CCCT", "AACG"],
    ["2", "AACG", "TATG"],
    ["2a", "AACG", "AAAC"],
    ["2b", "AAAC", "CTGA"],
    ["2ab", "AACG", "CTGA"],
    ["2c", "CTGA", "AAGA"],
    ["2d", "AAGA", "TATG"],
    ["2cd", "CTGA", "TATG"],
    ["2cd34", "CTGA", "GCTG"],
    ["3", "T", "GGATCC"],         # The YTK framework uses the ATG start-codon as part of the overhang. It is therefore not extra added here
    ["3a", "T", "GGTTCT"],        # The YTK framework uses the ATG start-codon as part of the overhang. It is therefore not extra added here
    ["3b", "TTCT", "GGATCC"],     # An extra GG is added before the 3' site, per the YTK framework
    ["4", "ATCC", "GCTG"],
    ["4a", "ATCC", "TGGC"],
    ["4b", "TGGC", "GCTG"],
    ["5", "GCTG", "TACA"],
    ["6", "TACA", "GAGT"],
    ["7", "GAGT", "CCGA"],
    ["8", "CCGA", "CCCT"],
    ["8a", "CCGA", "CAAT"],
    ["8b", "CAAT", "CCCT"],
    ["678", "TACA", "CCCT"]]),
    columns=["Type", "5' overhang", "3' overhang"])



def DNA_typer(csv_path):
    
    # Step 1: Load the .csv file
        
    # Determine if .csv uses , or ; as separator
    with open(csv_path, newline='') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read())
    
    if dialect.delimiter == ',':
        df = pd.read_csv(csv_path)            # Import the csv with a comma as the separator
    elif dialect.delimiter == ';':
        df = pd.read_csv(csv_path, sep=';')   # Import the csv with a semicolon as the separator
    
    
    # Step 2: Adding the new overhangs
    for index, row in df.iterrows():
        # We start by extrancting the overhangs for the type as strings
        five_overhang  = overhangs[overhangs['Type'] == row['Type']]["5' overhang"].astype(str).values.flatten()[0]
        three_overhang = overhangs[overhangs['Type'] == row['Type']]["3' overhang"].astype(str).values.flatten()[0]
    
        # The new sequence can then be constructed before it is saved
        new_sequence = start_constant + five_overhang + row['Sequence'] + three_overhang + end_constant
        df.loc[index, "New_seq"] = new_sequence
    
    # We can now rename the old sequence to "Part_sequence" and "New_seq" to "Sequence", as we will need that later
    typed_df = df.rename({'Sequence': 'Part_sequence', 'New_seq': 'Sequence'}, axis='columns')
    
    
    # Step 3: Save the outputs
    out_path = Path(csv_path).with_name(Path(csv_path).stem + "_typed.csv")
        
    # Output file will have either ; or , depending on the input
    if dialect.delimiter == ',':
        typed_df.to_csv(out_path, index=False, sep=",")
    elif dialect.delimiter == ';':
        typed_df.to_csv(out_path, index=False, sep=";")

    return out_path