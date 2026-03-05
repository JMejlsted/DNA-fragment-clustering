from __future__ import annotations

import numpy as np
import pandas as pd
import random
import distance
from typing import Callable, List, Set

# ----------------------------------------------------------------------
# Defines the cut-sites and lengths of the fragments accordingly

BsmBI      = "CGTCTC";  BsmBI_rev  = "GAGACG"
BsaI       = "GGTCTC";  BsaI_rev   = "GAGACC"
BbsI       = "GAAGACT"; BbsI_rev   = "AGTCTTC"
BspMI      = "ACCTGCTTA"; BspMI_rev = "TAAGCAGGT"

NUCL       = ["A", "T", "C", "G"]
MIN_LENGTH = 301
MAX_LENGTH = 499


# ----------------------------------------------------------------------

def nth_repl(s: str, old: str, new: str, n: int) -> str:
    """Replace the nâ€‘th occurrence of *old* in *s* by *new*."""
    find = s.find(old)
    i = int(find != -1)
    while find != -1 and i != n:
        find = s.find(old, find + 1)
        i += 1
    if i == n and i <= len(s.split(old)) - 1:
        return s[:find] + new + s[find + len(old) :]
    return s

# ----------------------------------------------------------------------
# Pipeline functions

def compute_distance_matrix(
    fragments: np.ndarray,
    progress: Callable[[float], None] | None = None,
) -> np.ndarray:
    """Compute Levenshtein similarity matrix, emitting â‰¤100 progress updates."""
    n = len(fragments)
    matrix = np.empty((n, n), dtype=float)
    step = max(1, n // 100)
    for i, f2 in enumerate(fragments):
        if i % step == 0 and progress:
            progress((i / n) * 80 + 2)
        matrix[i] = [-distance.levenshtein(f1, f2) for f1 in fragments]
    if progress:
        progress(90)
    return matrix


def group_fragments(
    df: pd.DataFrame,
    *,
    max_length: int = MAX_LENGTH,
    progress: Callable[[float], None] | None = None,
) -> pd.DataFrame:
    fragments = df["Sequence"].tolist()
    used_fragments: Set[str] = set()
    group_labels = np.empty(len(fragments))
    group_counter = 0
    iter_counter = 0

    while len(used_fragments) < len(fragments):
        group: List[str] = []
        total_length = 0
        clusters: Set[int] = set()

        for idx, row in df.iterrows():
            frag = row["Sequence"]
            cluster = row["Cluster"]
            if frag in used_fragments:
                continue
            frag_length = len(frag)
            if (
                len(group) < 3 and total_length + frag_length <= max_length and cluster not in clusters
            ):
                group.append(frag)
                clusters.add(cluster)
                total_length += frag_length
                used_fragments.add(frag)
                group_labels[idx] = group_counter
            elif len(group) == 0 and frag_length >= 500:
                group_labels[idx] = group_counter
                group_counter += 1
                used_fragments.add(frag)
        group_counter += 1
        iter_counter += 1
        if progress:
            progress(90 + (len(used_fragments) / len(fragments)) * 5)  # 90â€“95 %
        if iter_counter > len(fragments) * 3:
            break
    out = df.copy()
    out["Group"] = group_labels.astype(int)
    return out


def apply_aggressive_grouping(df1: pd.DataFrame) -> pd.DataFrame:
    df = df1.copy()
    singleton_groups = [g for g, sub in df.groupby("Group") if len(sub) == 1]
    for i, gid in enumerate(singleton_groups):
        if i % 3 == 0:
            new = gid
        else:
            df.loc[df["Group"] == gid, "Group"] = new
    return df


def replace_cut_sites_and_pad(grouped_df: pd.DataFrame) -> pd.DataFrame:
    df = grouped_df.copy()
    for idx, row in df.iterrows():
        group, sequence, name, length = row
        new_seq = nth_repl(sequence, BsmBI, BbsI, 2)
        new_seq = nth_repl(new_seq, BsmBI_rev, BbsI_rev, 2)
        new_seq = nth_repl(new_seq, BsmBI, BspMI, 2)
        new_seq = nth_repl(new_seq, BsmBI_rev, BspMI_rev, 2)
        if len(new_seq) < 300:
            padding = "".join(random.choices(NUCL, k=MIN_LENGTH - len(new_seq)))
            for cut in (
                BsmBI, BsaI, BbsI, BspMI,
                BsmBI_rev, BsaI_rev, BbsI_rev, BspMI_rev,
            ):
                padding = padding.replace(cut, "ATCCGATGGTC")
            new_seq += padding
        df.at[idx, "Sequence"] = new_seq
        df.at[idx, "Length"] = len(new_seq)
    return df
