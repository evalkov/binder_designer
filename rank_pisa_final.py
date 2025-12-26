#!/usr/bin/env python3
"""
Ranking + heatmap for CCP4 PISA contact summary.

Expected header (your file):
binder,bsa_score,salt_bridges,h_bonds,interface_count,interface_area_A2,
dg_dissociation,solvation_energy_gain,overlap,specificity,interface_residue_count,
pct_polar,pct_hydrophobic

Interpretation:
- dg_dissociation: more negative is more favorable (convert to good_dg = -dg_dissociation)
- solvation_energy_gain: more negative is more favorable (good_solv = -solvation_energy_gain)
- overlap: Yes/No, penalize strongly (or filter)
- specificity: -log10(pvalue), higher is better

Outputs:
- ranked_all_latest.csv
- ranked_no_overlap_latest.csv
- top_candidates_heatmap.png
- top_candidates_heatmap.pdf
"""

from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def zscore(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    mu = np.nanmean(s.values)
    sd = np.nanstd(s.values, ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - mu) / sd


def minmax(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    mn = np.nanmin(s.values)
    mx = np.nanmax(s.values)
    if not np.isfinite(mn) or not np.isfinite(mx) or mx == mn:
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - mn) / (mx - mn)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="contacts.csv")
    ap.add_argument("-o", "--outdir", default="rank_out", help="output directory")
    ap.add_argument("--top_n", type=int, default=20, help="top N non-overlap candidates for heatmap")
    ap.add_argument("--overlap_penalty_z", type=float, default=0.80,
                    help="penalty applied to composite score when overlap==Yes (z-score units)")
    args = ap.parse_args()

    inp = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(inp)

    # Validate minimally
    required = ["binder", "overlap"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Missing column '{c}'. Found: {list(df.columns)}")

    # Coerce numeric columns (exactly from your header)
    for c in [
        "bsa_score","salt_bridges","h_bonds","interface_count","interface_area_A2",
        "dg_dissociation","solvation_energy_gain","specificity","interface_residue_count",
        "pct_polar","pct_hydrophobic"
    ]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Overlap flag
    df["overlap_flag"] = df["overlap"].astype(str).str.strip().str.lower().isin(["yes","y","true","1"])

    # Derived metrics
    df["interactions"] = df["salt_bridges"].fillna(0) + df["h_bonds"].fillna(0)
    df["good_dg"] = -df["dg_dissociation"]                 # higher is better
    df["good_solv"] = -df["solvation_energy_gain"]         # higher is better
    df["dg_per_100A2"] = df["good_dg"] / (df["interface_area_A2"] / 100.0).replace(0, np.nan)

    # Size proxy: blend BSA + interface area
    df["size_avg"] = (minmax(df["bsa_score"]) + minmax(df["interface_area_A2"])) / 2.0

    # Composite score (z-scored, weighted)
    weights = {
        "size_avg": 0.25,
        "interactions": 0.20,
        "good_dg": 0.25,
        "good_solv": 0.10,
        "specificity": 0.15,
        "dg_per_100A2": 0.05,
    }
    wsum = sum(weights.values())
    df["score"] = sum((w/wsum) * zscore(df[k]) for k, w in weights.items())

    # Overlap penalty
    df.loc[df["overlap_flag"], "score"] -= float(args.overlap_penalty_z)

    # Rank
    ranked_all = df.sort_values("score", ascending=False).reset_index(drop=True)
    ranked_all["rank"] = np.arange(1, len(ranked_all) + 1)

    ranked_no_overlap = ranked_all[~ranked_all["overlap_flag"]].reset_index(drop=True)
    ranked_no_overlap["rank_no_overlap"] = np.arange(1, len(ranked_no_overlap) + 1)

    ranked_all.to_csv(outdir / "ranked_all_latest.csv", index=False)
    ranked_no_overlap.to_csv(outdir / "ranked_no_overlap_latest.csv", index=False)

    # Heatmap for top N non-overlap candidates
    top = ranked_no_overlap.head(int(args.top_n)).copy()

    heat_cols = [
        "bsa_score","salt_bridges","h_bonds","interactions","interface_count","interface_area_A2",
        "dg_dissociation","solvation_energy_gain","specificity",
        "interface_residue_count","pct_polar","pct_hydrophobic","dg_per_100A2","score"
    ]
    heat_cols = [c for c in heat_cols if c in top.columns]

    # Standardize each metric across the selected top set for comparability
    M = top[heat_cols].apply(pd.to_numeric, errors="coerce")
    Mz = M.apply(zscore, axis=0).fillna(0.0)

    plt.figure(figsize=(max(10, len(heat_cols)*0.7), max(6, len(top)*0.35)))
    plt.imshow(Mz.to_numpy(), aspect="auto", interpolation="nearest")
    plt.yticks(np.arange(len(top)), top["binder"].tolist(), fontsize=8)
    plt.xticks(np.arange(len(heat_cols)), heat_cols, rotation=45, ha="right", fontsize=8)
    plt.colorbar(label="z-score (within top set)")
    plt.title(f"Top {len(top)} candidates (overlap==No): standardized metrics heatmap")
    plt.tight_layout()

    plt.savefig(outdir / "top_candidates_heatmap.png", dpi=250)
    plt.savefig(outdir / "top_candidates_heatmap.pdf")
    plt.close()

    print(f"Done. Wrote outputs to: {outdir.resolve()}")
    print("Top 10 (no-overlap):")
    show = ["rank_no_overlap","binder","score","bsa_score","interface_area_A2",
            "interactions","dg_dissociation","solvation_energy_gain","specificity","overlap"]
    show = [c for c in show if c in ranked_no_overlap.columns]
    print(ranked_no_overlap.head(10)[show].to_string(index=False))


if __name__ == "__main__":
    main()

