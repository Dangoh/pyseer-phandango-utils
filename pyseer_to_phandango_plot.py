#!/usr/bin/env python3
"""
pyseer_to_phandango_plot.py

Convert Pyseer output (e.g., SNP GWAS from --lmm or default linear model) into a
Phandango-compatible '.plot' track file.

Phandango GWAS/track format expected (tab-delimited):
#CHR    SNP     BP      minLOG10(P)  log10(p)    r^2

This script:
  - Reads Pyseer output (tab-delimited, header required).
  - Extracts genomic position (BP) from the 'variant' column by splitting on a delimiter
    (default: underscore) and taking a specified field index (default: 2nd field).
    Example variant: "AE017143.1_12345_A_T" -> BP = 12345 (field index 1).
  - Uses a p-value column (auto-detect or specified) to compute -log10(p).
  - Writes a .plot file suitable for Phandango.

Typical usage:
  python pyseer_to_phandango_plot.py --pyseer infection_SNPs_lmm.txt --out infection_SNPs_lmm.plot

Author: (you)
"""

import argparse
import math
import sys
from typing import Optional, List


DEFAULT_HEADER = "#CHR\tSNP\tBP\tminLOG10(P)\tlog10(p)\tr^2"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert Pyseer output to Phandango .plot track."
    )
    p.add_argument("--pyseer", required=True, help="Pyseer output TSV (with header).")
    p.add_argument("--out", required=True, help="Output .plot file.")
    p.add_argument(
        "--chr-index",
        default="26",
        help="Value to put in #CHR column (default: 26).",
    )
    p.add_argument(
        "--variant-col",
        default="variant",
        help="Column name containing variant IDs (default: variant).",
    )
    p.add_argument(
        "--pcol",
        default=None,
        help=(
            "P-value column name to use. If not set, will auto-detect from common "
            "Pyseer columns (lrt-pvalue, pvalue, p-value, wald-pvalue, etc.)."
        ),
    )
    p.add_argument(
        "--variant-delim",
        default="_",
        help="Delimiter used in variant IDs (default: '_').",
    )
    p.add_argument(
        "--bp-field-index",
        type=int,
        default=1,
        help=(
            "0-based index of the field in the split variant ID that contains BP. "
            "Default 1 (second field). Example: contig_12345_A_T => BP field is 1."
        ),
    )
    p.add_argument(
        "--snp-name",
        default=".",
        help="Value to put in SNP column (default: '.').",
    )
    p.add_argument(
        "--r2",
        default="0",
        help="Value to put in r^2 column (default: 0).",
    )
    p.add_argument(
        "--skip-nonpositive-p",
        action="store_true",
        help="Skip rows with p<=0 instead of erroring.",
    )
    p.add_argument(
        "--allow-missing-bp",
        action="store_true",
        help="Skip rows where BP cannot be parsed (default: error).",
    )
    return p.parse_args()


def detect_pcol(header: List[str]) -> Optional[str]:
    """Try to detect the most likely p-value column."""
    candidates = [
        "lrt-pvalue",
        "pvalue",
        "p-value",
        "pval",
        "wald-pvalue",
        "score-pvalue",
        "filter-pvalue",
    ]
    lower = {h.lower(): h for h in header}
    for c in candidates:
        if c in lower:
            return lower[c]
    # Heuristic: any column containing 'pvalue'
    for h in header:
        if "pvalue" in h.lower() or h.lower() in ("p", "pval"):
            return h
    return None


def safe_neglog10(p: float) -> float:
    # Guard against underflow; Phandango can handle large values,
    # but p can be extremely small. If p==0, -log10 is inf; we prevent that.
    if p <= 0.0:
        return float("inf")
    return -math.log10(p)


def main() -> None:
    args = parse_args()

    # Read input
    with open(args.pyseer, "r", encoding="utf-8") as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]

    if not lines:
        raise SystemExit(f"[ERROR] Input file is empty: {args.pyseer}")

    header = lines[0].split("\t")
    col_index = {name: i for i, name in enumerate(header)}

    if args.variant_col not in col_index:
        raise SystemExit(
            f"[ERROR] Variant column '{args.variant_col}' not found. "
            f"Available columns: {header}"
        )

    pcol = args.pcol
    if pcol is None:
        pcol = detect_pcol(header)

    if pcol is None or pcol not in col_index:
        raise SystemExit(
            "[ERROR] Could not determine p-value column. "
            "Specify explicitly with --pcol. "
            f"Available columns: {header}"
        )

    var_i = col_index[args.variant_col]
    p_i = col_index[pcol]

    out_rows = []
    out_rows.append(DEFAULT_HEADER)

    bad_bp = 0
    bad_p = 0

    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) <= max(var_i, p_i):
            continue

        variant = parts[var_i]
        p_str = parts[p_i]

        # Parse BP from variant ID
        fields = variant.split(args.variant_delim)
        try:
            bp = fields[args.bp_field_index]
        except IndexError:
            bad_bp += 1
            if args.allow_missing_bp:
                continue
            raise SystemExit(
                f"[ERROR] Could not parse BP: variant '{variant}' does not have "
                f"field index {args.bp_field_index} when split by '{args.variant_delim}'."
            )

        # Validate BP numeric
        try:
            bp_int = int(bp)
        except ValueError:
            bad_bp += 1
            if args.allow_missing_bp:
                continue
            raise SystemExit(
                f"[ERROR] BP field '{bp}' parsed from variant '{variant}' is not an integer."
            )

        # Parse p-value
        try:
            p_val = float(p_str)
        except ValueError:
            bad_p += 1
            if args.skip_nonpositive_p:
                continue
            raise SystemExit(
                f"[ERROR] P-value '{p_str}' in column '{pcol}' is not numeric (variant={variant})."
            )

        if p_val <= 0.0:
            bad_p += 1
            if args.skip_nonpositive_p:
                continue
            # Pyseer can output extremely small p-values but not negative; p==0 can happen due to rounding.
            # Represent as a large value rather than failing.
            neglog = float("inf")
        else:
            neglog = safe_neglog10(p_val)

        # Phandango expects both minLOG10(P) and log10(p); many pipelines duplicate the same value
        # to satisfy schema.
        out_rows.append(
            "\t".join(
                [
                    str(args.chr_index),
                    str(args.snp_name),
                    str(bp_int),
                    str(neglog),
                    str(neglog),
                    str(args.r2),
                ]
            )
        )

    with open(args.out, "w", encoding="utf-8") as out:
        out.write("\n".join(out_rows) + "\n")

    sys.stderr.write(
        f"[INFO] Wrote {len(out_rows)-1} variants to {args.out} "
        f"(pcol='{pcol}', variant_col='{args.variant_col}').\n"
    )
    if bad_bp:
        sys.stderr.write(f"[WARN] Skipped/failed BP parsing for {bad_bp} rows.\n")
    if bad_p:
        sys.stderr.write(f"[WARN] Encountered non-positive/non-numeric p for {bad_p} rows.\n")


if __name__ == "__main__":
    main()

