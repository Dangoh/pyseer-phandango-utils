# Pyseer → Phandango utilities

Utilities to convert **Pyseer GWAS output** into **Phandango-compatible `.plot` tracks**
for genome-wide visualisation of association signals along a reference genome.

This repository contains lightweight **post-processing scripts** used to annotate and
visualise bacterial GWAS results. The scripts **do not perform statistical testing**
and **do not modify association results**, but reformat outputs for interpretation
and plotting.

---

## Background

Pyseer outputs genome-wide association results as tab-delimited tables, whereas
Phandango expects GWAS or score tracks in a specific `.plot` format: #CHR SNP BP minLOG10(P) log10(p) r^2


This utility converts Pyseer SNP-level association output (from linear models or
linear mixed models) into a format that can be loaded directly into Phandango
alongside a reference genome annotation and phylogeny.

---

## What the script does

`pyseer_to_phandango_plot.py` performs the following steps:

1. Reads a Pyseer result file (tab-delimited with a header).
2. Extracts genomic positions (BP) from the `variant` field by splitting on a delimiter
   (default: `_`).
3. Reads a specified p-value column (auto-detected or user-defined).
4. Computes `-log10(p)` for each variant.
5. Writes a Phandango-compatible `.plot` file.

The script is intended for **SNP-based GWAS outputs** from Pyseer
(default linear model or `--lmm`).

---

## Requirements

- Python ≥ 3.8  
- No external Python dependencies

---

## Input format

### Pyseer output

The input file must be a tab-delimited Pyseer output file with a header, typically
containing:

- `variant` — SNP identifier (e.g. `AE017143.1_12345_A_T`)
- a p-value column (e.g. `lrt-pvalue`)

Example variant format:

contig_position_reference_alternative
AE017143.1_12345_A_T


In this example, the genomic position (`BP`) is `12345`.

---

## Usage

## Basic usage

```bash
python pyseer_to_phandango_plot.py \
  --pyseer infection_SNPs_lmm.txt \
  --out infection_SNPs_lmm.plot
```
This uses default assumptions:
variant column name
_ as delimiter
second field (index 1) as genomic position
auto-detection of the p-value column

# Reproducing the original bash workflow (recommended)
The following command reproduces the original shell-based pipeline used in the analysis:

python pyseer_to_phandango_plot.py \
  --pyseer infection_SNPs_lmm.txt \
  --out infection_SNPs_lmm.plot \
  --chr-index 26 \
  --variant-col variant \
  --pcol lrt-pvalue \
  --variant-delim "_" \
  --bp-field-index 1 \
  --skip-nonpositive-p
  
This corresponds to:
extracting BP using cut -d"_" -f2,
using column 4 (lrt-pvalue) as the p-value,
computing -log10(p),
assigning constant values for CHR, SNP, and r^2.

# Variant formats with different delimiters
If variants use a different delimiter or layout:
```bash
python pyseer_to_phandango_plot.py \
  --pyseer results.txt \
  --out results.plot \
  --variant-delim "|" \
  --bp-field-index 1
```

# Specifying a p-value column explicitly
If auto-detection fails:
```bash
python pyseer_to_phandango_plot.py \
  --pyseer results.txt \
  --out results.plot \
  --pcol lrt-pvalue
```

# Output

The script writes a single output file:
*.plot — Phandango-compatible GWAS track
This file can be loaded directly into Phandango together with:
* a reference genome annotation (GFF or GenBank),
* a phylogenetic tree,
* associated metadata.

# Notes and safeguards

Extremely small p-values may appear as 0 in Pyseer output due to rounding.
These are handled safely and converted to large -log10(p) values.
Rows with non-parsable genomic positions or invalid p-values can be skipped using:
--skip-nonpositive-p
--allow-missing-bp
The script does not alter statistical results and should be considered a
visualisation utility only.

# Example workflow context

This script was used to visualise SNP-level GWAS results generated with Pyseer
(using both linear models and linear mixed models) in Phandango, alongside
recombination-filtered phylogenies and reference genome annotations.

# License
MIT License

# Citation
