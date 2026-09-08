# Genomic Signatures of Assortative Mating: Cross-Trait GPD Matrix

Pipeline for detecting **cross-trait genomic signatures of assortative
mating** — testing whether gametic phase disequilibrium (GPD) between one
trait's polygenic score and a *different* trait's polygenic score points to
cross-trait assortment, beyond what within-trait GPD or ordinary genetic
correlation would predict.

This repository accompanies the manuscript *"Genomic Signatures of
Assortative Mating Across Psychiatric and Related Traits in the iPSYCH
Cohort"* and extends the companion repository
[`gpd-within-trait-am`](../gpd-within-trait-am), which this pipeline
depends on for its inputs.

## Background

The within-trait GPD estimate ([see the companion repo](../gpd-within-trait-am))
tests whether a trait's own odd- and even-chromosome polygenic scores are
correlated — a signature of people assorting on that trait specifically.
The **cross-trait** extension asks a broader question: does assortment on
trait A also leave a detectable genomic signature involving trait B's
polygenic score? This can reveal assortative-mating structure that isn't
visible from within-trait estimates or standard genetic correlation alone
(e.g. two traits that are not strongly genetically correlated, but on which
people nonetheless co-assort).

## Method summary

For every **ordered pair** of traits (i, j) in the trait manifest:

1. Take trait *i*'s odd-chromosome PGS and trait *j*'s even-chromosome PGS
   (computed upstream in the within-trait pipeline), restricted to unrelated
   samples.
2. Fit two regressions:
   - `PGS_i,odd  ~ PGS_j,even + 20 PCs(computed from even chromosomes)`
   - `PGS_i,even ~ PGS_j,odd  + 20 PCs(computed from odd chromosomes)`
3. Retain the estimate from whichever regression has the larger predictor
   variance (larger standard error), as in the within-trait method.
4. Repeat for all N² ordered pairs, producing a full trait × trait matrix.
   The matrix is **directional and not necessarily symmetric** — i-on-j and
   j-on-i are different regressions.
5. Apply multiple-testing correction across all pairs tested.

This stage is **O(N²)** in the number of traits and is the slowest stage in
the overall pipeline — with ~30 traits, that's ~900 regressions.

## Repository structure

```
.
├── scripts/
│   └── 01_cross_trait_gpd_matrix.R     # pairwise cross-trait GPD regression
├── plots/
│   └── 01_heatmap.R                    # visualize the matrix as a heatmap
├── config/
│   └── trait_manifest_example.tsv      # example trait manifest (fictional identifiers)
├── results/                            # pipeline output (not tracked; see .gitignore)
├── environment.yml
├── LICENSE
└── README.md
```

## Requirements

- R ≥ 4.2, with `data.table`, `dplyr`, `ggplot2`
- Output from the [within-trait GPD pipeline](../gpd-within-trait-am)
  (specifically, each trait's `<trait_id>_pgs.txt` file and the
  `odd_chr.eigenvec` / `even_chr.eigenvec` PC files)

```bash
conda env create -f environment.yml
conda activate gpd-cross-trait
```

## Usage

```bash
# 1. Compute the full cross-trait GPD matrix
Rscript scripts/01_cross_trait_gpd_matrix.R

# 2. Visualize as a heatmap
Rscript plots/01_heatmap.R
```

Both scripts read the trait list from `config/trait_manifest.tsv` — add or
remove traits there, not in the scripts.

## Data availability

This repository documents methodology, not results reproduction. All paths
and trait identifiers here are placeholders; the underlying iPSYCH data are
held under a restricted-access data use agreement and are not distributed.

## Citation

If you use this pipeline, please cite the manuscript (citation to be added
upon publication) and the original within-trait method:

> Yengo L, et al. (2018). Imprint of assortative mating on the human genome. *PNAS*.

## License

Released under the MIT License — see [LICENSE](LICENSE).
