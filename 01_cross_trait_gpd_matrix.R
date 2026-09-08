#!/usr/bin/env Rscript
#
# 01_cross_trait_gpd_matrix.R
#
# Extends the within-trait GPD regression to every ordered pair of traits:
# regresses trait i's odd-chromosome PGS on trait j's even-chromosome PGS
# (and vice versa), to detect cross-trait genomic signatures of assortative
# mating. Depends on PGS files produced by the within-trait GPD pipeline
# (gpd-within-trait-am/scripts/04_merge_odd_even_pgs.R).
#
# Usage:
#   Rscript scripts/01_cross_trait_gpd_matrix.R
#
# Output:
#   <results_dir>/cross_trait_gpd_matrix.tsv
#
# Note: this is O(N^2) in the number of traits. With ~30 traits, expect
# ~900 regressions; runtime scales accordingly.

suppressMessages({
  library(data.table)
  library(dplyr)
})

# --- config ------------------------------------------------------------
project_root   <- "/path/to/project"
pgs_root       <- file.path(project_root, "data/pgs")            # from within-trait pipeline
pc_dir         <- file.path(project_root, "work")                 # from within-trait pipeline
unrelated_list <- file.path(project_root, "config/unrelated_ancestry_matched_samples.txt")
trait_manifest <- "config/trait_manifest.tsv"
results_dir    <- "results"

dir.create(results_dir, showWarnings = FALSE)

traits <- fread(trait_manifest, header = TRUE)
traits$trait_id <- as.character(traits$trait_id)

unrel    <- fread(unrelated_list, header = FALSE)
pcs_odd  <- fread(file.path(pc_dir, "odd_chr.eigenvec"),  header = FALSE)
pcs_even <- fread(file.path(pc_dir, "even_chr.eigenvec"), header = FALSE)
pc_cols  <- paste0("V", 3:22)

# --- pre-load and cache each trait's PGS (restricted to unrelated) --------
pgs_cache <- list()
for (i in seq_len(nrow(traits))) {
  trait_id <- traits$trait_id[i]
  pgs_file <- file.path(pgs_root, trait_id, paste0(trait_id, "_pgs.txt"))
  if (!file.exists(pgs_file)) {
    message("Note: no PGS file for ", trait_id, " — pairs involving it will be skipped")
    next
  }
  pgs_cache[[trait_id]] <- fread(pgs_file, header = TRUE) %>% filter(IID %in% unrel$V2)
}

out <- data.frame()

for (i in seq_len(nrow(traits))) {
  trait_i <- traits$trait_id[i]
  if (is.null(pgs_cache[[trait_i]])) next
  pgs_i <- pgs_cache[[trait_i]]

  pgs_i_odd_pc  <- merge(pgs_i[, c("IID", "pgs_odd")],  pcs_even, by.x = "IID", by.y = "V2")
  pgs_i_even_pc <- merge(pgs_i[, c("IID", "pgs_even")], pcs_odd,  by.x = "IID", by.y = "V2")

  for (j in seq_len(nrow(traits))) {
    trait_j <- traits$trait_id[j]
    if (is.null(pgs_cache[[trait_j]])) next
    pgs_j <- pgs_cache[[trait_j]]

    df_odd  <- merge(pgs_i_odd_pc,  pgs_j[, c("IID", "pgs_even")], by = "IID")
    df_even <- merge(pgs_i_even_pc, pgs_j[, c("IID", "pgs_odd")],  by = "IID")

    model_odd  <- lm(paste("pgs_odd ~ pgs_even +",  paste(pc_cols, collapse = " + ")), data = df_odd)
    model_even <- lm(paste("pgs_even ~ pgs_odd +",  paste(pc_cols, collapse = " + ")), data = df_even)

    coef_odd  <- summary(model_odd)$coefficients["pgs_even", , drop = FALSE]
    coef_even <- summary(model_even)$coefficients["pgs_odd", , drop = FALSE]

    v_odd  <- data.frame(model = paste0(trait_i, "_odd~", trait_j, "_even"),
                          trait_1 = trait_i, trait_2 = trait_j, as.data.frame(coef_odd))
    v_even <- data.frame(model = paste0(trait_i, "_even~", trait_j, "_odd"),
                          trait_1 = trait_i, trait_2 = trait_j, as.data.frame(coef_even))

    v <- if (v_even$Std..Error > v_odd$Std..Error) v_even else v_odd
    out <- rbind(out, v)
  }
  message("Completed cross-trait row for ", trait_i)
}

names(out) <- c("model", "trait_1", "trait_2", "estimate", "se", "t_value", "p_value")

# Bonferroni correction across all ordered pairs tested
n_pairs <- nrow(out)
out$p_bonferroni <- pmin(out$p_value * n_pairs, 1)
out$significant  <- out$p_bonferroni < 0.05

fwrite(out, file.path(results_dir, "cross_trait_gpd_matrix.tsv"), quote = FALSE, sep = "\t")
message("Done. Results written to ", file.path(results_dir, "cross_trait_gpd_matrix.tsv"))
