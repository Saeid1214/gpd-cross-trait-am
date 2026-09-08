#!/usr/bin/env Rscript
#
# 01_heatmap.R
#
# Visualizes the cross-trait GPD matrix (scripts/01_cross_trait_gpd_matrix.R)
# as a heatmap. Cells significant after Bonferroni correction are outlined
# and labeled in bold; non-significant cells are shown in muted gray.
#
# Usage:
#   Rscript plots/01_heatmap.R
#
# Output:
#   <plot_dir>/cross_trait_gpd_heatmap.pdf
#   <plot_dir>/cross_trait_gpd_heatmap.png

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

results_dir <- "results"
plot_dir    <- "plots/output"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

trait_labels <- fread("config/trait_manifest.tsv", header = TRUE) %>%
  select(trait_id, label)

cross <- fread(file.path(results_dir, "cross_trait_gpd_matrix.tsv"), header = TRUE)

data1 <- merge(cross, trait_labels, by.x = "trait_1", by.y = "trait_id") %>%
  rename(label_1 = label)
data <- merge(data1, trait_labels, by.x = "trait_2", by.y = "trait_id") %>%
  rename(label_2 = label)

# substitute your own preferred display order for trait labels
desired_order <- sort(unique(c(data$label_1, data$label_2)))
data$label_1 <- factor(data$label_1, levels = desired_order)
data$label_2 <- factor(data$label_2, levels = desired_order)

data$sig_est    <- ifelse(data$significant, data$estimate, NA_real_)
data$nonsig_est <- ifelse(!data$significant, data$estimate, NA_real_)

plot <- ggplot(data, aes(x = label_2, y = label_1, fill = estimate)) +
  geom_tile(color = NA) +
  geom_tile(data = subset(data, significant), color = "black", fill = NA, linewidth = 0.7) +
  scale_fill_gradient2(
    low = "orange2", mid = "white", high = "darkolivegreen",
    midpoint = 0, limits = c(-0.25, 0.25), na.value = "gray", guide = "colourbar"
  ) +
  geom_text(aes(label = round(sig_est, 3)), fontface = "bold", size = 3, na.rm = TRUE) +
  geom_text(aes(label = round(nonsig_est, 3)), size = 3, color = "gray70", na.rm = TRUE) +
  theme(
    text = element_text(size = 11, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, vjust = 1, hjust = 0),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5, linetype = "solid"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.title = element_text(size = 10, hjust = 0, vjust = 0.5),
    legend.text = element_text(size = 7, hjust = 0, vjust = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.background = element_rect(fill = "white"),
    plot.margin = unit(c(0.2, 1.7, 0.2, 0.2), "cm")
  ) +
  scale_x_discrete(expand = c(0, 0), limits = rev, position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "white", colour = NA), panel.grid = element_blank())

ggsave(file.path(plot_dir, "cross_trait_gpd_heatmap.pdf"), plot, width = 14, height = 12, device = cairo_pdf)
ggsave(file.path(plot_dir, "cross_trait_gpd_heatmap.png"), plot, width = 14, height = 12, device = cairo_pdf)

message("Done. Heatmap written to ", plot_dir)
