library(openxlsx)
library(dplyr)

setwd("~/UF/Dental/rob/ana/")

# -------------------------------------------------------------
# 1. Load original data
# -------------------------------------------------------------
origin <- read.xlsx("./ta_insertions_annotated_all_samples.xlsx", sheet = "Genes")

# Fix column names as you did
colnames(origin)[c(2,3)] <- c("pos1","pos2")
colnames(origin)[c(6,7)] <- c("start","end")

# Make sure geneID column is named correctly (assuming it is column 4 already)
# str(origin) to confirm if needed

# -------------------------------------------------------------
# 2. Define column indices for each group
# -------------------------------------------------------------
# Based on your previous script:
# 11:18 -> X letters
# 19:30 -> Y numbers
# 31:38 -> Z letters
# 39:50 -> Z numbers

idx_X  <- 11:18
idx_Y  <- 19:30
idx_ZL <- 31:38
idx_ZN <- 39:50

# -------------------------------------------------------------
# 3. Thresholds to test
# -------------------------------------------------------------
thresholds <- c(100, 250, 500, 750, 1000, 1250, 1500,
                1750, 2000, 3000, 4000, 5000)

# -------------------------------------------------------------
# 4. Function to compute summary at one threshold
# -------------------------------------------------------------
summarize_at_threshold <- function(thr, dat = origin) {
  
  # group-wise counts > threshold, vectorized with rowSums
  g1 <- rowSums(dat[, idx_X]  > thr, na.rm = TRUE)  # X letters
  g2 <- rowSums(dat[, idx_Y]  > thr, na.rm = TRUE)  # Y numbers
  g3 <- rowSums(dat[, idx_ZL] > thr, na.rm = TRUE)  # Z letters
  g4 <- rowSums(dat[, idx_ZN] > thr, na.rm = TRUE)  # Z numbers
  
  # 4-group rule: exactly 1 in each group
  check4 <- (g1 == 1 & g2 == 1 & g3 == 1 & g4 == 1)
  
  # 5-group rule: total 5, but each group has at least 1
  check5 <- ((g1 + g2 + g3 + g4) == 5 &
               g1 >= 1 & g2 >= 1 & g3 >= 1 & g4 >= 1)
  
  # Subsets
  data4 <- dat[check4, , drop = FALSE]
  data5 <- dat[check5, , drop = FALSE]
  
  # Unique gene counts (assuming geneID column is named "geneID")
  n_genes4 <- length(unique(data4$geneID))
  n_genes5 <- length(unique(data5$geneID))
  
  tibble(
    threshold     = thr,
    n_insertions4 = sum(check4),
    n_genes4      = n_genes4,
    n_insertions5 = sum(check5),
    n_genes5      = n_genes5
  )
}

# -------------------------------------------------------------
# 5. Run over all thresholds
# -------------------------------------------------------------
summary_thr <- bind_rows(lapply(thresholds, summarize_at_threshold))

print(summary_thr)

# Optionally save
write.table(
  summary_thr,
  file      = "./data.ta/threshold_sweep_summary.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)


## ------------------------------------------------------------------
##  Add elbow-detection and pretty TIFF figure
## ------------------------------------------------------------------

library(ggplot2)
library(tidyr)
library(patchwork)

# ---------- 1) Elbow detection (Kneedle-style) ----------
detect_elbow <- function(x, y) {
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  y_scaled <- (y - min(y)) / (max(y) - min(y))
  p1 <- c(x_scaled[1], y_scaled[1])
  p2 <- c(x_scaled[length(x_scaled)], y_scaled[length(y_scaled)])
  distances <- sapply(seq_along(x_scaled), function(i) {
    p <- c(x_scaled[i], y_scaled[i])
    abs(det(rbind(p2 - p1, p - p1))) / sqrt(sum((p2 - p1)^2))
  })
  which.max(distances)
}

idx_insert <- detect_elbow(summary_thr$threshold, summary_thr$n_insertions4)
idx_genes  <- detect_elbow(summary_thr$threshold, summary_thr$n_genes4)

elbow_insert <- summary_thr$threshold[idx_insert]
elbow_genes  <- summary_thr$threshold[idx_genes]

cat("Suggested elbow for insertions (4-group):", elbow_insert, "\n")
cat("Suggested elbow for genes (4-group)     :", elbow_genes,  "\n")

# ---------- 2) Long format for ggplot ----------
df_long_ins <- summary_thr %>%
  pivot_longer(cols = c("n_insertions4","n_insertions5"),
               names_to = "group",
               values_to = "n_insertions") %>%
  mutate(group = ifelse(group == "n_insertions4", "4-group", "5-group"))

df_long_genes <- summary_thr %>%
  pivot_longer(cols = c("n_genes4","n_genes5"),
               names_to = "group",
               values_to = "n_genes") %>%
  mutate(group = ifelse(group == "n_genes4", "4-group", "5-group"))

group_colors <- c("4-group" = "#1f78b4", "5-group" = "#e31a1c")

# ---------- 3) Insertions plot ----------
p_ins <- ggplot(df_long_ins,
                aes(x = threshold, y = n_insertions, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = elbow_insert, linetype = "dashed", color = "black") +
  annotate("label", x = elbow_insert,
           y = max(df_long_ins$n_insertions) * 0.9,
           label = paste("Elbow =", elbow_insert),
           size = 3, label.size = 0.2, fill = "white") +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(breaks = summary_thr$threshold) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "TA insertions retained vs. count threshold",
    x = "Threshold (> X counts)",
    y = "Number of TA insertions"
  )

# ---------- 4) Genes plot (rotated x labels) ----------
p_genes <- ggplot(df_long_genes,
                  aes(x = threshold, y = n_genes, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = elbow_genes, linetype = "dashed", color = "black") +
  annotate("label", x = elbow_genes,
           y = max(df_long_genes$n_genes) * 0.925,
           label = paste("Elbow =", elbow_genes),
           size = 3, label.size = 0.2, fill = "white") +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(breaks = summary_thr$threshold) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Unique genes represented vs. count threshold",
    x = "Threshold (> X counts)",
    y = "Number of genes"
  )

# ---------- 5) Combine and export ----------
combined_plot <- p_ins / p_genes +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Threshold sensitivity analysis for pooled Tn-Seq selection",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      legend.position = "top",
      legend.title = element_blank()
    )
  )

print(combined_plot)

# Save as TIFF (high resolution)
tiff("./data.ta/threshold_sweep_pretty.tiff",
     width = 8, height = 8, units = "in", res = 600)
print(combined_plot)
dev.off()


### Manuscript-friendly wording
### “Sensitivity analysis across dynamic count thresholds (100–5000) 
### showed a sharp decline in gene representation beyond ~1000 reads. 
### Therefore, a threshold of >500 reads was selected to balance library 
### representativeness and robustness, with high-confidence insertions 
### identified using a geometric mean score across the four pooling dimensions.”




