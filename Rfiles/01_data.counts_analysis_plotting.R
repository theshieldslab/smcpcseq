### This script will write several filtered dataset at different filtering parameters
### These were later plotted by script 02_plotting_counts_threshold.R
require(dplyr)
require(tidyr)
require(stringr)
require(purrr)
require(openxlsx)

setwd("~/UF/Dental/rob/ana/")

samples <- read.delim("samples.txt", header=F)
colnames(samples) <- "locus_tag"

data.counts <- read.delim("./counts/htseq_counts_all_samples.txt", header=F, sep= ' ')
colnames(data.counts) <- c("locus_tag", samples[,1])
data.counts <- data.counts[, c(1, 2:9, 13:21, 10:12, 34:41, 25:33, 22:24)]
data.counts <- data.counts[-which(substr(data.counts[,1],1,2)=="__"),]


### ---------------------------
### 1. Define the factor levels that exist in YOUR data
###    (based on your column names and what we inspected)
### ---------------------------

# X group (8 levels)
X_levels  <- c("XA","XB","XC","XD","XE","XF","XG","XH")

# Y group (12 levels)
Y_levels  <- c("Y1","Y2","Y3","Y4","Y5","Y6","Y7","Y8","Y9","Y10","Y11","Y12")

# ZL = Z letter group (8 levels)
ZL_levels <- c("ZA","ZB","ZC","ZD","ZE","ZF","ZG","ZH")

# ZN = Z number group
# Important: "Z4" and "Z5" are merged in your file header as "Z4Z5".
# So right now we treat "Z4Z5" as a single observed condition.
# That gives us 11 total unique ZN conditions.
ZN_levels <- c("Z1","Z2","Z3","Z4","Z5","Z6","Z7","Z8","Z9","Z10","Z11","Z12")

### ---------------------------
### 2. Build the full combination grid
###    This gives all combos of:
###    X (8) x Y (12) x ZL (8) x ZN (11) = 8448 rows
### ---------------------------
all_combos <- expand.grid(
  X  = X_levels,
  Y  = Y_levels,
  ZL = ZL_levels,
  ZN = ZN_levels,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

### ---------------------------
### 3. Reshape your raw data.counts into long format
###    Assumptions:
###    - data.counts exists in your environment
###    - first column is locus_tag (gene ID like "SMU_1006")
###    - remaining columns are the condition columns:
###         XA, XB, ..., XH,
###         Y1...Y12,
###         ZA...ZH,
###         Z1, Z2, Z3, Z4Z5, Z6...Z12
### ---------------------------

counts_long <- data.counts %>%
  pivot_longer(
    cols      = -locus_tag,
    names_to  = "condition",
    values_to = "count"
  ) %>%
  mutate(
    # group = first letter of condition ("X","Y","Z")
    group = str_sub(condition, 1, 1),
    
    # We'll label which subgroup each condition belongs to:
    subgroup = dplyr::case_when(
      condition %in% X_levels      ~ "X",
      condition %in% Y_levels      ~ "Y",
      condition %in% ZL_levels     ~ "ZL",
      condition %in% ZN_levels     ~ "ZN",
      TRUE                         ~ "OTHER"   # safety net
    )
  )

# After this step, counts_long looks like:
# locus_tag | condition | count | group | subgroup
# SMU_1006  | XA        | 269   | "X"   | "X"
# SMU_1006  | Y10       | 111474| "Y"   | "Y"
# SMU_1006  | ZA        | 137   | "Z"   | "ZL"
# SMU_1006  | Z10       | 48600 | "Z"   | "ZN"
# etc.

### ---------------------------
### 4. Define a helper that, for ONE gene, attaches its observed counts
###    to all 8*12*8*11 = 8448 theoretical combinations
### ---------------------------

expand_gene <- function(gene_df) {
  
  gene_id <- unique(gene_df$locus_tag)
  
  # Build lookup tables for each subgroup (X, Y, ZL, ZN)
  x_lookup <- gene_df %>%
    filter(subgroup == "X") %>%
    select(condition, count)
  
  y_lookup <- gene_df %>%
    filter(subgroup == "Y") %>%
    select(condition, count)
  
  zl_lookup <- gene_df %>%
    filter(subgroup == "ZL") %>%
    select(condition, count)
  
  zn_lookup <- gene_df %>%
    filter(subgroup == "ZN") %>%
    select(condition, count)
  
  # Join these lookups into the full combo grid
  gene_expanded <- all_combos %>%
    # Attach X_count for that gene
    left_join(x_lookup,  by = c("X"  = "condition")) %>%
    rename(X_count = count) %>%
    # Attach Y_count
    left_join(y_lookup,  by = c("Y"  = "condition")) %>%
    rename(Y_count = count) %>%
    # Attach ZL_count
    left_join(zl_lookup, by = c("ZL" = "condition")) %>%
    rename(ZL_count = count) %>%
    # Attach ZN_count
    left_join(zn_lookup, by = c("ZN" = "condition")) %>%
    rename(ZN_count = count) %>%
    # Add gene ID and placeholder for a true joint measurement
    mutate(
      locus_tag   = gene_id,
      joint_count = NA_real_
    ) %>%
    select(
      locus_tag,
      X, Y, ZL, ZN,
      X_count, Y_count, ZL_count, ZN_count,
      joint_count
    )
  
  # Explanation:
  # - X_count is this gene's observed count in that X condition alone.
  # - Y_count is this gene's observed count in that Y condition alone.
  # - ZL_count is this gene's observed count in that Z-letter condition alone.
  # - ZN_count is this gene's observed count in that Z-number condition alone.
  # - joint_count stays NA because you do NOT have direct experimental data
  #   for that exact (X,Y,ZL,ZN) four-way combination.
  
  return(gene_expanded)
}

### ---------------------------
### 5. Apply that helper to every gene
###    Result: one giant data frame with 8448 rows per gene
### ---------------------------

expanded_all_genes <- counts_long %>%
  group_by(locus_tag) %>%
  group_split() %>%
  map_df(expand_gene)


### ---------------------------
### 6. Quick QC examples
### ---------------------------

# How many unique combos per gene?
expanded_all_genes %>%
  filter(locus_tag == first(locus_tag)) %>%
  nrow()
# Expect 8448 with current levels (8 * 12 * 8 * 11)

# Peek at one gene, first few rows:
expanded_all_genes %>%
  filter(locus_tag == first(locus_tag)) %>%
  head()

# Check that we didn't lose anything weird
summary(expanded_all_genes)


### 
### 7.
### 

# Choose a similarity tolerance
# Example: CV < 0.2 (20% variability)
cv_threshold <- 0.2
count_threshold <- 500

expanded_all_genes <- expanded_all_genes %>%
  rowwise() %>%
  mutate(
    # compute mean and sd across the four counts
    mean_count = mean(c(X_count, Y_count, ZL_count, ZN_count), na.rm = TRUE),
    sd_count   = sd(c(X_count, Y_count, ZL_count, ZN_count), na.rm = TRUE),
    
    # coefficient of variation
    cv = ifelse(mean_count > 0, sd_count / mean_count, NA_real_),
    
    # TRUE/FALSE if all 4 are within the allowed CV
    similar_by_cv = !is.na(cv) & cv <= cv_threshold,
    
    # also check absolute range if you prefer Option A:
    range_count = max(c(X_count, Y_count, ZL_count, ZN_count), na.rm = TRUE) -
      min(c(X_count, Y_count, ZL_count, ZN_count), na.rm = TRUE),
    similar_by_range = range_count <= cv_threshold * mean_count  # ≤ 20 % spread
  ) %>%
  ungroup()


for (cv_threshold in c(0.05,0.1,0.15,0.2,0.25)) {
  
  for (count_threshold in c(500,1000,1500,2000,2500)) {
    
    expanded_all_genes2 <- expanded_all_genes[expanded_all_genes[,11]>count_threshold & expanded_all_genes[,13]<cv_threshold,]
    write.xlsx(expanded_all_genes2, paste0("./output/ta_insertions_per_gene_cv-",cv_threshold,"_counts-",count_threshold,"_npositions_",nrow(expanded_all_genes2),".xlsx") )
  }
}


