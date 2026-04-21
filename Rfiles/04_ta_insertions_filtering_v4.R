require(openxlsx)
require(dplyr)

# -------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------
thr        <- 500     # count threshold
max_gap_bp <- 50      # max distance (bp) between contiguous TA sites within same groups

# 1. Define the count columns exactly as in your dataset
count_cols <- c(
  # X pools
  "XA","XB","XC","XD","XE","XF","XG","XH",
  # Y pools
  paste0("Y", 1:12),
  # Z row pools
  "ZA","ZB","ZC","ZD","ZE","ZF","ZG","ZH",
  # Z column pools
  paste0("Z", 1:12)
)

# Convenience subsets
x_cols  <- c("XA","XB","XC","XD","XE","XF","XG","XH")
y_cols  <- paste0("Y", 1:12)
zl_cols <- c("ZA","ZB","ZC","ZD","ZE","ZF","ZG","ZH")
zn_cols <- paste0("Z", 1:12)

setwd("~/UF/Dental/rob/ana/")

### Reading dataset
#origin <- read.xlsx("./ta_insertions_annotated_all_samples.xlsx", sheet="Genes")
#colnames(origin)[c(2,3)] <- c("pos1","pos2")
#colnames(origin)[c(6,7)] <- c("start","end")

### Groups columns (for sanity check)
origin[1,11:18]
origin[1,19:30]
origin[1,31:38]
origin[1,39:50]

data.ta <- data.frame(origin, check01=0, check02=0)

### Loop for count 4 and 5
for ( i.row in 1:nrow(data.ta)) {
  
  if( length(which(data.ta[i.row,11:18]>thr))==1 &
      length(which(data.ta[i.row,19:30]>thr))==1 &
      length(which(data.ta[i.row,31:38]>thr))==1 &
      length(which(data.ta[i.row,39:50]>thr))==1 ) {
    data.ta[i.row,51] <- 1
  }
  
  g1 <- length(which(data.ta[i.row,11:18] > thr))  # Xletters
  g2 <- length(which(data.ta[i.row,19:30] > thr))  # Ynumbers
  g3 <- length(which(data.ta[i.row,31:38] > thr))  # Zletters
  g4 <- length(which(data.ta[i.row,39:50] > thr))  # Znumbers
  
  if( (g1 + g2 + g3 + g4) == 5  && g1 >= 1 && g2 >= 1 && g3 >= 1 && g4 >= 1 ) {
    data.ta[i.row,52] <- 1
  }
}

data.ta.4 <- data.ta[data.ta[,51]==1,]
data.ta.5 <- data.ta[data.ta[,52]==1,]

# -------------------------------------------------------------------
# Add relative position of TA within gene: (pos1 - start) / (end - start)
# -------------------------------------------------------------------
data.ta.4$rel_pos <- (data.ta.4$pos1 - data.ta.4$start) / (data.ta.4$end - data.ta.4$start)
data.ta.5$rel_pos <- (data.ta.5$pos1 - data.ta.5$start) / (data.ta.5$end - data.ta.5$start)

# Plate mapping: (Zletter, Znumber) -> plate_pos 1..96
z_rows <- c("A","B","C","D","E","F","G","H")
z_cols <- 1:12

plate_map <- expand.grid(
  Zletter = z_rows,
  Znumber = z_cols
)
# Ensure row-wise order: A1..A12, B1..B12, ..., H1..H12
plate_map <- plate_map[order(plate_map$Zletter, plate_map$Znumber), ]
plate_map$plate_pos <- seq_len(nrow(plate_map))  # 1..96

### FINAL FILTERING FOR 4 and 5 

# + ---------------------------------------------------------------------------------------------------- + 
# Extract Zletter index (which of the 8 letter pools was >thr)
data.ta.4$Zletter <- z_rows[
  apply(data.ta.4[,31:38], 1, function(x) which(x > thr))
]
# Extract Znumber index (which of the 12 number pools was >thr)
data.ta.4$Znumber <- apply(data.ta.4[,39:50], 1, function(x) which(x > thr))

# NEW: Xlet/Ynum for 4-group (always exactly one hit per group by design)
data.ta.4$Xlet <- z_rows[
  apply(data.ta.4[,11:18], 1, function(x) which(x > thr))
]
data.ta.4$Ynum <- apply(data.ta.4[,19:30], 1, function(x) which(x > thr))

# ---- Add plate position using the mapping ----
data.ta.4 <- merge(
  data.ta.4,
  plate_map,
  by = c("Zletter", "Znumber"),
  all.x = TRUE
)

cols <- setdiff(names(data.ta.4), c("Zletter", "Znumber", "plate_pos"))
data.ta.4 <- data.ta.4[, c(cols, "Zletter", "Znumber", "plate_pos")]
data.ta.4 <- data.ta.4[with(data.ta.4, order(pos1)),]

### For each gene, pick the TA insertion with the largest total counts
data.ta.4.best <- data.ta.4 %>%
  mutate(
    total_reads = rowSums(across(all_of(count_cols)), na.rm = TRUE)
  ) %>%
  group_by(geneID) %>%
  slice_max(total_reads, n = 1, with_ties = FALSE) %>%  # one best TA per gene
  ungroup()

data.ta.4.best <- data.ta.4.best[with(data.ta.4.best, order(start)),]

# + ---------------------------------------------------------------------------------------------------- +
# Helper for 5-group confidence: TRUE if group is ambiguous (top two hits too similar)
ambig_flag_fun <- function(x, thr) {
  idx <- which(x > thr)
  if (length(idx) <= 1) return(FALSE)   # 0 or 1 hit -> not ambiguous
  vals <- x[idx]
  if (length(vals) < 2) return(FALSE)
  vals <- sort(vals, decreasing = TRUE)
  ratio <- vals[1] / vals[2]
  # ambiguous if top two counts are too similar (ratio < 2)
  return(ratio < 2)
}

## ---- Decode Zletter for data.ta.5 ----
## 1 hit  > thr  -> use that one
## 2 hits > thr  -> keep the one with the larger count
data.ta.5$Zletter <- apply(data.ta.5[, 31:38], 1, function(x) {
  idx <- which(x > thr)
  
  if (length(idx) == 1) {
    return(z_rows[idx])
  } else if (length(idx) == 2) {
    # keep the one with the larger count
    idx_keep <- idx[which.max(x[idx])]
    return(z_rows[idx_keep])
  } else {
    return(NA_character_)  # safety fallback
  }
})

## ---- Decode Znumber for data.ta.5 ----
## same logic but we keep the numeric index 1..12
data.ta.5$Znumber <- apply(data.ta.5[, 39:50], 1, function(x) {
  idx <- which(x > thr)
  
  if (length(idx) == 1) {
    return(idx)
  } else if (length(idx) == 2) {
    idx_keep <- idx[which.max(x[idx])]
    return(idx_keep)
  } else {
    return(NA_integer_)    # safety fallback
  }
})

## ---- NEW: Xlet/Ynum for 5-group with tie-breaking ----
data.ta.5$Xlet <- apply(data.ta.5[,11:18], 1, function(x) {
  idx <- which(x > thr)
  if (length(idx) == 1) {
    z_rows[idx]
  } else if (length(idx) == 2) {
    idx_keep <- idx[which.max(x[idx])]
    z_rows[idx_keep]
  } else {
    NA_character_
  }
})

data.ta.5$Ynum <- apply(data.ta.5[,19:30], 1, function(x) {
  idx <- which(x > thr)
  if (length(idx) == 1) {
    idx
  } else if (length(idx) == 2) {
    idx_keep <- idx[which.max(x[idx])]
    idx_keep
  } else {
    NA_integer_
  }
})

## ---- NEW: Confidence metric + ambiguity flag for 5-group ----
ambig_X  <- apply(data.ta.5[,11:18],  1, ambig_flag_fun, thr = thr)
ambig_Y  <- apply(data.ta.5[,19:30],  1, ambig_flag_fun, thr = thr)
ambig_ZL <- apply(data.ta.5[,31:38],  1, ambig_flag_fun, thr = thr)
ambig_ZN <- apply(data.ta.5[,39:50],  1, ambig_flag_fun, thr = thr)

data.ta.5$ambiguous_5group <- ambig_X | ambig_Y | ambig_ZL | ambig_ZN

data.ta.5$conf_5group <- pmax(
  0,
  1 - 0.25 * (
    as.integer(ambig_X) +
      as.integer(ambig_Y) +
      as.integer(ambig_ZL) +
      as.integer(ambig_ZN)
  )
)

## ---- Add plate position using the mapping ----
data.ta.5 <- merge(
  data.ta.5,
  plate_map,
  by = c("Zletter", "Znumber"),
  all.x = TRUE
)

## ---- Put Zletter, Znumber, plate_pos at the end (like we did for data.ta.4) ----
cols <- setdiff(names(data.ta.5), c("Zletter", "Znumber", "plate_pos"))
data.ta.5 <- data.ta.5[, c(cols, "Zletter", "Znumber", "plate_pos")]

data.ta.5 <- data.ta.5[with(data.ta.5, order(pos1)),]

### For each gene, pick the TA insertion with the largest total counts
data.ta.5.best <- data.ta.5 %>%
  mutate(
    total_reads = rowSums(across(all_of(count_cols)), na.rm = TRUE)
  ) %>%
  group_by(geneID) %>%
  slice_max(total_reads, n = 1, with_ties = FALSE) %>%  # one best TA per gene
  ungroup()

data.ta.5.best <- data.ta.5.best[with(data.ta.5.best, order(start)),]

# -------------------------------------------------------------------
# NEW: Find contiguous cross-gene TA clusters within same (Xlet,Ynum,Zletter,Znumber)
# starting from data.ta.4 and data.ta.5
# -------------------------------------------------------------------
find_neighbor_insertions <- function(df, max_gap = 50) {
  # df must contain: insertion, geneID, pos1, Xlet, Ynum, Zletter, Znumber
  
  df_work <- df[order(df$Xlet, df$Ynum, df$Zletter, df$Znumber, df$pos1), ]
  
  if (nrow(df_work) < 2) return(character(0))
  
  # group by plate well (Xlet,Ynum,Zletter,Znumber)
  grp <- interaction(df_work$Xlet, df_work$Ynum,
                     df_work$Zletter, df_work$Znumber,
                     drop = TRUE)
  
  keep_idx <- logical(nrow(df_work))
  
  for (g in split(seq_len(nrow(df_work)), grp)) {
    if (length(g) < 2) next
    
    # walk through this group to build contiguous clusters
    current_cluster <- g[1]
    for (k in 2:length(g)) {
      i_prev <- g[k-1]
      i_curr <- g[k]
      if (df_work$pos1[i_curr] - df_work$pos1[i_prev] <= max_gap) {
        # still contiguous in genomic space
        current_cluster <- c(current_cluster, i_curr)
      } else {
        # cluster ended; check if cross-gene
        if (length(unique(df_work$geneID[current_cluster])) > 1) {
          keep_idx[current_cluster] <- TRUE
        }
        current_cluster <- g[k]
      }
    }
    # last cluster in this group
    if (length(unique(df_work$geneID[current_cluster])) > 1) {
      keep_idx[current_cluster] <- TRUE
    }
  }
  
  df_work$insertion[keep_idx]
}

neighbors_ins4 <- find_neighbor_insertions(data.ta.4, max_gap = max_gap_bp)
neighbors_ins5 <- find_neighbor_insertions(data.ta.5, max_gap = max_gap_bp)

# + ---------------------------------------------------------------------------------------------------- +
### EXCEL WRITING (reorder only at export so rel_pos is column 11)
# + ---------------------------------------------------------------------------------------------------- +

# Front columns: rel_pos immediately after annotation (col 11)
front10 <- c("insertion","pos1","pos2","geneID","name",
             "start","end","length","strand","annotation")
front_with_rel <- c(front10, "rel_pos")

# Build clean export data frames, dropping check01/check02/total_reads
df4      <- subset(data.ta.4,      select = -c(check01, check02))
df4_best <- subset(data.ta.4.best, select = -c(check01, check02, total_reads))
df5      <- subset(data.ta.5,      select = -c(check01, check02))
df5_best <- subset(data.ta.5.best, select = -c(check01, check02, total_reads))

reorder_for_export <- function(df) {
  keep_front <- front_with_rel[front_with_rel %in% names(df)]
  
  # tail block: plate-related + confidence columns, in this order if present
  tail_cols <- c(
    "Xlet","Ynum","Zletter","Znumber","plate_pos",
    "ambiguous_5group","conf_5group"
  )
  tail_cols <- tail_cols[tail_cols %in% names(df)]
  
  # everything else (non-front, non-count, non-tail)
  middle <- setdiff(names(df), c(keep_front, count_cols, tail_cols))
  
  df[, c(keep_front, count_cols, middle, tail_cols)]
}

df4      <- reorder_for_export(df4)
df4_best <- reorder_for_export(df4_best)
df5      <- reorder_for_export(df5)
df5_best <- reorder_for_export(df5_best)

# Neighbor tables: same columns as df4 / df5, just subset by insertion ID
neighbors4 <- df4[df4$insertion %in% neighbors_ins4, ]
neighbors5 <- df5[df5$insertion %in% neighbors_ins5, ]

# Create a blank workbook
OUT <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(OUT, "4 across groups")
addWorksheet(OUT, "best 4 for each gene")
addWorksheet(OUT, "5 across groups")
addWorksheet(OUT, "best 5 for each gene")
addWorksheet(OUT, "neighbors 4")
addWorksheet(OUT, "neighbors 5")

# Write the data to the sheets
writeData(OUT, sheet="4 across groups",      x = df4)
writeData(OUT, sheet="best 4 for each gene", x = df4_best)
writeData(OUT, sheet="5 across groups",      x = df5)
writeData(OUT, sheet="best 5 for each gene", x = df5_best)
writeData(OUT, sheet="neighbors 4",          x = neighbors4)
writeData(OUT, sheet="neighbors 5",          x = neighbors5)

# Setting columns width (first up to 55 columns)
for (nm in c("4 across groups","best 4 for each gene",
             "5 across groups","best 5 for each gene",
             "neighbors 4","neighbors 5")) {
  df_name <- switch(
    nm,
    "4 across groups"       = "df4",
    "best 4 for each gene"  = "df4_best",
    "5 across groups"       = "df5",
    "best 5 for each gene"  = "df5_best",
    "neighbors 4"           = "neighbors4",
    "neighbors 5"           = "neighbors5"
  )
  df <- get(df_name)
  ncols <- ncol(df)
  max_cols <- min(55, ncols)
  
  # base widths: default 7
  widths <- rep(7, max_cols)
  if (max_cols >= 1)  widths[1]  <- 14          # insertion
  if (max_cols >= 2)  widths[2:9] <- 9          # pos1,pos2,geneID,name,start,end,length,strand
  if (max_cols >= 10) widths[10] <- 75          # annotation
  if (max_cols >= 11) widths[11] <- 9           # rel_pos
  
  # plate + confidence columns get width 8 if within first max_cols
  tail_for_width <- c("Xlet","Ynum","Zletter","Znumber","plate_pos",
                      "ambiguous_5group","conf_5group")
  idx_tail <- match(tail_for_width, names(df))
  idx_tail <- idx_tail[!is.na(idx_tail) & idx_tail <= max_cols]
  if (length(idx_tail) > 0) {
    widths[idx_tail] <- 8
  }
  
  setColWidths(
    OUT,
    sheet  = nm,
    cols   = 1:max_cols,
    widths = widths
  )
}

# Reorder worksheets (keep main four first, then neighbors)
worksheetOrder(OUT) <- c(1,2,3,4,5,6)

### Conditional Formatting: counts > thr
blueStyle <- createStyle(bgFill = "#80D1FF")

df_list <- list(
  "4 across groups"      = df4,
  "best 4 for each gene" = df4_best,
  "5 across groups"      = df5,
  "best 5 for each gene" = df5_best,
  "neighbors 4"          = neighbors4,
  "neighbors 5"          = neighbors5
)

# Center both horizontally and vertically
centerHV <- createStyle(
  halign = "center",
  valign = "center"
)
# Left horizontal only (used to override column 10)
centerH <- createStyle(
  halign = "left"
)

for (sheet in names(df_list)) {
  df <- df_list[[sheet]]
  # count columns by name
  idx_counts <- match(count_cols, names(df))
  idx_counts <- idx_counts[!is.na(idx_counts)]
  
  if (length(idx_counts) > 0) {
    conditionalFormatting(
      OUT,
      sheet = sheet,
      cols  = idx_counts,
      rows  = 2:(nrow(df) + 1),
      type  = "expression",
      style = blueStyle,
      rule  = paste0(">", thr)
    )
  }
  
  ### Center all
  addStyle(
    OUT,
    sheet = sheet,
    style = centerHV,
    rows = 1:(nrow(df)+1),   # include header row
    cols = 1:ncol(df),
    gridExpand = TRUE,
    stack = TRUE
  )
  ### Recenter column 10 (left horizontal)
  addStyle(
    OUT,
    sheet = sheet,
    style = centerH,
    rows = 1:(nrow(df)+1),   # include header row
    cols = 10,
    gridExpand = TRUE,
    stack = TRUE
  )
  
  ### Conditional Formatting for the relative position of the TA insertion
  conditionalFormatting(
    OUT,
    sheet = sheet,
    cols = c(11),
    rows = 2:(nrow(df)+1),  # exclude header row
    type = "colorScale",
    style = c("green", "yellow", "red"),
    rule = quantile(df[,11], c(0, 0.5, 1), na.rm = TRUE)  # explicitly set scale values
  )
}

# + ---------------------------------------------------------------------------------------------------- +
### Define color of the groups columns (X, Y, Z-letter, Z-number)
# + ---------------------------------------------------------------------------------------------------- +
header_colors <- list(
  X  = "#DAEEF4",
  Y  = "#B7DEE8",
  ZL = "#FDE9D9",
  ZN = "#FBD5B4"
)

for (sheet in names(df_list)) {
  df <- df_list[[sheet]]
  
  x_idx  <- match(x_cols,  names(df)); x_idx  <- x_idx[!is.na(x_idx)]
  y_idx  <- match(y_cols,  names(df)); y_idx  <- y_idx[!is.na(y_idx)]
  zl_idx <- match(zl_cols, names(df)); zl_idx <- zl_idx[!is.na(zl_idx)]
  zn_idx <- match(zn_cols, names(df)); zn_idx <- zn_idx[!is.na(zn_idx)]
  
  if (length(x_idx) > 0) {
    addStyle(OUT, sheet = sheet,
             style = createStyle(fgFill = header_colors$X),
             rows = 1, cols = min(x_idx):max(x_idx),
             gridExpand = TRUE, stack = TRUE)
  }
  if (length(y_idx) > 0) {
    addStyle(OUT, sheet = sheet,
             style = createStyle(fgFill = header_colors$Y),
             rows = 1, cols = min(y_idx):max(y_idx),
             gridExpand = TRUE, stack = TRUE)
  }
  if (length(zl_idx) > 0) {
    addStyle(OUT, sheet = sheet,
             style = createStyle(fgFill = header_colors$ZL),
             rows = 1, cols = min(zl_idx):max(zl_idx),
             gridExpand = TRUE, stack = TRUE)
  }
  if (length(zn_idx) > 0) {
    addStyle(OUT, sheet = sheet,
             style = createStyle(fgFill = header_colors$ZN),
             rows = 1, cols = min(zn_idx):max(zn_idx),
             gridExpand = TRUE, stack = TRUE)
  }
}

# + ---------------------------------------------------------------------------------------------------- +
### Define color of the plate-related columns (rel_pos, Xlet, Ynum, Zletter, Znumber, plate_pos, ambig/conf)
# + ---------------------------------------------------------------------------------------------------- +
greenStyle <- createStyle(fgFill = "#DAE4C0")
plate_cols <- c("rel_pos","Xlet","Ynum","Zletter","Znumber","plate_pos")

for (sheet in names(df_list)) {
  df <- df_list[[sheet]]
  p_idx <- match(plate_cols, names(df))
  p_idx <- p_idx[!is.na(p_idx)]
  
  if (length(p_idx) > 0) {
    addStyle(
      OUT,
      sheet = sheet,
      style = greenStyle,
      rows = 1,
      cols = p_idx,
      gridExpand = TRUE,
      stack = TRUE
    )
  }
  # For 5-group-related sheets, also shade ambiguous/conf columns (wherever they are)
  if (sheet %in% c("5 across groups","best 5 for each gene","neighbors 5")) {
    extra_idx <- match(c("ambiguous_5group","conf_5group"), names(df))
    extra_idx <- extra_idx[!is.na(extra_idx)]
    if (length(extra_idx) > 0) {
      addStyle(
        OUT,
        sheet = sheet,
        style = greenStyle,
        rows = 1,
        cols = extra_idx,
        gridExpand = TRUE,
        stack = TRUE
      )
    }
  }
}

# Export the file
saveWorkbook(OUT,
             paste0("./data.ta/ta_insertions_annotated_all_samples_filtered_",
                    thr,"_v4_neighbors.xlsx"),
             overwrite=TRUE)
