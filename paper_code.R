# Paper: Heinze. Murdoch et al. 2026 - 
# Script summary: Publication script combing all code needed.

# Loading packages
library(tidyr)
library(tidyverse)
library(readxl)
library(writexl)
library(dplyr)
library(stringr)
library(here)
library(ggalluvial)
library(ggplot2)
library(pheatmap)
library(tibble)
library(matrixStats)
library(reshape2)
library(RColorBrewer)
library(glmnet)
library(caret)
library(mclust)   
library(irr)       
library(gridExtra)
library(vegan)     
library(limma)
library(ggrepel)
library(ggplotify)
library(rcompanion)
library(forcats)
library(patchwork)
library(ggsignif)
library(scales)
library(purrr)
library(VennDiagram)
library(ggVennDiagram)
library(fgsea)
library(msigdbr)
library(lubridate)
library(ggpubr)

# X) ----- SMALL HELPER FUNCTIONS -----
# Helper to ensure folder exists ----
dir_create_safe <- function(path) if(!dir.exists(path)) dir.create(path, recursive=TRUE)

outdir <- "analysis_outputs"
dir_create_safe(outdir)
dir_create_safe(file.path(outdir,"plots"))
dir_create_safe(file.path(outdir,"tables"))

# To generate alluvial plot per comparison ----
make_alluvial_plot <- function(df, site_x, site_y, title_letter, output_name) {
  message("Processing ", site_x, " → ", site_y)
  
  combined_df <- df %>% filter(new_sample_site_s1 == site_x, new_sample_site_s2 == site_y) %>% select(Patient_ID, pred_s1, pred_s2)
  
  if (nrow(combined_df) == 0) {message("No pairs found for ", site_x, " → ", site_y)
    return(NULL)}
  
  # Contingency table
  contingency_matrix <- table(combined_df$pred_s1, combined_df$pred_s2)
  data_long <- as.data.frame(as.table(contingency_matrix))
  colnames(data_long) <- c("pred_s1", "pred_s2", "Freq")
  
  # Calculate each PrOTYPE Frequency per cohort
  total_cases <- sum(data_long$Freq)
  data_long <- data_long %>% mutate(Percent = 100 * Freq / total_cases)
  
  # Axis labeling logic
  x_labels <- if (site_x == site_y) {c(paste0(site_x, " one"), paste0(site_y, " two"))
  } else {c(site_x, site_y)}
  
  # Plot updated
  plot <- ggplot(data_long, aes(axis1 = pred_s1, axis2 = pred_s2, y = Freq)) +
    geom_alluvium(aes(fill = pred_s1), alpha = 0.7, width = 1/12) +
    geom_stratum(width = 1/12, fill = "grey90", color = "grey40") +
    geom_text(stat = "stratum", aes(label = paste0("\n", scales::percent(after_stat(prop), accuracy = 0.1))), 
              size = 3.5, color = "black", fontface = "bold") +
    scale_x_discrete(limits = x_labels, expand = c(0.15, 0.05)) +
    scale_fill_manual(values = adnexa_colors, name = "Subtype") +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), legend.position = "right", axis.title = element_blank(),
          axis.ticks = element_blank(), axis.text.y = element_blank(),
          plot.title = element_text(size = 15, face = "bold")) +
    labs(title = title_letter,
         #subtitle = "Alluvial flow based on anatomical laterality", 
         y = "Number of Paired Samples")

  ggsave(paste0(output_name, "_alluvial.png"), plot = plot, width = 7, height = 5, units = "in", dpi = 300)
  
  return(plot)
}

# ---------------------------------------------------------------------

# A) ----- DATA PREPARATION -----
# Loading data frames if needed - pulled in from private git
df_PrOTYPE_combined <- read.csv(here::here("output/PrOTYPE_combined_20260218.csv"))
scaled_dataframe <- df_PrOTYPE_combined %>%
  select(14:68) %>% #gene expression only; different range since different dataframe w/ diff. order of columns
  scale()
# Clipping values to set boundaries
clip_values <- function(dat, hi = 2, lo = -2) {
  dat0 <- dat
  dat0[dat0 > hi] <- hi
  dat0[dat0 < lo] <- lo
  return(dat0)
}
dat_scaled_clipped <- clip_values(scaled_dataframe, hi = 2, lo = -2)
# Transform format of dataframe for heatmaps
all.matrix <- t(dat_scaled_clipped)

## Further preparing frame 
#expected objects: df_PrOTYPE_combined and all.matrix (genes x samples)
if(!exists("df_PrOTYPE_combined") || !exists("all.matrix")) stop("Load 'df_PrOTYPE_combined' and 'all.matrix' before running script.")

df.all <- df_PrOTYPE_combined #creating shorthand version
# Check matrix orientation: rows = genes, cols = samples
if(!is.matrix(all.matrix)) all.matrix <- as.matrix(all.matrix)

# Make sure samples align, set colnames of matrix from 'df.all' if missing
if(is.null(colnames(all.matrix))) {
  if("SampleID_clean" %in% colnames(df.all)) {
    colnames(all.matrix) <- df.all$SampleID_clean
  } else stop("Matrix has no colnames and 'df' lacks SampleID_clean — set sample names")
}
# Align ordering: ensure colnames(all.matrix) match df.all$SampleID_clean
if(!all(colnames(all.matrix) %in% df.all$SampleID_clean)) {
  stop("Some sample names in 'all.matrix' are not found in 'df.all'/'df_PrOTYPE_combined$SampleID_clean'. Correct the naming.")
}
# Reorder 'df.all' dataframe to match matrix columns
df.all <- df.all[match(colnames(all.matrix), df.all$SampleID_clean), ]
rownames(df.all) <- df.all$SampleID_clean

# Infer patient id column
patient_col <- "Patient_ID"
# Infer SampleID
sample_col <- "SampleID_clean"
# Check pairing counts
pairs_count <- table(df.all[[patient_col]])
pairs_num <- as.numeric(pairs_count) #patient counts as numbers
message("Number of patients: ", length(unique(df.all[[patient_col]])), 
        "; sample counts per patient (summary): ", paste(summary(pairs_num), collapse = ", "))
# Seeing how many are of each pair
counts_table <- sort(table(pairs_num), decreasing = TRUE)
message("Patient sample-count distribution (samples per patient):", paste((counts_table), collapse = ", " ))# 2 > 135 patients, 3 > 3 patients
# ----------------------------------------------------

# B) ----- QC & exploratory (### SUPPL. FIGURE 2A & 2B) ----- 
# PCA (samples) 
pca <- prcomp(t(all.matrix), center=TRUE, scale.=TRUE)
# Calculate variance
var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
# Building graph data frame
pc_df <- data.frame(Sample = colnames(all.matrix), PC1 = pca$x[,1], PC2 = pca$x[,2],
                    pred = df.all$pred,
                    #subtype = df.all$subtype, #ONCE C1.MES to C5.PRO values are available for all
                    #entropy = df.all$entropy, 
                    site = df.all$new_sample_site, patient = df.all[[patient_col]])
# Create label for each plot
pc1_lab <- paste0("PC1 (", round(var_explained[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_explained[2], 1), "%)")
# Create PCA plots
pca1 <- ggplot(pc_df, aes(PC1, PC2, color = pred)) + geom_point(size=2) + theme_minimal() + 
  labs(x= pc1_lab, y=pc2_lab) + ggtitle("PCA colored by PrOTYPE subtype")
pca2 <- ggplot(pc_df, aes(PC1, PC2, color = site, shape = pred)) + geom_point(size=2) + theme_minimal() + 
  labs(x= pc1_lab, y=pc2_lab) + ggtitle("PCA colored by site, shaped by PrOTYPE")

# Assess variance per gene
gene_var <- rowVars(all.matrix)
top_var_genes <- names(sort(gene_var, decreasing=TRUE))[1:50]
annotation_col <- data.frame(pred = df.all$pred)
rownames(annotation_col) <- df.all$SampleID_clean
# ------------------------------------------------------

# C) ----- Paired limma for site effect  (### SUPPL. TABLE S3) -----
## Goal: per-gene site effect controlling for patient pairing (using paired limma)
# Require site variable: assume df.all$new_sample_site and that "adnexa" is reference or set comparison "other vs adnexa"
if(!("new_sample_site" %in% colnames(df.all))) stop("'df.all' must contain new_sample_site for site comparisons.")

## Create design matrix: will do model ~ 0 + site with blocking on patient using duplicateCorrelation
# Only considering everything versus adnexa
sitef <- factor(df.all$new_sample_site)
design <- model.matrix(~0 + sitef)
colnames(design) <- levels(sitef)
block <- factor(df.all[[patient_col]])

## Limma pipeline, using normalized data directly
# Duplicate correlation
corfit <- duplicateCorrelation(all.matrix, design, block=block)
message("DupCor consensus: ", corfit$consensus)

fit <- lmFit(all.matrix, design, block=block, correlation = corfit$consensus)

# Define contrasts: contrast each non-adnexa vs adnexa (if adnexa exists), else pick first level as reference
ref <- if("adnexa" %in% levels(sitef)) "adnexa" else levels(sitef)[1]
other_levels <- setdiff(levels(sitef), ref)
contrasts_list <- paste0(other_levels, "_vs_", ref)
contrast_expr <- paste0(other_levels, "-", ref)

# Keep only contrasts where both groups exist
site_counts <- table(sitef)
valid_contrasts <- sapply(other_levels, function(x) sum(sitef == x) > 0 & sum(sitef == ref) > 0)
other_levels <- other_levels[valid_contrasts]

if(length(other_levels) == 0){
  stop("No valid site contrasts with at least one sample in each group.")
}
contrast_expr <- paste0(other_levels, "-", ref)
contrast_names <- paste0(other_levels, "_vs_", ref)
cm <- makeContrasts(contrasts = contrast_expr, levels = design)
colnames(cm) <- contrast_names

fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
# Save top tables for each contrast
limma_results <- list()
for(i in seq_along(contrast_names)){
  ct <- contrast_names[i]
  top <- topTable(fit2, coef=i, number=Inf, adjust.method="BH")
  top$gene <- rownames(top)
  limma_results[[ct]] <- top
  
  fname <- paste0("limma_site_", ct, ".csv")
  write.csv(top, file.path(outdir, "tables", fname), row.names = FALSE)
}

# Save summary of DE counts
de_summary <- sapply(limma_results, function(tt) sum(tt$adj.P.Val < 0.05, na.rm=TRUE))
write.csv(data.frame(contrast=names(de_summary), nDE=de_summary), file.path(outdir,"tables","limma_site_DE_counts.csv"))
# -----------------------------------------------------

# D) ----- DATA PREP - Subtype concordance and agreement metrics -----
# Ensure subtype = factor
df.all$pred <- factor(df.all$pred)
## Build pairs: for each patient, find the samples and compare PrOTYPE across pairs
pair_counts <- table(df.all[[patient_col]]) # can also be used from above
pair_ids <- names(pair_counts[pair_counts == 2])
message("Number of patients with exactly 2 samples: ", length(pair_ids))
pairs_df <- df.all %>% filter((.data[[patient_col]] %in% pair_ids)) %>% arrange(.data[[patient_col]])

# Creating wide table, FORMAT: patient x (subtype_sample1, subtype_sample2) 
meta_cols <- c("sample", "SampleID_clean", "sample.origin", "sample.site", "new_sample_site",  
               "primary_site_flag", "secondary_site_flag", "Patient_ID", "pred",
               "C1.MES", "C2.IMM", "C4.DIF", "C5.PRO", "entropy", "pred_primary",
               "shifted", "Dataset")
# Identify expression columns (everything in df.all except meta_cols)
expr_cols <- setdiff(colnames(df.all), meta_cols)

##### STEPWISE #####
# 1) Make expression dataframe: samples x genes, with SampleID_clean column
expr_df <- df.all[, c(sample_col, expr_cols)] #sample_col from line 92
# Make meta dataframe:
meta_df <- df.all[, setdiff(colnames(df.all), expr_cols)]  #metadata dataframe
# Ensure SampleID_clean column type matches meta
expr_df[[sample_col]] <- as.character(expr_df[[sample_col]])
meta_df[[sample_col]] <- as.character(meta_df[[sample_col]])
meta_df[[patient_col]] <- as.character(meta_df[[patient_col]])

# 2) Split metadata into primary and secondary by primary_site_flag
primary_meta <- meta_df %>% filter(primary_site_flag == TRUE)
secondary_meta <- meta_df %>% filter(primary_site_flag == FALSE)

# 3) Join expressions to each meta table
primary_join <- primary_meta %>% left_join(expr_df, by = sample_col)
secondary_join <- secondary_meta %>% left_join(expr_df, by = sample_col)

# 4) Detect and handle duplicates: >1 primary / >1 secondary per patient
dup_prim <- primary_join %>% group_by(!!sym(patient_col)) %>% summarise(n = n(), .groups = "drop") %>% filter(n > 1)
dup_sec <- secondary_join %>% group_by(!!sym(patient_col)) %>% summarise(n = n(), .groups = "drop") %>% filter(n > 1)
# Checking specimen numbers per case and simply sticking to one per patient
if(nrow(dup_prim) > 0) {
  warning("Some patients have >1 primary sample. Keeping the first row per patient. See dup_prim for list.")
  print(dup_prim)
  # Keep first by SampleID_clean (or some other deterministic rule)
  primary_join <- primary_join %>% arrange(!!sym(patient_col), !!sym(sample_col)) %>%
    group_by(!!sym(patient_col)) %>% slice(1) %>% ungroup()
}
if(nrow(dup_sec) > 0) {
  warning("Some patients have >1 secondary sample. Keeping the first row per patient. See dup_sec for list.")
  print(dup_sec)
  secondary_join <- secondary_join %>% arrange(!!sym(patient_col), !!sym(sample_col)) %>%
    group_by(!!sym(patient_col)) %>% slice(1) %>% ungroup()
}

# 5) Rename expression columns to have suffixes _s1 and _s2
primary_expr_cols_s1 <- paste0(expr_cols, "_s1")
secondary_expr_cols_s2 <- paste0(expr_cols, "_s2")
primary_join <- primary_join %>% rename_at(vars(all_of(expr_cols)), ~ primary_expr_cols_s1)
secondary_join <- secondary_join %>% rename_at(vars(all_of(expr_cols)), ~ secondary_expr_cols_s2)

# 6) Also rename metadata columns to suffix _s1/_s2 to keep them distinct
meta_cols_keep <- c(sample_col, "pred", "new_sample_site", "sample.origin", "sample.site", "entropy")
# intersect(meta_cols_keep, colnames(meta_df))  # safety check
primary_join <- primary_join %>% rename_at(vars(all_of(meta_cols_keep)), ~ paste0(., "_s1"))
secondary_join <- secondary_join %>% rename_at(vars(all_of(meta_cols_keep)), ~ paste0(., "_s2"))

# 7) Join primary and secondary by patient id (full join to keep unmatched)
wide_pairs <- full_join(primary_join, secondary_join, by = "Patient_ID", suffix = c("_s1", "_s2"))
#### ----- ####

# Find all genes that are usable and are in s1/s2 pairs
gene_cols_s1 <- grep("_s1$", colnames(wide_pairs), value = TRUE)
gene_names <- sub("_s1$", "", gene_cols_s1)
gene_cols_s2 <- paste0(gene_names, "_s2")
# Only keep genes that have both s1 and s2 columns >>> since we have meta labeled as _s1/_s2
valid_genes <- expr_cols
# Compute differences for all genes
for (g in valid_genes) {
  col1 <- paste0(g, "_s1")
  col2 <- paste0(g, "_s2")
  diff_col <- paste0(g, "_diff")
  wide_pairs[[diff_col]] <- wide_pairs[[col2]] - wide_pairs[[col1]]
}
#To check for many genes diff were calcultated
cat("Calculated paired differences for", length(valid_genes), "genes.\n")
summary(wide_pairs[[paste0(valid_genes[1], "_diff")]])
# --------------------------------------------------

# E) ----- LATERAL ASSESSMENT (FIGURE 1B, 1C & 1D, TABLE in FIGURE 1A) -----
# 1) Prepare data
# Extract adnexa only samples
adnexa_lr <- wide_pairs %>% filter(new_sample_site_s1 == "adnexa",
                                   new_sample_site_s2 == "adnexa") %>%
  # Extract laterality from sample_site columns
  filter(!is.na(sample.site_s1), !is.na(sample.site_s2))

adnexa_lr_ordered <- adnexa_lr %>% mutate(pred_left = case_when(sample.site_s1 == "left" ~ pred_s1,
                                                                sample.site_s2 == "left" ~ pred_s2,
                                                                TRUE ~ NA_character_),
                                          pred_right = case_when(sample.site_s1 == "right" ~ pred_s1,
                                                                 sample.site_s2 == "right" ~ pred_s2,
                                                                 TRUE ~ NA_character_))
subtype_levels <- c("C1.MES", "C2.IMM", "C4.DIF", "C5.PRO")
adnexa_lr_ordered <- adnexa_lr_ordered %>% mutate(
  pred_left  = factor(pred_left,  levels = subtype_levels),
  pred_right = factor(pred_right, levels = subtype_levels)) %>%
  filter(!is.na(pred_left), !is.na(pred_right))

# Contingency table
cont_lr <- table(adnexa_lr_ordered$pred_left,adnexa_lr_ordered$pred_right)

lr_long <- as.data.frame(cont_lr)
colnames(lr_long) <- c("Left_Adnexa", "Right_Adnexa", "Freq")

# 2) Visualize with plot
# Alluvial plot 
adnexa_colors <- c("C1.MES" = "#FF1F5B", "C2.IMM" = "#00CD6C", "C4.DIF" = "#009ADE", "C5.PRO" = "#AF58BA")

lr_alluvial_plot <- ggplot(adnexa_lr_ordered, aes(axis1 = pred_left, axis2 = pred_right)) +
  geom_alluvium(aes(fill = pred_left), alpha = 0.7, width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey90", color = "grey40") +
  geom_text(stat = "stratum", aes(label = paste0("\n", scales::percent(after_stat(prop), accuracy = 0.1))),
            size = 3.5, color = "black", fontface = "bold") +
  scale_x_discrete(limits = c("Left Adnexa", "Right Adnexa"), expand = c(0.15, 0.05)) +
  scale_fill_manual(values = adnexa_colors, name = "Subtype") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), legend.position = "right", axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text.y = element_blank(),
        plot.title = element_text(size = 15, face = "bold")) +
  labs(title = "Laterality Concordance in Adnexa", 
       # subtitle = "Alluvial flow based on anatomical laterality",
       y = "Number of Paired Samples")

# 3) Quantify laterality stability (data FIGURE 1A)
## Kappa
kappa_lr <- kappa2(adnexa_lr_ordered[, c("pred_left", "pred_right")])
## Chi-square
tab_lr <- table(adnexa_lr_ordered$pred_left, adnexa_lr_ordered$pred_right)
chisq_lr <- chisq.test(tab_lr)
# Copy into here and filter R/L data  -- #need to use PrOTYPE combined from 20260218 to have sample.site
df_lr <- df_PrOTYPE_combined %>% filter(sample.site %in% c("left", "right")) %>% filter(!is.na(pred))

# Define right/left concordance per patient
lr_collapsed <- df_lr %>% filter(sample.site %in% c("left", "right")) %>%
  group_by(Patient_ID, sample.site, pred) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(Patient_ID, sample.site) %>% slice_max(n, with_ties = FALSE) %>%   # dominant subtype
  ungroup()

lr_concordance <- lr_collapsed %>% pivot_wider(names_from = sample.site, values_from = pred) %>%
  filter(!is.na(left) & !is.na(right)) %>%
  mutate(lr_status = ifelse(left == right, "concordant", "discordant"))
# Check discordant/cordant numbers
table(lr_concordance$lr_status)

df_lr_entropy <- df_lr %>% inner_join(lr_concordance %>% select(Patient_ID, lr_status), by = "Patient_ID")

## Statistical test: entropy vs concordance
wilcox.test(entropy ~ lr_status, data = df_lr_entropy)
## Side-specific entropy comparison (left vs right)
wilcox.test(entropy ~ sample.site, data = df_lr_entropy)

# 4) Visulize stats with concorfance
# Entropy by concordance
p1 <- ggplot(df_lr_entropy, aes(x = lr_status, y = entropy, fill = lr_status)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8, linewidth = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, shape = 21, stroke = 0.3, color = "black") +
  scale_fill_manual(values = c("concordant" = "#F0F0F0", "discordant" = "#E69F00")) +
  theme_classic() +
  labs(title = "Subtype entropy in left/right adnexal samples", x = "Left–Right Concordance",
       y = "Subtype Entropy") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(size = 14, face = "bold", hjust = 0),
        axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 11, color = "black"),
        axis.line = element_line(linewidth = 0.5), panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(min(df_lr_entropy$entropy) * 0.95, max(df_lr_entropy$entropy) * 1.05))

# By right/left entropy
p2 <- ggplot(df_lr_entropy, aes(x = sample.site, y = entropy, fill = sample.site)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8, linewidth = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, shape = 21, stroke = 0.3, color = "black") +
  scale_fill_manual(values = c("left" = "#7570B3", "right" = "#1B9E77")) +
  theme_classic() +
  labs(title = "Subtype Entropy by Laterality", x = "Side", y = "Subtype Entropy") +
  theme(legend.position = "none", plot.title = element_text(size = 14, face = "bold", hjust = 0),
        axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 11, color = "black"),
        axis.line = element_line(linewidth = 0.5), panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(min(df_lr_entropy$entropy) * 0.95, max(df_lr_entropy$entropy) * 1.05))
# --------------------------------------------------

# CSD) Replicate Analysis (SUPPL. FIGURE S3A, S3B & S3C) ------
# -------------------------------------------

# F) ----- ALLULIVAL PLOT FOR NON-ADNEXA (FIGURE 2C, 2D & 2E) -----
if (!require(ggalluvial)) install.packages("ggalluvial")
# Prepare data and load if needed
# Make sure wide_pairs exists and contains the expected variables
if(!exists("wide_pairs")) stop("Load 'wide_pairs' before running script.")
# Optional: Double-check the setup of the loaded dataframe
message("Loaded dataset dimensions: ", nrow(wide_pairs), " samples and ", ncol(wide_pairs), " columns.")

# Create plot for each site pair (relies on helper function from above)
omentum_plot        <- make_alluvial_plot(wide_pairs, "adnexa", "omentum", "  Omentum", "omentum")
peritoneum_plot     <- make_alluvial_plot(wide_pairs, "adnexa", "peritoneum", "  Peritoneum", "peritoneum")
uterus_plot         <- make_alluvial_plot(wide_pairs, "adnexa", "uterus", "  Uterus", "uterus")
ln_plot             <- make_alluvial_plot(wide_pairs, "adnexa", "lymphnode", "  Lymph Node", "lymphnode")
# ---------------------------------------------------------

# G) ---- ENTROPY ASSOCIATION  w/ anatom. sites/shifting pattern (FIGURE 2A & 2B) -----
# 1) Defining variable and data prep
# Non-adnexal metastatic site with molecular subtype shift
meta_shift <- final_combined[
  final_combined$secondary_site != "adnexa" & final_combined$primary_pred != final_combined$secondary_pred,]

# Non-adnexal metastatic site without molecular subtype shift
meta_no_shift <- final_combined[
  final_combined$secondary_site != "adnexa" & final_combined$primary_pred == final_combined$secondary_pred,]

# Double adnexal with molecular subtype shift
adnexa_shift <- final_combined[
  final_combined$secondary_site == "adnexa" & final_combined$primary_pred != final_combined$secondary_pred,]

# Double adnexal without subtype shift
adnexa_no_shift <- final_combined[
  final_combined$secondary_site == "adnexa" &  final_combined$primary_pred == final_combined$secondary_pred,]

# Making a sublist of all shift cases
shift <- rbind(meta_shift, adnexa_shift)
saveRDS(shift, file = "shift.rds")

# 2) Visualize in plot
#### Function to create improved chi-squared plot ####
create_chi_plot <- function(data, title_suffix, filename) {
  # Get observed counts and chi-sq p-value
  observed <- table(data$second_match)
  if (!"TRUE" %in% names(observed)) observed["TRUE"] <- 0
  if (!"FALSE" %in% names(observed)) observed["FALSE"] <- 0
  
  n <- sum(observed)
  expected <- c("Match" = 0.25 * n, "No Match" = 0.75 * n)
  chi_result <- chisq.test(x = observed, p = expected / n)
  
  # Create data frame
  plot_data <- data.frame(
    Status = factor(c("Match", "No Match"), levels = c("Match", "No Match")),
    Observed = c(observed["TRUE"], observed["FALSE"]),
    Expected = expected
  ) %>%
    pivot_longer(cols = c("Observed", "Expected"), names_to = "Type", values_to = "Count")
  
  # Format p-value
  p_label <- case_when(
    chi_result$p.value < 0.001 ~ "*** (p < 0.001)",
    chi_result$p.value < 0.01  ~ "** (p < 0.01)", 
    chi_result$p.value < 0.05  ~ "* (p < 0.05)",
    chi_result$p.value < 0.1   ~ "† (p < 0.1)",
    TRUE                       ~ "ns (p ≥ 0.1)"
  )
  
  # Create plot
  ggplot(plot_data, aes(x = Status, y = Count, fill = Type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
    # geom_text(aes(label = paste0(Count, " (", percent(Count/n, accuracy = 0.1), ")")), 
    #          position = position_dodge(width = 0.8), vjust = -0.3, size = 4) +
    geom_signif(
      comparisons = list(c("Match", "No Match")),
      annotations = p_label,
      y_position = max(plot_data$Count) * 1.05,
      tip_length = 0.01,
      textsize = 5
    ) +
    scale_fill_manual(values = c("Observed" = "#56B4E9", "Expected" = "grey75"), 
                      name = "Frequency") +
    labs(
      title = paste("Prediction Stability:", title_suffix),
      # subtitle = paste("χ² test vs 25% expected match rate | n =", n),
      x = "Prediction Match Status",
      y = "Number of Samples"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 16),
      # plot.subtitle = element_text(size = 12),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = 12)
    )
}
#### --- end function --- ####

# Create Chi-plot for non-adnexal samples
meta_plot <- create_chi_plot(meta_shift, "Non-Adnexal Sites", "meta_chi_plot.png")

# Adding into the boxplot with anatom. sites
data_sorted <- df_PrOTYPE_combined %>%mutate(new_sample_site = factor(new_sample_site, 
                                                                      levels = c("adnexa", "omentum", "peritoneum", "uterus", "lymphnode", "unk")))  # Order for plot

boxplot_entropy1 <- ggplot(data_sorted, aes(x = new_sample_site, y = entropy)) +
  geom_boxplot(alpha = 0.8, width = 0.6, 
               outlier.shape = NA, linewidth = 0.7, median.linewidth = 0.7, fatten = 0) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, stroke = 0.3, shape = 21, color = "black") +
  labs(title = "Entropy by Anatomical Site", 
       subtitle = "Kruskal-Wallis χ²=24.69, p<0.001", x = "Anatomical Site", y = "Subtype Entropy") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(size = 14, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(linewidth = 0.5),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(), panel.spacing = unit(0.5, "lines")
  ) 
#+ coord_cartesian(ylim = c(0, quantile(data_sorted$entropy, 0.98))) # Remove extreme outliers
# ------------------------------------------------------

# H) ----- Entropy associations w/ genes (SUPPL. FIGURE 4A, 4B & 4C) -----
## Creating 3-Per-gene correlation with entropy (vectorized)
# Prepare a vector of entropy values matching the sample order in all.matrix
entropy_vec_long <- c(wide_pairs$entropy_s1, wide_pairs$entropy_s2)
sample_ids_long <- c(wide_pairs$SampleID_clean_s1, wide_pairs$SampleID_clean_s2)
# ID mismatched samples
missing_samples <- setdiff(colnames(all.matrix), sample_ids_long)

# Reorder to match all.matrix column names
entropy_vec_ordered <- entropy_vec_long[match(colnames(all.matrix), sample_ids_long)]
# Check
stopifnot(length(entropy_vec_ordered) == ncol(all.matrix))

# Compute correlations
gene_rho <- apply(all.matrix, 1, function(g) cor(g, entropy_vec_ordered, method="spearman", use="complete.obs"))
gene_p <- apply(all.matrix, 1, function(g) cor.test(g, entropy_vec_ordered, method="spearman")$p.value)

corr_tab <- data.frame(gene = rownames(all.matrix), rho = gene_rho, p = gene_p, p_adj = p.adjust(gene_p, method="BH"))

## Visualize 
# 1 - Volcano plot
corr_tab$neg_log10_p <- -log10(corr_tab$p_adj)
top_hits <- corr_tab[order(corr_tab$p_adj), ][1:10, ]

ggplot(corr_tab, aes(x = rho, y = neg_log10_p)) + geom_point(aes(color = p_adj < 0.05), size = 2.8, alpha = 0.85) +
  geom_text_repel(data = top_hits, aes(label = gene), size = 3.4, box.padding = 0.4, point.padding = 0.25,
                  segment.color = "grey70",  segment.size = 0.3, max.overlaps = Inf, show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray55", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey80", linewidth = 0.5) +
  scale_color_manual(values = c("gray70", "firebrick"), labels = c(`FALSE` = "FDR >= 0.05", `TRUE` = "FDR < 0.05")) +
  labs(x = "Spearman rho", y = expression(-log[10](adjusted~p)),
       color = NULL, title = "Gene-Entropy Correlation (Spearman)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 11),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        legend.box = "vertical") +
  guides(color = guide_legend(title = NULL))


# Ranked barplot of top correlated genes ---
top_genes <- corr_tab[order(corr_tab$p_adj), ][1:21, ]
# Remove entropy from scale
top_genes <- corr_tab %>% arrange(p_adj) %>% slice_head(n = 21) %>% filter(gene != "entropy")
# 2 - Boxplot
ggplot(top_genes, aes(x = reorder(gene, rho), y = rho, fill = rho)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.3) +
  coord_flip() +
  scale_fill_gradient2(low = "#4C78A8", mid = "grey90", high = "#C73E4B", midpoint = 0) +
  labs(x = NULL, y = "Spearman rho", title = "Top 20 Entropy-Associated Genes") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5))

# 3 - Heatmap of top correlated genes
top_genes <- corr_tab[order(corr_tab$p_adj), "gene"][1:25]
expr_sub <- all.matrix[top_genes, , drop = FALSE]  # top 25 genes
# Removing NA lables before plotting
valid_samples <- !is.na(entropy_vec_ordered)
expr_sub_valid <- expr_sub[, valid_samples, drop = FALSE]
entropy_vec_valid <- entropy_vec_ordered[valid_samples]
annot_df <- data.frame(entropy = entropy_vec_valid)
rownames(annot_df) <- colnames(expr_sub_valid)

HT <- pheatmap(expr_sub_valid, annotation_col = annot_df,
               scale = "row", clustering_distance_rows = "correlation", clustering_distance_cols = "euclidean",
               show_rownames = TRUE, show_colnames = FALSE,
               main = "Top 25 Genes by Correlation with Entropy")
# -----------------------------------------------

# I) ----- VARIANCE PARTITION (patient vs site vs residual) (SUPPL. TABLE S4) -----
# Using simple per-gene random effects via lmer or ANOVA style variance components
if(requireNamespace("lme4", quietly=TRUE)) {
  library(lme4)
  # Prepare output table to store prroportion of variances
  varcomp <- data.frame(gene = rownames(all.matrix), patient = NA, site = NA, residual = NA)
  # Creating variance per gene across factors, looping for each gene
  for(i in seq_len(nrow(all.matrix))) {
    g <- as.numeric(all.matrix[i, ])
    dat_g <- data.frame(expr = g, patient = factor(df.all[[patient_col]]), site = factor(df.all$new_sample_site))
    # Scale expr for stability
    dat_g$expr <- scale(dat_g$expr)
    # Fit linear mixed model w/ Expression = overall mean + patient effect + site effect + residual error
    fit_lmm <- try(lmer(expr ~ (1|patient) + (1|site), data = dat_g), silent = TRUE)
    if(inherits(fit_lmm, "merMod")) {
      # Looping extracting variance components
      vc <- as.data.frame(VarCorr(fit_lmm)) 
      # Usual as vc rows: patient, site, residual
      var_pat <- vc$vcov[vc$grp == "patient"]
      var_site <- vc$vcov[vc$grp == "site"] 
      var_res <- vc$vcov[vc$grp == "Residual"]
      total <- var_pat + var_site + var_res
      varcomp$patient[i] <- var_pat/total
      varcomp$site[i] <- var_site/total
      varcomp$residual[i] <- var_res/total
    }
  } # Saving in tables folder as one csv file
  # write.csv(varcomp, file.path(outdir,"tables","variance_partition_by_gene.csv"), row.names=FALSE)
}
# -------------------------------------------------------

# J) ----- VENN DIAGRAM (FIGURE 3A) -----
# If needed import files from output folder
limma_dir <- here::here(outdir, 'tables')

# List all limma site CSV files
limma_files <- list.files(path = limma_dir, pattern = "^limma_site_.*\\.csv$", full.names = TRUE)

# Extract contrast names from filenames
contrast_names <- gsub("^limma_site_|\\.csv$", "", basename(limma_files))

# Read files into a named list
limmaGenes_Reimported <- setNames(lapply(limma_files, read.csv, stringsAsFactors = FALSE), contrast_names)

# using non-adjusted pValues 
pValue_cutoff <- 0.05
# Need to filter our general overview "DE_counts' files
sig_gene_lists <- imap(limmaGenes_Reimported, ~ {
  # Skip "DE_counts" or any table without adj. p-values
  if (!"adj.P.Val" %in% colnames(.x)) {
    message("skipping ", .y, "(no adj. P.Val column)")
    return(NULL)}
  .x %>% filter(P.Value < pValue_cutoff) %>% pull(gene) %>% unique()})

# Remove contrasts with no significant genes
sig_gene_lists <- sig_gene_lists[lengths(sig_gene_lists) > 0] 
sig_gene_long <- imap_dfr(sig_gene_lists, ~ data.frame(contrast = .y, gene = .x, stringsAsFactors = FALSE))

venn.diagram(x = sig_gene_lists, filename = NULL, fill = RColorBrewer::brewer.pal(length(sig_gene_lists), "Set2"),
             alpha = 0.5, cex = 1.2, cat.cex = 1.2)

# Generate Annotated Venn Diagram
# Computing overlap by hand
A <- sig_gene_lists[[1]]
B <- sig_gene_lists[[2]]
C <- sig_gene_lists[[3]]

overlap_AB <- intersect(A, B)
overlap_AC <- intersect(A, C)
overlap_BC <- intersect(B, C)
overlap_ABC <- Reduce(intersect, list(A, B, C))

overlap_AB_only <- setdiff(intersect(A, B), overlap_ABC)
overlap_AC_only <- setdiff(intersect(A, C), overlap_ABC)
overlap_BC_only <- setdiff(intersect(B, C), overlap_ABC)

# Clean color version
venn.plot <- venn.diagram(x = sig_gene_lists[c(1,2,3)],
                          filename = NULL, fill = c("#8DA0CB", "#66C2A5", "#FC8D62"),
                          alpha = 0.45, cex = 1.4, cat.cex = 1.2,  cat.fontface = "bold",  fontfamily = "sans",  lwd = 1.2,
                          col = "grey40")
grid::grid.draw(venn.plot)

# Adding labels manually
grid::grid.text(paste(overlap_AB_only, collapse = ", "), x = 0.75, y = 0.45, gp = grid::gpar(fontsize =10, fontfamily = "sans")) 
grid::grid.text(paste(overlap_AC_only, collapse = ", "), x = 0.25, y = 0.45, gp = grid::gpar(fontsize =10, fontfamily="sans"))
grid::grid.text(paste(overlap_ABC, collapse = ", "), x = 0.5, y = 0.45, gp = grid::gpar(fontsize =10, fontfamily="sans"))
# ------------------------------------------

# K) ----- GENE SET ENRICHMENT ANALYSIS (FIGURE 3B) ------
# Load if needed: 
expr_matrix <- all.matrix
#### Function to automate the output with GSEA ####
run_gsea_analysis <- function(target_value, group_column, metadata, expr_matrix, hallmark_sets, outdir) {
  # Target_value individual variable found in group_column
  # group_column can be whatever you want as discriminating factor
  # Metadata & expr_matrix -  those are the input dataframes used
  message("Running GSEA for: ", target_value, " (group column: ", group_column, ")")
  # 1. Create group factor
  metadata$group <- ifelse(metadata[[group_column]] == target_value, "target", "other") 
  group_factor <- factor(metadata$group)
  
  # 2. Differential expression (limma)
  design <- model.matrix(~0 + group_factor)
  colnames(design) <- levels(group_factor)
  fit <- lmFit(expr_matrix, design)
  contrast <- makeContrasts(target_vs_other = target - other, levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # 3. Extract differential expression results 
  de_results <- topTable(fit2, coef = "target_vs_other", number = Inf, sort.by = "t")
  
  # 4. Prepare ranked gene list for GSEA
  ranked_genes <- de_results$t
  names(ranked_genes) <- rownames(de_results)
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # 5. Run fgsea
  set.seed(123)
  fgsea_res <- fgsea(pathways = hallmark_sets, stats = ranked_genes, nperm = 10000)
  fgsea_res <- fgsea_res %>% arrange(padj)
  
  # 6. Barplot of NES 
  fgsea_plot_data <- fgsea_res %>% arrange(desc(NES))
  BP_plot <- ggplot(fgsea_plot_data, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
    geom_bar(stat = "identity") + coord_flip() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(
      x = "Pathway",
      y = "Normalized Enrichment Score (NES)",
      title = paste("GSEA for", target_value)
    ) +
    theme_minimal()
  
  # Make sure output dirs exist
  dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  
  # Save barplot
  barplot_file <- file.path(outdir, "plots", paste0("GSEA_for_", target_value, "_all.bmp"))
  ggsave(barplot_file, plot = BP_plot, width = 10, height = 7, units = "in", dpi = 300)
  
  # 7. Leading Edge Dotplot
  leading_edge_long <- fgsea_res %>%
    select(pathway, leadingEdge) %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ","))) %>%
    unnest(leadingEdge = strsplit(leadingEdge, ","))
  
  gene_counts <- leading_edge_long %>%
    group_by(leadingEdge) %>%
    summarise(pathway_count = n(), .groups = "drop")
  
  leading_edge_long <- leading_edge_long %>% left_join(gene_counts, by = "leadingEdge")
  
  DP_plot <- ggplot(leading_edge_long, aes(y = pathway, x = leadingEdge)) +
    geom_point(aes(size = pathway_count), alpha = 0.6, color = "steelblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = paste("Leading edge genes for", target_value), x = "Gene", y = "Pathway")
  
  # Save dotplot
  dotplot_file <- file.path(outdir, "plots", paste0("LeadingEdge_", gsub("\\.", "", target_value), ".bmp"))
  ggsave(dotplot_file, plot = DP_plot, width = 10, height = 7, units = "in", dpi = 300)
  
  # 8. Save results
  fgsea_res_to_save <- fgsea_res %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ",")))
  
  table_file <- file.path(outdir, "tables", paste0("GSEA_results_", target_value, ".csv"))
  write.csv(fgsea_res_to_save, table_file, row.names = FALSE)
  
  message("Finished GSEA for: ", target_value)
  
  return(list(
    fgsea_results = fgsea_res,
    barplot = BP_plot,
    dotplot = DP_plot
  ))
}
#### --- end of function --- ####

# GSEA looping > for multiple at once, e.g. anatomical site
# IT will give you an output in console about progression update
targets <- c("adnexa", "omentum", "uterus", "peritoneum") #Anatomical site variable defined
results <- lapply(targets, function(x)
  run_gsea_analysis(x, "new_sample_site", metadata, expr_matrix, hallmark_sets, outdir))

# Visualize & Combine Barplots
combined_plot_BP <- (results[[1]]$barplot + ggtitle(paste("contralateral ", targets[1])) | results[[2]]$barplot + ggtitle(paste(" ", targets[2]))) /
  (results[[3]]$barplot + ggtitle(paste(" ", targets[3])) | results[[4]]$barplot + ggtitle(paste(" ", targets[4])))
# -------------------------------------------------------------------------

# L) ----- COMPARISION DATA vs. Talhouk et al. (SUPPL. TABLE 1) -----
# Input Talhouk data
CRC <- matrix(c(429,447,550,314, 394,389,574,290, 296,91,37,21, 1119,924,1161,625), nrow = 4)
obs <- matrix(c(1119,924,1161,625, 89,83,53,54), nrow = 2, byrow = TRUE)

colnames(obs) <- c("C1.MES", "C2.IMM", "C4.DIF", "C5.PRO")
rownames(obs) <- c("Published", "This cohort")
chisq_result <- chisq.test(obs)

chi_contributions <- (obs - chisq_result$expected)^2 / chisq_result$expected
chi_table <- data.frame(Subtype = colnames(obs),
                        `Published (n=3829)` = obs[1,],
                        Expected_This = round(chisq_result$expected[1,], 1),
                        `χ² contrib` = round(chi_contributions[1,], 2),
                        This_Cohort_n282 = obs[2,],
                        Expected_Pub = round(chisq_result$expected[2,], 1))

# Adnexa only stats
adnexa <-matrix(c(430,447,550,313, 44,66,49,52), nrow = 2, byrow = TRUE)
colnames(adnexa) <- c("C1.MES", "C2.IMM", "C4.DIF", "C5.PRO")
rownames(adnexa) <- c("Published", "This cohort")
chisq_result2 <- chisq.test(adnexa)

chi_contributions2 <- (adnexa - chisq_result2$expected)^2 / chisq_result2$expected
chi_table <- data.frame(Subtype = colnames(adnexa),
                        `Published (n=1740)` = adnexa[1,],
                        Expected_This = round(chisq_result2$expected[1,], 1),
                        `χ² contrib` = round(chi_contributions2[1,], 2),
                        This_Cohort_n211 = adnexa[2,],
                        Expected_Pub = round(chisq_result2$expected[2,], 1))
# -------------------------------------------------------

# M) ----- ADDED-ON DESCRIPTIVES (TABLE 1) -----
# Overall confusion matrix Adnexa-NonAdnexa
df_pairwise <- df_PrOTYPE_combined %>%
  group_by(Patient_ID) %>%
  filter(any(new_sample_site == "adnexa") & any(new_sample_site != "adnexa")) %>%
  ungroup()

df_wide <- df_pairwise %>%
  mutate(site_group = ifelse(new_sample_site == "adnexa", "adnexa", "non_adnexa")) %>%
  select(Patient_ID, site_group, pred) %>%
  distinct() %>%
  pivot_wider(names_from = site_group, values_from = pred)
colnames(df_wide)
conf_mat <- table(Adnexa = df_wide$adnexa, Non_Adnexa = df_wide$non_adnexa)

# Overall confusion matrix Adnexa-Adnexa
df_pairAdnexa <- df_PrOTYPE_combined %>%
  group_by(Patient_ID) %>%
  filter(n() >= 2, all(new_sample_site == "adnexa")) %>%
  ungroup()

df_lateral <- df_pairAdnexa %>%
  filter(sample.site %in% c("left", "right")) %>%
  mutate(max_prob = pmax(C1.MES, C2.IMM, C4.DIF, C5.PRO, na.rm = TRUE)) %>%
  group_by(Patient_ID, sample.site) %>%
  slice_max(max_prob, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Patient_ID, sample.site, pred) %>%
  pivot_wider(names_from = sample.site, values_from = pred)

colnames(df_lateral)
conf_mat <- table(Right = df_lateral$right, Left  = df_lateral$left)
# -----------------------------------------------------

