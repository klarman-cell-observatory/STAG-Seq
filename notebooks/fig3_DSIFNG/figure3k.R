# --- 1. Load Required Libraries ---
# Ensure you have these libraries installed:
# install.packages("tidyverse")
# install.packages("ggrepel")

library(tidyverse) # For ggplot2, dplyr, readr, etc.
library(ggrepel)   # For non-overlapping text labels

# --- 2. Configuration ---
# --- A. Define Genotypes and File Paths ---
# Path to the WIDE precomputed CSV file that contains all comparisons
input_wide_csv_file <- "/mnt/data/project/25_04_29_Figure3_reanalysis/src/analysis_results_v9/deg_vs_AAV_control_wide_table/all_pairwise_deg_stats_vs_AAV_control_wide.csv" # <-- UPDATE THIS PATH

# Define the exact genotype names from the wide table's columns you want to compare
gof1_genotype <- "STAT1_C324R_hom_pure"
gof2_genotype <- "STAT1_D165G_hom_pure"

# Define the condition and reference used for DE analysis (to build column names)
condition <- "IFNG"
reference_genotype <- "AAV_control"

# Set the output file path for the plot
output_dir <- "./Figure3k"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_plot_file <- file.path(output_dir, paste0(gof1_genotype, "_vs_", gof2_genotype, "_LogFC_Scatter.pdf"))

# --- B. Define Plotting Parameters ---
de_logfc_threshold <- 0.25
p_adj_cutoff <- 0.05

# --- C. Define Genes to Label ---
genes_to_label <- c(
    "SOCS3", "IFIT2", "RSAD2", "STING1", "PIM1", "SOCS1", "IFIT3", "CCL8",
    "MX2", "CXCL9", "OAS3", "IL1B", "IL18", "TNFRSF11A", "C5AR1", "MYC", "THBD",
    "IER3", "IRF1", "OASL", "IFITM2", "IFITM3", "ISG15", "IFI44L", "IFI35",
    "EPSTI1", "IFITM1", "STAT2", "IFI44", "CSNK1G2", "MMP9", "EIF4B", "LTA4H",
    "CIITA", "HLA-DPA1", "HLA-DRB5", "HLA-DRB1", "HLA-DRA", "HLA-DQB1",
    "HLA-DQB2", "APOE", "FABP5", "SORT1", "CYP27A1", "NPC1", "CDIPT", "UGCG",
    "SPHK1", "PTGR1"
)


# --- 3. Load and Prepare Data from Wide Table ---
print(paste("Loading wide data from:", input_wide_csv_file))
wide_data <- read_csv(input_wide_csv_file, col_types = cols(.default = "c"))

# Dynamically construct the column names to extract
gof1_logfc_col <- paste(condition, gof1_genotype, "vs", reference_genotype, "logFC", sep = "_")
gof1_pval_col <- paste(condition, gof1_genotype, "vs", reference_genotype, "pval_adj", sep = "_")
gof2_logfc_col <- paste(condition, gof2_genotype, "vs", reference_genotype, "logFC", sep = "_")
gof2_pval_col <- paste(condition, gof2_genotype, "vs", reference_genotype, "pval_adj", sep = "_")

# Select and rename columns to the generic format used for plotting
plot_data <- wide_data %>%
  select(
    gene = names, # Assuming the gene name column is 'names'
    logfc_gof1 = all_of(gof1_logfc_col),
    p_adj_gof1 = all_of(gof1_pval_col),
    logfc_gof2 = all_of(gof2_logfc_col),
    p_adj_gof2 = all_of(gof2_pval_col)
  ) %>%
  # Convert data to numeric and handle potential NA values from selection
  mutate(across(c(logfc_gof1, p_adj_gof1, logfc_gof2, p_adj_gof2), as.numeric)) %>%
  drop_na() %>%
  # Add columns for plotting aesthetics (color, size, and labels)
  mutate(
    color_category = case_when(
      p_adj_gof1 < p_adj_cutoff & p_adj_gof2 < p_adj_cutoff ~ "Significant in Both",
      p_adj_gof1 < p_adj_cutoff ~ paste("Significant in", gof1_genotype, "Only"),
      p_adj_gof2 < p_adj_cutoff ~ paste("Significant in", gof2_genotype, "Only"),
      TRUE ~ "Not Significant"
    ),
    `-log10_p_adj_max` = -log10(pmin(p_adj_gof1, p_adj_gof2, na.rm = TRUE)),
    label = if_else(gene %in% genes_to_label, gene, "")
  )

# Define the order and colors for the significance categories
category_levels <- c(
  "Significant in Both",
  paste("Significant in", gof1_genotype, "Only"),
  paste("Significant in", gof2_genotype, "Only"),
  "Not Significant"
)
plot_data$color_category <- factor(plot_data$color_category, levels = category_levels)

# Split data into significant and non-significant groups for layered plotting
data_sig <- plot_data %>% filter(color_category != "Not Significant")
data_insig <- plot_data %>% filter(color_category == "Not Significant")


# --- 4. Generate Scatter Plot with ggplot2 ---
print("Generating LogFC-LogFC scatter plot with ggplot2 and ggrepel...")

# Define the color palette.
color_palette <- c(
  "Significant in Both" = "green",
  setNames("red", category_levels[2]),
  setNames("blue", category_levels[3]),
  "Not Significant" = "gray"
)

# Create the plot using a layered approach
logfc_plot <- ggplot() +
  geom_vline(xintercept = c(-de_logfc_threshold, de_logfc_threshold), linetype = "dashed", color = "lightgray") +
  geom_hline(yintercept = c(-de_logfc_threshold, de_logfc_threshold), linetype = "dashed", color = "lightgray") +
  geom_vline(xintercept = 0, color = "gray") +
  geom_hline(yintercept = 0, color = "gray") +
  geom_point(data = data_insig, aes(x = logfc_gof1, y = logfc_gof2), color = "lightgray", size = 1, alpha = 0.5) +
  geom_point(data = data_sig, aes(x = logfc_gof1, y = logfc_gof2, color = color_category, size = `-log10_p_adj_max`), alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +
  geom_label_repel(
    data = plot_data,
    aes(x = logfc_gof1, y = logfc_gof2, label = label),
    box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf, segment.color = 'black',
    segment.size = 0.5, min.segment.length = 0, force = 2, force_pull = 0.1,
    size = 2.5, fontface = "bold", label.r = unit(0.15, "lines"),
    label.size = NA, fill = alpha(c("white"), 0.8)
  ) +
  
  # --- 5. Customize Plot Appearance ---
  coord_cartesian(ylim = c(-1.5, 2)) +
  theme_bw() +
  scale_color_manual(values = color_palette, name = "Significance", drop = FALSE) +
  scale_size_continuous(range = c(2, 8), name = "-log10(Adjusted P-value)") +
  labs(
    x = paste("LogFC (", gof1_genotype, " vs WT)"),
    y = paste("LogFC (", gof2_genotype, " vs WT)"),
    title = paste("Differential Gene Expression:", gof1_genotype, "vs", gof2_genotype),
    subtitle = paste(condition, "Stimulated Condition")
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

# Print the plot to the RStudio viewer (if using).
print(logfc_plot)

# --- 6. Save the Plot to a File ---
ggsave(
  output_plot_file,
  plot = logfc_plot,
  width = 12,
  height = 9,
  units = "in",
  dpi = 300
)

print(paste("--- Plot saved successfully to:", output_plot_file, "---"))
