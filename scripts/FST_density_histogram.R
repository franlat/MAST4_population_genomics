# ============================================================
# FST histograms for MAST-4 species (A, B, C, E)
# - Reads pairwise FST matrices
# - Subsets by station list
# - Converts negatives/NA to 0
# - Plots histograms per species + combined density
# ============================================================

library(tidyverse)
library(ggplot2)
library(ggpubr)

# -----------------------------
# User parameters
# -----------------------------
bins_hist   <- 200
x_limits    <- c(0, 0.6)
y_limits_h  <- c(0, 75)

# Genetic differentiation thresholds (Hartl & Clark)
fst_thresholds <- c(0.05, 0.15, 0.25)
fst_colors     <- c("green", "orange", "red")

# Replace station labels (keeps your original convention)
station_pattern <- ".SUR.GGMM."
station_repl    <- "RA"

# -----------------------------
# Input files
# -----------------------------
fst_file_a <- "Documents/Biocluster/PGs/SNPs/POGENOM/MAST4/TA.MAST4A.noN.c10.s4.fst.txt"
fst_file_b <- "Documents/Biocluster/PGs/SNPs/POGENOM/MAST4/TA.MAST4B.noN.c10.s4.fst.txt"
fst_file_c <- "Documents/Biocluster/PGs/SNPs/POGENOM/MAST4/TA.MAST4C.noN.c10.s4.fst.txt"
fst_file_e <- "Documents/Biocluster/PGs/SNPs/POGENOM/MAST4/TA.MAST4E.noN.c10.s4.fst.txt"

# -----------------------------
# Station lists (Tara surface metagenomes)
# -----------------------------
sampfil_a <- c("TA_SUR_GGMM_102","TA_SUR_GGMM_11","TA_SUR_GGMM_111","TA_SUR_GGMM_112","TA_SUR_GGMM_124","TA_SUR_GGMM_131","TA_SUR_GGMM_132","TA_SUR_GGMM_135","TA_SUR_GGMM_136","TA_SUR_GGMM_137","TA_SUR_GGMM_138","TA_SUR_GGMM_142","TA_SUR_GGMM_143","TA_SUR_GGMM_145","TA_SUR_GGMM_146","TA_SUR_GGMM_147","TA_SUR_GGMM_148","TA_SUR_GGMM_149","TA_SUR_GGMM_150","TA_SUR_GGMM_151","TA_SUR_GGMM_152","TA_SUR_GGMM_16","TA_SUR_GGMM_18","TA_SUR_GGMM_20","TA_SUR_GGMM_22","TA_SUR_GGMM_23","TA_SUR_GGMM_25","TA_SUR_GGMM_30","TA_SUR_GGMM_32","TA_SUR_GGMM_4","TA_SUR_GGMM_5","TA_SUR_GGMM_51","TA_SUR_GGMM_6","TA_SUR_GGMM_64","TA_SUR_GGMM_65","TA_SUR_GGMM_66","TA_SUR_GGMM_68","TA_SUR_GGMM_7","TA_SUR_GGMM_72","TA_SUR_GGMM_76","TA_SUR_GGMM_78","TA_SUR_GGMM_80","TA_SUR_GGMM_81","TA_SUR_GGMM_83","TA_SUR_GGMM_9","TA_SUR_GGMM_93","TA_SUR_GGMM_95","TA_SUR_GGMM_96","TA_SUR_GGMM_97","TA_SUR_GGMM_98")

sampfil_b <- c("TA_SUR_GGMM_124","TA_SUR_GGMM_130","TA_SUR_GGMM_131","TA_SUR_GGMM_132","TA_SUR_GGMM_138","TA_SUR_GGMM_41","TA_SUR_GGMM_42","TA_SUR_GGMM_43","TA_SUR_GGMM_45","TA_SUR_GGMM_46","TA_SUR_GGMM_51")

sampfil_c <- c("TA_SUR_GGMM_100","TA_SUR_GGMM_102","TA_SUR_GGMM_106","TA_SUR_GGMM_109","TA_SUR_GGMM_11","TA_SUR_GGMM_110","TA_SUR_GGMM_113","TA_SUR_GGMM_123","TA_SUR_GGMM_124","TA_SUR_GGMM_125","TA_SUR_GGMM_128","TA_SUR_GGMM_129","TA_SUR_GGMM_130","TA_SUR_GGMM_131","TA_SUR_GGMM_132","TA_SUR_GGMM_136","TA_SUR_GGMM_137","TA_SUR_GGMM_138","TA_SUR_GGMM_139","TA_SUR_GGMM_142","TA_SUR_GGMM_143","TA_SUR_GGMM_32","TA_SUR_GGMM_34","TA_SUR_GGMM_36","TA_SUR_GGMM_38","TA_SUR_GGMM_39","TA_SUR_GGMM_41","TA_SUR_GGMM_42","TA_SUR_GGMM_43","TA_SUR_GGMM_45","TA_SUR_GGMM_46","TA_SUR_GGMM_5","TA_SUR_GGMM_51","TA_SUR_GGMM_52","TA_SUR_GGMM_58","TA_SUR_GGMM_64","TA_SUR_GGMM_65","TA_SUR_GGMM_7","TA_SUR_GGMM_72","TA_SUR_GGMM_9")

sampfil_e <- c("TA_SUR_GGMM_135","TA_SUR_GGMM_145","TA_SUR_GGMM_146","TA_SUR_GGMM_147","TA_SUR_GGMM_148","TA_SUR_GGMM_149","TA_SUR_GGMM_150","TA_SUR_GGMM_151","TA_SUR_GGMM_152","TA_SUR_GGMM_6","TA_SUR_GGMM_66","TA_SUR_GGMM_81","TA_SUR_GGMM_82","TA_SUR_GGMM_83","TA_SUR_GGMM_89","TA_SUR_GGMM_93")

# -----------------------------
# Helper functions
# -----------------------------

# Convert station IDs from the list to match matrix colnames/rownames
to_matrix_ids <- function(x) gsub("_", "\\.", x)

# Convert matrix row/col station labels to shorter plot labels (RAxx convention)
to_plot_station <- function(x) gsub(station_pattern, station_repl, x)

# Read, subset, clean, and convert to long format
read_fst_long <- function(fst_file, sample_ids, mast4_label) {
  mat_ids <- to_matrix_ids(sample_ids)

  fst_tb <- read.table(fst_file, header = TRUE, stringsAsFactors = TRUE) %>%
    filter(rownames(.) %in% mat_ids) %>%
    select(all_of(mat_ids))

  fst_df <- fst_tb %>%
    mutate(
      Sample = to_plot_station(rownames(fst_tb)),
      MAST4  = mast4_label
    )

  # Replace negative values and NAs (keeps original behavior)
  fst_df[fst_df < 0] <- 0
  fst_df[is.na(fst_df)] <- 0

  fst_df %>%
    pivot_longer(
      cols = -c(MAST4, Sample),
      names_to  = "Stations",
      values_to = "Fst"
    ) %>%
    mutate(Stations = to_plot_station(Stations))
}

# Plot histogram for one species
plot_fst_hist <- function(fst_long, fill_color) {
  ggplot(fst_long, aes(x = Fst)) +
    geom_histogram(alpha = 0.3, color = NA, fill = fill_color, bins = bins_hist) +
    geom_vline(xintercept = fst_thresholds, color = fst_colors) +
    theme(legend.position = "none") +
    coord_cartesian(xlim = x_limits, ylim = y_limits_h)
}

# -----------------------------
# Build long tables + plots
# -----------------------------
fst_a_long <- read_fst_long(fst_file_a, sampfil_a, "A")
fst_b_long <- read_fst_long(fst_file_b, sampfil_b, "B")
fst_c_long <- read_fst_long(fst_file_c, sampfil_c, "C")
fst_e_long <- read_fst_long(fst_file_e, sampfil_e, "E")

p_a <- plot_fst_hist(fst_a_long, "deepskyblue1")
p_b <- plot_fst_hist(fst_b_long, "burlywood4")
p_c <- plot_fst_hist(fst_c_long, "chartreuse3")
p_e <- plot_fst_hist(fst_e_long, "darkmagenta")

# 2x2 panel
ggarrange(
  p_a, p_b, p_c, p_e,
  labels = c("MAST4A", "MAST4B", "MAST4C", "MAST4E"),
  ncol = 2, nrow = 2
)

# -----------------------------
# Combined density plot
# -----------------------------
all_fst <- bind_rows(fst_a_long, fst_b_long, fst_c_long, fst_e_long)

p_all <- ggplot(all_fst, aes(x = Fst)) +
  geom_density(alpha = 0.3, color = NA, fill = "black") +
  geom_vline(xintercept = fst_thresholds, color = fst_colors) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 6.5))

p_all
