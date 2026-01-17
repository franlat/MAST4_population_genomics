# ============================================================
# FST dendrogram + abundance map + environmental variables
# and PERMANOVA tests explaining FST differentiation
#
# Workflow:
# 1) Read metadata + counts, normalize abundance
# 2) Read intradiversity (pi) and SNP counts
# 3) Read FST matrix, subset to selected stations, clean negatives/NA
# 4) Build UPGMA dendrogram (average linkage)
# 5) Plot dendrogram + intradiversity + temperature + world map
# 6) Run PERMANOVA on FST distances (excluding one sample with NAs)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(phangorn)
  library(ggdendro)
  library(cowplot)
  library(maps)
  library(zoo)   # for na.locf (used in distance matrix step)
})

# -----------------------------
# Arguments (CLI-friendly)
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 9) {
  stop(
    paste0(
      "Usage:\n",
      "  Rscript script.R <tarameta> <bincounts_dir> <binsize> <taracounts> <bcf_input> <binname> <preout> <fst_input> <intradiv>\n"
    ),
    call. = FALSE
  )
}

tarameta   <- args[1]
bincounts  <- args[2]
binsize    <- args[3]
taracounts <- args[4]
bcf_input  <- args[5]
binname    <- args[6]
preout     <- args[7]
fst_input  <- args[8]
intradiv   <- args[9]

# ---- Example local paths (optional; keep commented in repo) ----
# tarameta   <- "/Volumes/My Passport for Mac/Back-up/Biocluster/PGs/SNPs/profiles/ggmm.meta.v2.txt"
# bincounts  <- "/Volumes/My Passport for Mac/Back-up/Biocluster/PGs/SNPs/counts/"
# binsize    <- "/Volumes/My Passport for Mac/Back-up/Biocluster/PGs/SNPs/reads/mast.gsize"
# taracounts <- "/Volumes/My Passport for Mac/Back-up/Biocluster/PGs/SNPs/reads/tara.ggmm.counts"
# bcf_input  <- "/Volumes/My Passport for Mac/Back-up/Biocluster/PGs/SNPs/reads/MAST4A/TA.MAST4A.fb.noN.stats.psc"
# binname    <- "MAST4A"
# fst_input  <- "/Volumes/My Passport for Mac/Back-up/Biocluster/PGs/SNPs/POGENOM/MAST4/TA.MAST4A.noN.c10.s4.fst.txt"
# intradiv   <- "/Volumes/My Passport for Mac/Back-up/Biocluster/PGs/SNPs/POGENOM/MAST4/TA.MAST4A.noN.c10.s4.intradiv.txt"

# Sample to exclude due to missing values (as in your original script)
sample_drop <- "TA11_MS"

# FST reference thresholds shown in dendrogram (dotted lines)
fst_thresholds <- c(0.25, 0.15, 0.05)
fst_thr_colors <- c("red", "orange", "green")

# -----------------------------
# Helper functions
# -----------------------------

# Merge multiple tab-delimited count files into one table
multmerge <- function(mypath, ptrn) {
  filenames <- list.files(path = mypath, pattern = ptrn, full.names = TRUE)
  datalist <- lapply(filenames, function(x) read.table(file = x, header = TRUE, sep = "\t"))
  Reduce(function(x, y) merge(x, y), datalist)
}

# Pretty formatting for numeric axis labels (keeps max decimals observed)
pretty_zero <- function(x) {
  max_dec <- max(nchar(str_extract(x, "\\.[0-9]+")), na.rm = TRUE) - 1
  formatC(x, replace.zero = TRUE, zero.print = "0", digits = max_dec, format = "f", preserve.width = TRUE)
}

# Dendrogram + environmental panels + map
plot_dendro_panels <- function(df, dend, title_text) {

  dend_data <- dendro_data(dend, type = "rectangle")

  segment_data <- with(segment(dend_data),
                       data.frame(x = x, y = y, xend = xend, yend = yend))

  sample_pos <- with(dend_data$labels,
                     data.frame(x_center = x, Samples = as.character(label), height = 0.5))

  # Join sample positions to the plotting data (one row per sample)
  ptdf <- df %>%
    left_join(sample_pos, by = "Samples") %>%
    mutate(
      ymin = x_center - height / 2,
      ymax = x_center + height / 2,
      t_center  = Temperature / 2,
      pi_center = Intra_pi / 2
    )

  axis_limits <- with(sample_pos,
                      c(min(x_center - 0.5 * height), max(x_center + 0.5 * height))) +
    0.1 * c(-1, 1)

  # --- Temperature panel (line + points) ---
  plt_temp <- ggplot(ptdf, aes(x = x_center)) +
    geom_point(aes(y = Temperature), colour = "#277da1") +
    geom_line(aes(y = Temperature), colour = "#277da1") +
    scale_y_continuous(expand = c(0, 0), labels = pretty_zero) +
    scale_x_continuous(
      breaks = sample_pos$x_center,
      labels = sample_pos$Samples,
      limits = axis_limits,
      expand = c(0, 0)
    ) +
    labs(y = "Temperature (°C)", x = "") +
    coord_cartesian(ylim = c(-1, 35)) +
    theme_bw() +
    theme(
      plot.margin = unit(c(0.2, 0.05, 0.2, 0.05), "cm"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4.5),
      panel.grid.minor = element_blank()
    )

  # --- Intra-diversity (pi) panel (line + points) ---
  plt_pi <- ggplot(ptdf, aes(x = x_center)) +
    geom_point(aes(y = Intra_pi), colour = "#277da1") +
    geom_line(aes(y = Intra_pi), colour = "#277da1") +
    scale_y_continuous(expand = c(0, 0), labels = pretty_zero) +
    scale_x_continuous(
      breaks = sample_pos$x_center,
      labels = rep("", nrow(sample_pos)),
      limits = axis_limits,
      expand = c(0, 0)
    ) +
    labs(y = "Intradiversity (π)", x = "") +
    coord_cartesian(ylim = c(0, 0.0007)) +
    theme_bw() +
    theme(
      plot.margin = unit(c(0.2, 0.05, 0.2, 0.05), "cm"),
      axis.text.x = element_blank(),
      panel.grid.minor = element_blank()
    )

  # --- World map panel (points sized by abundance, colored by cluster) ---
  world <- map_data("world") %>%
    filter(long >= -165, long <= 75, lat >= -70, lat <= 50)

  plt_map <- ggplot(ptdf, aes(Longitude, Latitude)) +
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "white", fill = "lightgray", linewidth = 0.1
    ) +
    geom_point(size = 0.3, alpha = 1) +
    geom_point(aes(size = Abundance, colour = as.factor(Cluster_ID)), alpha = 0.3) +
    labs(x = "", y = "") +
    theme_bw() +
    theme(
      axis.ticks = element_blank(),
      legend.position = "right",
      plot.margin = unit(c(0.2, -0.3, 0.2, 0.2), "cm")
    )

  # --- Dendrogram panel ---
  plt_dendr <- ggplot(segment_data) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_hline(yintercept = fst_thresholds, linetype = "dotted", colour = fst_thr_colors) +
    scale_y_continuous(expand = c(0, 0), labels = pretty_zero) +
    scale_x_continuous(
      breaks = sample_pos$x_center,
      labels = rep("", nrow(sample_pos)),
      limits = axis_limits,
      expand = c(0, 0)
    ) +
    labs(x = "", y = "FST distance") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = unit(c(0.2, -0.3, 0.2, 0.2), "cm")
    )

  title <- ggdraw() + draw_label(title_text, fontface = "bold")

  comb <- plot_grid(
    plt_dendr, plt_pi, plt_temp, plt_map,
    axis = "rl", align = "v",
    ncol = 1, rel_heights = c(0.2, 0.2, 0.3, 0.5)
  )

  plot_grid(title, comb, ncol = 1, align = "v", rel_heights = c(0.05, 1))
}

# -----------------------------
# Select samples for each bin (species)
# -----------------------------
sampfil <- switch(
  binname,
  "MAST4A" = c("TA_SUR_GGMM_102","TA_SUR_GGMM_11","TA_SUR_GGMM_111","TA_SUR_GGMM_112","TA_SUR_GGMM_124","TA_SUR_GGMM_131","TA_SUR_GGMM_132","TA_SUR_GGMM_135","TA_SUR_GGMM_136","TA_SUR_GGMM_137","TA_SUR_GGMM_138","TA_SUR_GGMM_142","TA_SUR_GGMM_143","TA_SUR_GGMM_145","TA_SUR_GGMM_146","TA_SUR_GGMM_147","TA_SUR_GGMM_148","TA_SUR_GGMM_149","TA_SUR_GGMM_150","TA_SUR_GGMM_151","TA_SUR_GGMM_152","TA_SUR_GGMM_16","TA_SUR_GGMM_18","TA_SUR_GGMM_20","TA_SUR_GGMM_22","TA_SUR_GGMM_23","TA_SUR_GGMM_25","TA_SUR_GGMM_30","TA_SUR_GGMM_32","TA_SUR_GGMM_4","TA_SUR_GGMM_5","TA_SUR_GGMM_51","TA_SUR_GGMM_6","TA_SUR_GGMM_64","TA_SUR_GGMM_65","TA_SUR_GGMM_66","TA_SUR_GGMM_68","TA_SUR_GGMM_7","TA_SUR_GGMM_72","TA_SUR_GGMM_76","TA_SUR_GGMM_78","TA_SUR_GGMM_80","TA_SUR_GGMM_81","TA_SUR_GGMM_83","TA_SUR_GGMM_9","TA_SUR_GGMM_93","TA_SUR_GGMM_95","TA_SUR_GGMM_96","TA_SUR_GGMM_97","TA_SUR_GGMM_98"),
  "MAST4B" = c("TA_SUR_GGMM_124","TA_SUR_GGMM_130","TA_SUR_GGMM_131","TA_SUR_GGMM_132","TA_SUR_GGMM_138","TA_SUR_GGMM_41","TA_SUR_GGMM_42","TA_SUR_GGMM_43","TA_SUR_GGMM_45","TA_SUR_GGMM_46","TA_SUR_GGMM_51"),
  "MAST4C" = c("TA_SUR_GGMM_100","TA_SUR_GGMM_102","TA_SUR_GGMM_106","TA_SUR_GGMM_109","TA_SUR_GGMM_11","TA_SUR_GGMM_110","TA_SUR_GGMM_113","TA_SUR_GGMM_123","TA_SUR_GGMM_124","TA_SUR_GGMM_125","TA_SUR_GGMM_128","TA_SUR_GGMM_129","TA_SUR_GGMM_130","TA_SUR_GGMM_131","TA_SUR_GGMM_132","TA_SUR_GGMM_136","TA_SUR_GGMM_137","TA_SUR_GGMM_138","TA_SUR_GGMM_139","TA_SUR_GGMM_142","TA_SUR_GGMM_143","TA_SUR_GGMM_32","TA_SUR_GGMM_34","TA_SUR_GGMM_36","TA_SUR_GGMM_38","TA_SUR_GGMM_39","TA_SUR_GGMM_41","TA_SUR_GGMM_42","TA_SUR_GGMM_43","TA_SUR_GGMM_45","TA_SUR_GGMM_46","TA_SUR_GGMM_5","TA_SUR_GGMM_51","TA_SUR_GGMM_52","TA_SUR_GGMM_58","TA_SUR_GGMM_64","TA_SUR_GGMM_65","TA_SUR_GGMM_7","TA_SUR_GGMM_72","TA_SUR_GGMM_9"),
  # default: MAST4E
  c("TA_SUR_GGMM_135","TA_SUR_GGMM_145","TA_SUR_GGMM_146","TA_SUR_GGMM_147","TA_SUR_GGMM_148","TA_SUR_GGMM_149","TA_SUR_GGMM_150","TA_SUR_GGMM_151","TA_SUR_GGMM_152","TA_SUR_GGMM_6","TA_SUR_GGMM_66","TA_SUR_GGMM_81","TA_SUR_GGMM_82","TA_SUR_GGMM_83","TA_SUR_GGMM_89","TA_SUR_GGMM_93")
)

# -----------------------------
# Load metadata and auxiliary data
# -----------------------------

# Genome sizes for normalization
bsize <- read.table(file = binsize, header = TRUE, sep = "\t") %>%
  as_tibble() %>%
  arrange(MASTs)

# Tara mapping read counts (used for RPKG-style normalization)
tara_reads <- read.table(file = taracounts, header = TRUE, sep = "\t")

# Intra-diversity file
intrapi <- read.table(file = intradiv, header = TRUE, sep = "\t") %>%
  select(Sample, Intra_pi) %>%
  mutate(Name = gsub("TA.SUR.GGMM.", "TA_SUR_GGMM_", Sample)) %>%
  select(Name, Intra_pi)

# Tara metadata
tara_meta <- read.table(file = tarameta, header = TRUE)

# -----------------------------
# Read BIN counts and normalize abundance
# -----------------------------

# Merge all count files matching pattern (one table)
cb_counts <- multmerge(bincounts, "MAST.*counts") %>%
  filter(str_detect(Samples, "TA")) %>%
  mutate(Samples = gsub("TA_SUR_GGMM_", "", Samples))

# Join with Tara read counts
cbn_counts <- merge(cb_counts, tara_reads, by = "Samples")

# Normalize counts (kept as in your original code)
cbn_counts[, 2:5] <- (cbn_counts[, 2:5] * 1000) / bsize$basepairs[col(cbn_counts[, 2:5])]
cbn_counts[, 2:5] <- (cbn_counts[, 2:5] * 10^9) / (cbn_counts[, 6][row(cbn_counts[, 2:5])] * 100)

# Extract abundance of the target bin
bin_counts <- cbn_counts %>%
  as_tibble() %>%
  select(Samples, all_of(binname))

colnames(bin_counts) <- c("Samples", "Abundance")

# Merge metadata + intrapi + abundance, then subset to the station list
bbin <- tara_meta %>%
  left_join(intrapi, by = "Name") %>%
  left_join(bin_counts, by = "Samples") %>%
  filter(Name %in% sampfil) %>%
  transmute(
    Samples,
    Latitude, Longitude,
    Temperature, Density, Distance_coast, Salinity,
    Abundance,
    NO2, PO4, Si,
    Intra_pi
  )

# -----------------------------
# SNP count per sample (from BCF stats file)
# -----------------------------
tb_snp <- read.table(file = bcf_input, header = FALSE, sep = "\t") %>%
  mutate(SNP = V12 + V13) %>%
  select(V3, SNP) %>%
  mutate(V3 = gsub("\\.", "_", V3))

colnames(tb_snp) <- c("Name", "SNP_count")

snpcount <- merge(tara_meta, tb_snp, by = "Name") %>%
  mutate(Samples = labels) %>%
  select(SNP_count, Samples)

bbin <- merge(bbin, snpcount, by = "Samples")

# -----------------------------
# Read and clean FST matrix
# -----------------------------
fst_raw <- read.table(file = fst_input, header = TRUE, stringsAsFactors = TRUE)

# Subset rows/cols to the desired sample set (matrix uses dots instead of underscores)
sampfil_mat <- gsub("_", "\\.", sampfil)

fst_sub <- fst_raw %>%
  filter(rownames(fst_raw) %in% sampfil_mat) %>%
  select(all_of(sampfil_mat))

fst_sub <- fst_sub[sampfil_mat, ]

# Convert to Tara 'labels' for plotting (TAxx_?? style)
fst_names <- data.frame(Name = gsub("\\.", "_", rownames(fst_sub)))

fst_link <- merge(fst_names, tara_meta, by = "Name")

rownames(fst_sub) <- fst_link[match(fst_names$Name, fst_link$Name), ]$labels
colnames(fst_sub) <- fst_link[match(fst_names$Name, fst_link$Name), ]$labels

# Replace negatives and missing values
fst_sub[fst_sub < 0] <- 0
fst_sub[is.na(fst_sub)] <- 0

# -----------------------------
# Dendrogram (UPGMA / average linkage)
# -----------------------------
fst_dist <- as.dist(fst_sub)
fst_dist <- zoo::na.locf(fst_dist)  # kept from your original script intent

dend_upgma <- as.dendrogram(hclust(fst_dist, method = "average"))

# Cluster assignment (cut at h = 0.15)
cluster_data <- data.frame(Cluster_ID = cutree(dend_upgma, h = 0.15)) %>%
  rownames_to_column("Samples")

plot_data <- merge(bbin, cluster_data, by = "Samples")

# Plot (optionally save with pdf())
plot_dendro_panels(df = plot_data, dend = dend_upgma, title_text = binname)

# -----------------------------
# PERMANOVA (adonis2) on FST distances
# -----------------------------

# Keep only samples that are present in FST (and drop the problematic one)
fst_mat <- fst_sub
keep_samples <- setdiff(rownames(fst_mat), sample_drop)

fst_mat <- fst_mat[keep_samples, keep_samples]
fst_dist2 <- as.dist(fst_mat)

meta_perm <- bbin %>%
  filter(Samples %in% keep_samples) %>%
  select(Samples, Latitude, Longitude, Temperature, Density, Distance_coast, Salinity, Abundance) %>%
  column_to_rownames("Samples")

# Z-score transform predictors (as you did)
meta_z <- scale(data.matrix(meta_perm), center = TRUE, scale = TRUE) %>%
  as.data.frame()

permanova_mast <- adonis2(fst_dist2 ~ Temperature + Salinity + Density, data = meta_z)
print(permanova_mast)

# Optional: exclude all "_MS" samples (your second PERMANOVA block)
fst_noMS <- fst_sub %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  filter(!str_detect(Sample, "_MS")) %>%
  select(Sample, !matches("_MS")) %>%
  column_to_rownames("Sample")

meta_noMS <- bbin %>%
  filter(Samples %in% rownames(fst_noMS)) %>%
  select(Samples, Latitude, Longitude, Temperature, Density, Distance_coast, Salinity, Abundance) %>%
  column_to_rownames("Samples")

meta_noMS_z <- scale(data.matrix(meta_noMS), center = TRUE, scale = TRUE) %>%
  as.data.frame()

permanova_mast_noMS <- adonis2(as.dist(fst_noMS) ~ Temperature + Salinity + Density, data = meta_noMS_z)
print(permanova_mast_noMS)
