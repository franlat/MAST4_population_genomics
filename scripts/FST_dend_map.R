# ============================================================
# FST dendrogram + environmental variables + abundance map
#
# Steps:
# 1) Read CoverM output: RPKM and covered bases, compute horizontal coverage
# 2) Filter samples by minimum horizontal coverage
# 3) Read Tara metadata + intradiversity (pi) + FST matrix
# 4) Build UPGMA dendrogram from FST distances and assign station clusters
# 5) Plot dendrogram + intradiversity + temperature + world map
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggdendro)
  library(dendextend)
  library(cowplot)
  library(maps)
  library(zoo)       # zoo::na.locf
})

# -----------------------------
# Configuration
# -----------------------------
genome <- "MAST4C"

coverage_file   <- "./Biocluster/SAGs/POPGEN/FST_plot/MAST4.TARA.coverage.95id.80rcov.tsv"
fst_file        <- "./Biocluster/PGs/SNPs/POGENOM/MAST4/TA.MAST4C.noN.c10.s4.fst.txt"
intradiv_file   <- "./Biocluster/PGs/SNPs/POGENOM/MAST4/TA.MAST4C.noN.c10.s4.intradiv.txt"
tara_meta_file  <- "./Biocluster/SAGs/POPGEN/FST_plot/ggmm.meta.v2.txt"
gsize_file      <- "./Biocluster/SAGs/POPGEN/FST_plot/mast.gsize"

min_hcov <- 0.25       # minimum horizontal coverage threshold
fst_cut_h <- 0.15      # dendrogram cut height for station clusters

# Optional: if you prefer working from a fixed directory, uncomment:
# setwd("/Volumes/My Passport for Mac/Back-up")

# -----------------------------
# Helper functions
# -----------------------------

pretty_zero <- function(x) {
  max_dec <- max(nchar(str_extract(x, "\\.[0-9]+")), na.rm = TRUE) - 1
  formatC(x, replace.zero = TRUE, zero.print = "0",
          digits = max_dec, format = "f", preserve.width = TRUE)
}

read_genome_size <- function(gsize_path, genome_id) {
  gsize_tbl <- read.table(file = gsize_path, header = TRUE, sep = "\t") %>%
    as_tibble() %>%
    arrange(MASTs) %>%
    filter(MASTs == genome_id)

  as.numeric(gsize_tbl[[2]][1])
}

read_coverm_rpkm <- function(coverage_path, genome_id) {
  read.table(file = coverage_path, header = TRUE, sep = "\t") %>%
    select(Genome, contains("RPKM")) %>%
    column_to_rownames("Genome") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Samples") %>%
    mutate(Samples = gsub(".SUR.GGMM.1.clean.fastq.gz.RPKM", "", Samples)) %>%
    select(Samples, all_of(genome_id)) %>%
    rename(Abundance = all_of(genome_id))
}

read_coverm_hcov <- function(coverage_path, genome_id, genome_size) {
  read.table(file = coverage_path, header = TRUE, sep = "\t") %>%
    select(Genome, contains("Covered.Bases")) %>%
    column_to_rownames("Genome") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Samples") %>%
    mutate(Samples = gsub(".SUR.GGMM.1.clean.fastq.gz.Covered.Bases", "", Samples)) %>%
    select(Samples, all_of(genome_id)) %>%
    mutate(
      Covered_Bases = as.numeric(.data[[genome_id]]),
      hcov = Covered_Bases / as.numeric(genome_size)
    ) %>%
    select(Samples, hcov)
}

read_tara_meta <- function(meta_path) {
  read.table(file = meta_path, header = TRUE) %>%
    as_tibble() %>%
    mutate(Samples = gsub("_SUR_GGMM_", "RA", Name)) %>%
    select(-Name)
}

read_intradiv <- function(intradiv_path) {
  read.table(file = intradiv_path, header = TRUE, sep = "\t") %>%
    as_tibble() %>%
    select(Sample, Intra_pi) %>%
    mutate(
      Samples  = gsub(pattern = ".SUR.GGMM.", replacement = "RA", Sample),
      Intra_pi = as.numeric(Intra_pi) * 100000
    ) %>%
    select(Samples, Intra_pi)
}

read_fst_subset <- function(fst_path, samples_ra, tara_meta) {
  fst_g <- read.table(file = fst_path, header = TRUE, stringsAsFactors = TRUE)

  # Harmonize IDs to "RAxx" convention
  colnames(fst_g) <- gsub(pattern = ".SUR.GGMM.", replacement = "RA", colnames(fst_g))
  rownames(fst_g) <- gsub(pattern = ".SUR.GGMM.", replacement = "RA", rownames(fst_g))

  fst <- fst_g %>%
    filter(rownames(fst_g) %in% samples_ra) %>%
    select(all_of(samples_ra))

  # Map to Tara station labels for plot-friendly names
  fst_names <- data.frame(Samples = rownames(fst))
  fst_link  <- merge(fst_names, tara_meta, by = "Samples")

  rownames(fst) <- fst_link[match(fst_names$Samples, fst_link$Samples), ]$labels
  colnames(fst) <- fst_link[match(fst_names$Samples, fst_link$Samples), ]$labels

  fst[fst < 0] <- 0
  fst[is.na(fst)] <- 0

  fst
}

build_fst_dendrogram <- function(fst_mat) {
  dm <- zoo::na.locf(as.dist(fst_mat))
  as.dendrogram(hclust(dm, method = "average"))
}

plot_dendro_panels <- function(df, dend, title_text, fst_cut_h) {

  dend_data <- dendro_data(dend, type = "rectangle")
  seg <- with(segment(dend_data),
              data.frame(x = x, y = y, xend = xend, yend = yend))

  sample_pos <- with(dend_data$labels,
                     data.frame(x_center = x, Samples = as.character(label), height = 0.5))

  ptdf <- df %>%
    left_join(sample_pos, by = "Samples") %>%
    mutate(
      ymin = x_center - height / 2,
      ymax = x_center + height / 2
    )

  axis_limits <- with(sample_pos,
                      c(min(x_center - 0.5 * height), max(x_center + 0.5 * height))) +
    0.1 * c(-1, 1)

  # Dendrogram panel
  plt_dendr <- ggplot(seg) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_hline(yintercept = fst_cut_h, linetype = "dotted", colour = "orange") +
    scale_y_continuous(expand = c(0, 0), labels = pretty_zero) +
    scale_x_continuous(
      breaks = sample_pos$x_center,
      labels = rep("", nrow(sample_pos)),
      limits = axis_limits,
      expand = c(0, 0)
    ) +
    labs(x = "", y = "Average FST") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.text.y  = element_text(size = 4),
      axis.title.y = element_text(size = 4),
      plot.margin  = unit(c(0.2, 0, 0, 0.2), "cm")
    )

  # Intradiversity panel
  max_pi <- max(ptdf$Intra_pi, na.rm = TRUE)

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
    coord_cartesian(ylim = c(0, max_pi + 0.0001)) +
    theme_bw() +
    theme(
      plot.margin  = unit(c(0, 0, 0, 0.2), "cm"),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_text(size = 4),
      axis.title.y = element_text(size = 4),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    )

  # Temperature panel
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
      plot.margin  = unit(c(0, 0, 0, 0.2), "cm"),
      axis.title.y = element_text(size = 4),
      axis.text.y  = element_text(size = 4),
      axis.text.x  = element_text(angle = 90, vjust = 0.1, hjust = 0.25, size = 3),
      panel.grid.minor = element_blank()
    )

  # World map panel
  world <- map_data("world") %>%
    filter(long >= -165, long <= 90, lat >= -70, lat <= 50)

  plt_map <- ggplot(ptdf, aes(Longitude, Latitude)) +
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "white", fill = "lightgray", linewidth = 0.1
    ) +
    geom_point(fill = "black", size = 0.3, alpha = 1) +
    geom_point(aes(size = Abundance, colour = as.factor(Cluster_ID)), alpha = 0.3) +
    guides(
      colour = "none",
      size = guide_legend(direction = "horizontal")
    ) +
    theme_bw() +
    theme(
      axis.ticks = element_blank(),
      axis.text  = element_blank(),
      legend.position = c(0.75, 0.12),
      legend.background = element_rect(fill = "white", linewidth = 0.2, colour = "grey"),
      legend.title = element_text(size = 6),
      legend.text  = element_text(size = 6),
      plot.margin  = unit(c(0, 0, 0, 0.2), "cm")
    ) +
    labs(x = "", y = "")

  # Combine panels
  title <- ggdraw() + draw_label(title_text, fontface = "bold")
  comb <- plot_grid(
    plt_dendr, plt_pi, plt_temp, plt_map,
    axis = "rl", align = "v",
    ncol = 1,
    rel_heights = c(0.2, 0.2, 0.3, 0.5)
  )

  plot_grid(title, comb, ncol = 1, align = "v", rel_heights = c(0.05, 1))
}

# -----------------------------
# Run analysis
# -----------------------------

# Genome size
gsize <- read_genome_size(gsize_file, genome)

# Abundance (RPKM) and horizontal coverage
rpkm <- read_coverm_rpkm(coverage_file, genome)
hcov <- read_coverm_hcov(coverage_file, genome, gsize)

# Filter samples by horizontal coverage
sampfil <- hcov %>%
  filter(hcov >= min_hcov) %>%
  pull(Samples)

# Metadata and intradiversity
tara_meta <- read_tara_meta(tara_meta_file)
intrapi   <- read_intradiv(intradiv_file)

# FST matrix subset to filtered samples
fst <- read_fst_subset(fst_file, sampfil, tara_meta)

# Dendrogram and station clusters
dendUPGMA <- build_fst_dendrogram(fst)

cluster_data <- data.frame(Cluster_ID = cutree(dendUPGMA, h = fst_cut_h)) %>%
  rownames_to_column("labels") %>%
  rename(Samples = labels)

# Merge all information into one plotting table
plot_data <- rpkm %>%
  inner_join(tara_meta, by = "Samples") %>%
  inner_join(cluster_data, by = "Samples") %>%
  inner_join(intrapi, by = "Samples") %>%
  mutate(Samples = factor(labels, levels = labels(dendUPGMA)))

# Plot (single panel set for the selected genome)
final_plot <- plot_dendro_panels(
  df = plot_data,
  dend = dendUPGMA,
  title_text = genome,
  fst_cut_h = fst_cut_h
)

final_plot
