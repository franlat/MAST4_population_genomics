# ============================================================
# dN/dS heatmap (all genes) + gene clusters -> eggNOG pie charts
#
# Steps:
# 1) Load FST matrix and Tara metadata; build a station dendrogram and assign station clusters
# 2) Load dN/dS matrix (genes x samples), cap values for visualization, and plot a heatmap
# 3) Use k-means gene clusters from the heatmap to summarize eggNOG functional categories
# 4) Plot per-cluster eggNOG category composition as pie charts
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(pheatmap)
  library(dendextend)
  library(pals)    # cols25()
})

# -----------------------------
# Inputs
# -----------------------------
binname <- "MAST4E"

fst_input <- "Documents/Biocluster/PGs/SNPs/POGENOM/MAST4/TA.MAST4E.noN.c10.s4.fst.txt"
tarameta  <- "Documents/Biocluster/PGs/SNPs/profiles/ggmm.meta.v2.txt"
dnds_file <- "Documents/Biocluster/PGs/SNPs/SNPEFF/MAST4/dNdS/MAST4E.snpeff.dnds.txt"

gh_file     <- "Documents/Biocluster/SAGs/MAST4E/ANNOTATION_BBNORM/dbCAN/MAST4E_17SAGs_BBNorm_vs_dbCAN_FILTERED.tbl"
eggnog_file <- "Documents/Biocluster/SAGs/MAST4E/ANNOTATION_BBNORM/eggNOG/MAST4E_17SAGs_BBNorm_vs_eggNOG_FILTERED.tbl"

nog_annotations_tsv <- "~/Downloads/NOG.annotations.tsv"
nog_hierarchy_tsv   <- "~/Documents/Biocluster/SAGs/eggNOG_DB/eggNOG_hierarchy.txt"

# Parameters
fst_cut_h  <- 0.15   # height threshold to cluster stations from the FST dendrogram
k_gene     <- 50     # number of k-means gene clusters in the dN/dS heatmap
dnds_cap   <- 2      # cap dN/dS values at this threshold for visualization
heat_title <- "MAST-4E"

# -----------------------------
# Helper functions
# -----------------------------
to_matrix_ids <- function(x) gsub("_", "\\.", x)

cap_matrix <- function(m, cap = 2) {
  m[m > cap] <- cap
  m
}

# Convert station labels (e.g., TA11_MS) into the sample IDs used in the dN/dS table (e.g., TARA11)
tara_labels_from_labels <- function(labels_vec) {
  x <- str_split(labels_vec, "_", simplify = TRUE)[, 1]
  gsub("TA", "TARA", x)
}

# -----------------------------
# Sample lists per lineage
# -----------------------------
sampfil <- switch(
  binname,
  "MAST4A" = c("TA_SUR_GGMM_102","TA_SUR_GGMM_11","TA_SUR_GGMM_111","TA_SUR_GGMM_112","TA_SUR_GGMM_124","TA_SUR_GGMM_131","TA_SUR_GGMM_132","TA_SUR_GGMM_135","TA_SUR_GGMM_136","TA_SUR_GGMM_137","TA_SUR_GGMM_138","TA_SUR_GGMM_142","TA_SUR_GGMM_143","TA_SUR_GGMM_145","TA_SUR_GGMM_146","TA_SUR_GGMM_147","TA_SUR_GGMM_148","TA_SUR_GGMM_149","TA_SUR_GGMM_150","TA_SUR_GGMM_151","TA_SUR_GGMM_152","TA_SUR_GGMM_16","TA_SUR_GGMM_18","TA_SUR_GGMM_20","TA_SUR_GGMM_22","TA_SUR_GGMM_23","TA_SUR_GGMM_25","TA_SUR_GGMM_30","TA_SUR_GGMM_32","TA_SUR_GGMM_4","TA_SUR_GGMM_5","TA_SUR_GGMM_51","TA_SUR_GGMM_6","TA_SUR_GGMM_64","TA_SUR_GGMM_65","TA_SUR_GGMM_66","TA_SUR_GGMM_68","TA_SUR_GGMM_7","TA_SUR_GGMM_72","TA_SUR_GGMM_76","TA_SUR_GGMM_78","TA_SUR_GGMM_80","TA_SUR_GGMM_81","TA_SUR_GGMM_83","TA_SUR_GGMM_9","TA_SUR_GGMM_93","TA_SUR_GGMM_95","TA_SUR_GGMM_96","TA_SUR_GGMM_97","TA_SUR_GGMM_98"),
  "MAST4B" = c("TA_SUR_GGMM_124","TA_SUR_GGMM_130","TA_SUR_GGMM_131","TA_SUR_GGMM_132","TA_SUR_GGMM_138","TA_SUR_GGMM_41","TA_SUR_GGMM_42","TA_SUR_GGMM_43","TA_SUR_GGMM_45","TA_SUR_GGMM_46","TA_SUR_GGMM_51"),
  "MAST4C" = c("TA_SUR_GGMM_100","TA_SUR_GGMM_102","TA_SUR_GGMM_106","TA_SUR_GGMM_109","TA_SUR_GGMM_11","TA_SUR_GGMM_110","TA_SUR_GGMM_113","TA_SUR_GGMM_123","TA_SUR_GGMM_124","TA_SUR_GGMM_125","TA_SUR_GGMM_128","TA_SUR_GGMM_129","TA_SUR_GGMM_130","TA_SUR_GGMM_131","TA_SUR_GGMM_132","TA_SUR_GGMM_136","TA_SUR_GGMM_137","TA_SUR_GGMM_138","TA_SUR_GGMM_139","TA_SUR_GGMM_142","TA_SUR_GGMM_143","TA_SUR_GGMM_32","TA_SUR_GGMM_34","TA_SUR_GGMM_36","TA_SUR_GGMM_38","TA_SUR_GGMM_39","TA_SUR_GGMM_41","TA_SUR_GGMM_42","TA_SUR_GGMM_43","TA_SUR_GGMM_45","TA_SUR_GGMM_46","TA_SUR_GGMM_5","TA_SUR_GGMM_51","TA_SUR_GGMM_52","TA_SUR_GGMM_58","TA_SUR_GGMM_64","TA_SUR_GGMM_65","TA_SUR_GGMM_7","TA_SUR_GGMM_72","TA_SUR_GGMM_9"),
  c("TA_SUR_GGMM_135","TA_SUR_GGMM_145","TA_SUR_GGMM_146","TA_SUR_GGMM_147","TA_SUR_GGMM_148","TA_SUR_GGMM_149","TA_SUR_GGMM_150","TA_SUR_GGMM_151","TA_SUR_GGMM_152","TA_SUR_GGMM_6","TA_SUR_GGMM_66","TA_SUR_GGMM_81","TA_SUR_GGMM_82","TA_SUR_GGMM_83","TA_SUR_GGMM_89","TA_SUR_GGMM_93")
)

# -----------------------------
# Load tables
# -----------------------------
dnds <- read.table(dnds_file, header = TRUE)

tara_meta <- read.table(file = tarameta, header = TRUE) %>%
  mutate(Name = gsub("_", ".", Name))

gh_tax <- read.table(gh_file, header = FALSE) %>%
  transmute(
    Gene  = gsub(".t1", "", V3),
    GH_ID = gsub(".hmm", "", V1)
  )

eggnog_tax <- read.table(eggnog_file, header = FALSE) %>%
  transmute(
    Gene   = gsub(".t1", "", V3),
    NOG_ID = V1 %>%
      gsub("NOG\\.", "", .) %>%
      gsub("\\..*", "", .)
  )

all_tax <- full_join(gh_tax, eggnog_tax, by = "Gene")

# -----------------------------
# FST matrix -> dendrogram -> station clusters
# -----------------------------
fst_g <- read.table(file = fst_input, header = TRUE, stringsAsFactors = TRUE)

sampfil_mat <- to_matrix_ids(sampfil)

fst_o <- fst_g %>%
  filter(rownames(fst_g) %in% sampfil_mat) %>%
  select(all_of(sampfil_mat))

fst <- fst_o[sampfil_mat, ]

fst_names <- data.frame(Name = rownames(fst))
fst_link  <- merge(fst_names, tara_meta, by = "Name")

rownames(fst) <- fst_link[match(fst_names$Name, fst_link$Name), ]$labels
colnames(fst) <- fst_link[match(fst_names$Name, fst_link$Name), ]$labels

fst[fst < 0] <- 0
fst[is.na(fst)] <- 0

dm <- dendextend::na_locf(as.dist(fst))
dendUPGMA <- as.dendrogram(hclust(dm, method = "average"))

station_clusters <- data.frame(Cluster_ID = cutree(dendUPGMA, h = fst_cut_h)) %>%
  rownames_to_column("labels")

station_meta <- merge(station_clusters, tara_meta, by = "labels") %>%
  select(labels, Temperature, Latitude, Longitude, Cluster_ID)

# -----------------------------
# dN/dS heatmap (genes x samples)
# -----------------------------
# Column naming conventions can differ between tables; this block selects the
# sample set and then reorders columns based on station labels.

smp_ra <- gsub("_SUR_GGMM_", "RA", sampfil)

dnds_mat <- dnds %>%
  select(Gene, all_of(smp_ra)) %>%
  column_to_rownames("Gene") %>%
  as.matrix() %>%
  cap_matrix(dnds_cap)

tara_order <- tara_labels_from_labels(station_meta$labels)
dnds_mat <- dnds_mat[, tara_order, drop = FALSE]

breaks_list <- seq(0, dnds_cap, by = 0.05)
heat_cols <- colorRampPalette(c("#FFFFFF", "#F7EBEC", "#84C318", "#28AFB0", "#F6488D"))(length(breaks_list))

set.seed(42)
ph <- pheatmap(
  dnds_mat,
  color = heat_cols,
  breaks = breaks_list,
  main = heat_title,
  clustering_method = "average",
  clustering_distance_rows = "manhattan",
  cluster_cols = FALSE,
  kmeans_k = k_gene,
  fontsize_row = 10,
  cex = 0.8,
  cellwidth = 8,
  cellheight = 8,
  border_color = "white"
)

# -----------------------------
# Gene cluster ordering for plotting (precomputed)
# -----------------------------
gclust_order_a <- c(15,34,36,11,14,32,3,19,28,44,2,45,4,35,22,23,25,47,42,6,30,50,10,12,5,26,29,40,8,9,17,33,37,16,31,18,21,38,39,41,27,24,20,48,7,43,13,49,1,46)
gclust_order_b <- c(26,14,39,19,23,17,11,49,20,28,33,9,22,1,15,34,4,36,27,7,25,35,41,8,43,16,37,48,24,40,29,21,50,3,45,44,5,13,2,12,10,38,30,42,46,18,47,32,6,31)
gclust_order_c <- c(35,5,31,42,28,8,29,10,48,26,32,9,18,41,27,2,19,40,24,14,22,4,25,37,1,16,21,13,43,17,38,12,45,11,20,30,49,47,7,44,6,39,34,33,50,15,46,36,3,23)
gclust_order_e <- c(7,5,20,37,12,36,48,49,6,3,23,17,10,42,40,41,45,22,35,9,1,50,44,19,46,28,34,4,21,26,8,30,14,16,13,29,18,31,33,38,11,43,32,15,24,25,27,2,39,47)

gclust_order <- gclust_order_e

genes_cluster <- data.frame(
  Gene = names(ph$kmeans[[1]]$cluster),
  GCluster = factor(ph$kmeans[[1]]$cluster, levels = gclust_order)
)

gene_cluster_tax <- merge(genes_cluster, all_tax, by = "Gene")

# -----------------------------
# eggNOG category mapping (NOG -> letters -> category name)
# -----------------------------
eggNOG_df <- read.table(
  file = nog_annotations_tsv,
  quote = "",
  sep = "\t",
  header = FALSE
) %>% as.data.frame()

colnames(eggNOG_df) <- c("NOG", "ID", "V1", "V2", "Letter", "Description")

eggNOG_df <- eggNOG_df %>%
  select(NOG_ID = ID, Letter)

NOG_Letter <- eggNOG_df %>%
  mutate(Letter = strsplit(as.character(Letter), "")) %>%
  unnest(Letter)

NOG_Categories <- read.table(
  file = nog_hierarchy_tsv,
  header = TRUE,
  quote = "",
  sep = "\t"
)

gene_cluster_tax_letter <- merge(gene_cluster_tax, NOG_Letter, by = "NOG_ID", all.x = TRUE)
gene_cluster_tax_NOG <- merge(gene_cluster_tax_letter, NOG_Categories, by = "Letter", all.x = TRUE)

# -----------------------------
# Pie charts: eggNOG categories per gene cluster
# -----------------------------
cluster_counts <- gene_cluster_tax_NOG %>%
  group_by(GCluster, Name) %>%
  count() %>%
  ungroup() %>%
  filter(!is.na(Name), Name != "NA")

cp <- coord_polar(theta = "y")
cp$is_free <- function() TRUE

ggplot(cluster_counts, aes(x = " ", y = n, group = Name, fill = Name)) +
  geom_bar(width = 0.5, stat = "identity", color = "white") +
  cp +
  facet_wrap(. ~ GCluster, ncol = 10, nrow = 5, scales = "free") +
  scale_fill_manual(values = cols25(23)) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(2, "mm"),
    legend.text = element_text(size = 5),
    aspect.ratio = 1,
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# -----------------------------
# Simple summary statistics
# -----------------------------
station_clusters_h010 <- data.frame(Cluster_ID = cutree(dendUPGMA, h = 0.10)) %>%
  rownames_to_column("labels")

station_meta_h010 <- merge(station_clusters_h010, tara_meta, by = "labels") %>%
  select(labels, Temperature, Latitude, Longitude, Cluster_ID, Salinity)

cluster_2 <- station_meta_h010 %>% filter(Cluster_ID == 2)

mean(cluster_2$Salinity, na.rm = TRUE)
sd(cluster_2$Temperature, na.rm = TRUE)
max(fst)
