#!/usr/bin/env Rscript

library(argparse)
library(jsonlite)
library(zellkonverter)

parser <- ArgumentParser(description = "Benchmarking entrypoint")

# define arguments passed into Omnibenchmark
parser$add_argument(
  "--output_dir", "-o",
  dest = "output_dir", type = "character",
  help = "output directory where files will be saved",
  default = getwd(), required = TRUE
)
parser$add_argument(
  "--name", "-n",
  dest = "name", type = "character",
  help = "name of the module",
  required = TRUE
)
parser$add_argument(
  "--data.h5ad",
  dest = "data_path", type = "character",
  help = "input data h5ad path",
  required = TRUE
)
parser$add_argument(
  "--data.clusters_truth_num",
  dest = "clusters_truth_num_path", type = "character",
  help = "input true number of clusters path",
  required = TRUE
)
parser$add_argument(
  "--seed",
  dest = "seed", type = "integer",
  help = "seed for the random number generator",
  default = 123
)

args <- parser$parse_args()

# output paths
timings_path <- file.path(args$output_dir, paste0(args$name, ".timings.json"))
clustering_info_path <- file.path(args$output_dir, paste0(args$name, ".clustering_info.json"))
clusters_path <- file.path(args$output_dir, paste0(args$name, ".clusters.tsv"))
pca_path <- file.path(args$output_dir, paste0(args$name, ".pca.tsv"))

# read data
sce <- readH5AD(args$data_path, reader = "python")
cell_ids <- colnames(sce)
clusters_truth_num <- as.integer(readLines(args$clusters_truth_num_path, n = 1))
n_cluster <- clusters_truth_num

# time object to store time involved (in seconds) in each step
time <- list(
  gpu_load = NA_real_, find_mit_gene = NA_real_, filter = NA_real_, normalization = NA_real_,
  hvg = NA_real_, scaling = NA_real_, pca = NA_real_,
  t_sne = NA_real_, umap = NA_real_,
  louvain = NA_real_, leiden = NA_real_
)

clustering_info <- list(
  resolutions = list(louvain = NA_real_, leiden = NA_real_),
  num_runs = list(louvain = NA_integer_, leiden = NA_integer_)
)

# generate random clusters
set.seed(args$seed)
random_clusters <- sample(seq_len(n_cluster), size = length(cell_ids), replace = TRUE)

# since this method does not involve PCA,
# yet the metrics module in benchmark.yaml needs it as input for the metrics stage,
# create dummy PCA results matrix with 2 cols (arbitrary small number that looks like a table)
pca <- matrix(rnorm(length(cell_ids)*2), ncol = 2, dimnames = list(cell_ids, NULL))

# write outputs to files
write_json(
  time, timings_path,
  auto_unbox = TRUE, pretty = TRUE
)
write_json(
  clustering_info, clustering_info_path,
  auto_unbox = TRUE, pretty = TRUE
)

write.table(
  data.frame(
    cell_id = cell_ids,
    louvain = random_clusters, leiden = random_clusters # col names required by metrics module
  ),
  clusters_path,
  sep = "\t", quote = F, row.names = F
)
pca <- data.frame(
  cell_id = rownames(pca),
  pca
)
colnames(pca) <- c(
  "cell_id",
  paste0("PC", seq_len(ncol(pca) - 1))
)
write.table(
  pca, pca_path,
  sep = "\t", quote = F, row.names = F
)