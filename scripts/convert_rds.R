#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(sceasy)
})

parse_args <- function(args) {
  result <- list(assay = "RNA")
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!key %in% c("--input", "--output", "--assay")) {
      stop("Unknown argument: ", key)
    }
    if (i == length(args)) {
      stop("Missing value for ", key)
    }
    result[[sub("^--", "", key)]] <- args[[i + 1]]
    i <- i + 2
  }
  if (is.null(result$input) || is.null(result$output)) {
    stop("Usage: convert_rds.R --input INPUT.rds --output OUTPUT.h5ad [--assay RNA]")
  }
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!file.exists(args$input)) {
  stop("Input does not exist: ", args$input)
}
if (!grepl("\\.rds$", args$input, ignore.case = TRUE)) {
  stop("Input must have an .rds extension")
}
if (!grepl("\\.h5ad$", args$output, ignore.case = TRUE)) {
  stop("Output must have an .h5ad extension")
}

message("Reading ", args$input)
object <- readRDS(args$input)

# Seurat v5 uses layer=; older Seurat versions use slot=. Keep both paths so the
# converter can reproduce legacy objects without retaining normalized matrices.
counts <- tryCatch(
  GetAssayData(object, assay = args$assay, layer = "counts"),
  error = function(e) GetAssayData(object, assay = args$assay, slot = "counts")
)
if (nrow(counts) == 0 || ncol(counts) == 0) {
  stop("Selected count matrix is empty")
}

minimal <- CreateSeuratObject(counts = counts, meta.data = object[[]])
dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
message("Writing raw counts and observation metadata to ", args$output)
sceasy::convertFormat(
  minimal,
  from = "seurat",
  to = "anndata",
  outFile = args$output,
  assay = "RNA",
  main_layer = "counts"
)

