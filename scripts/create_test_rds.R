#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(Seurat)
})

parse_args <- function(args) {
  result <- list(project_root = ".", config = "config/pipeline.json", studies = "all", cells = 200L, seed = 42L)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!key %in% c("--project-root", "--config", "--studies", "--cells", "--seed")) {
      stop("Usage: create_test_rds.R [--project-root DIR] [--config config/pipeline.json] [--studies all] [--cells 200] [--seed 42]")
    }
    if (i == length(args)) stop("Missing value for ", key)
    result[[sub("^--", "", key)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  result$cells <- as.integer(result$cells)
  result$seed <- as.integer(result$seed)
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.na(args$cells) || args$cells < 1L) stop("--cells must be positive")
if (is.na(args$seed)) stop("--seed must be an integer")

root <- normalizePath(args$project_root, mustWork = TRUE)
pipeline_path <- file.path(root, args$config)
if (!file.exists(pipeline_path)) stop("Configuration does not exist: ", pipeline_path)
pipeline <- fromJSON(pipeline_path, simplifyVector = FALSE)
studies_path <- file.path(root, pipeline$studies_config)
document <- fromJSON(studies_path, simplifyVector = FALSE)
requested <- if (args$studies == "all") names(document$studies) else strsplit(args$studies, ",", fixed = TRUE)[[1]]
requested <- trimws(requested)
unknown <- setdiff(requested, names(document$studies))
if (length(unknown)) stop("Unknown studies: ", paste(unknown, collapse = ", "))

for (study_id in requested) {
  study <- document$studies[[study_id]]
  conversion <- study$conversion
  if (is.null(conversion)) next
  if (is.null(conversion$source) || is.null(conversion$test_source)) {
    stop(study_id, ": conversion studies require both source and test_source")
  }
  input <- file.path(root, conversion$source)
  output <- file.path(root, conversion$test_source)
  if (!file.exists(input)) stop(study_id, ": input does not exist: ", input)

  object <- readRDS(input)
  available <- Cells(object)
  set.seed(args$seed)
  selected <- if (length(available) > args$cells) sample(available, args$cells) else available
  fixture <- subset(object, cells = selected)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  saveRDS(fixture, output)
  message("Wrote ", study_id, ": ", length(selected), " cells to ", output)
}
