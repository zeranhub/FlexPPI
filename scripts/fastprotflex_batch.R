#!/usr/bin/env Rscript

# Batch runner for the published FastProtFlex all-at-once predictor.
#
# The feature construction and coefficients follow predict.r and
# FUNCTION_GDV.r in https://github.com/jure-praznikar/FastProtFlex.  The
# original code uses netdist::count_orbits_per_node; here orca::count4 is used
# directly.  Both expose the ORCA node-orbit columns 0--14, while the latter is
# available as a compiled CRAN package and is practical for repeatable batches.

suppressPackageStartupMessages(library(bio3d))
suppressPackageStartupMessages(library(orca))

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3) {
  stop(
    "Usage: fastprotflex_batch.R MANIFEST.tsv PREDICTIONS.csv TIMINGS.csv",
    call. = FALSE
  )
}

manifest_path <- arguments[[1]]
prediction_path <- arguments[[2]]
timing_path <- arguments[[3]]

beta <- c(
  -1.9063749, -0.9075780, 2.9751331, 2.9623554, -0.9232815,
  -0.1610720, 0.1171952, -0.4793347, 0.2411250, 0.3110215,
  0.4226776, -0.7755721, 0.5147505, -1.7265063, -1.4770230
)

manifest <- read.delim(
  manifest_path,
  header = TRUE,
  sep = "\t",
  quote = "\"",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_columns <- c("sample_id", "pdb_path")
missing_columns <- setdiff(required_columns, colnames(manifest))
if (length(missing_columns) > 0) {
  stop(paste("Manifest lacks columns:", paste(missing_columns, collapse = ", ")))
}

predict_one <- function(sample_id, pdb_path, chain_filter = "") {
  total_started <- proc.time()[["elapsed"]]
  pdb <- read.pdb(pdb_path, verbose = FALSE)
  pdb <- trim.pdb(pdb, atom.select(pdb, "noh"))
  pdb <- trim.pdb(pdb, atom.select(pdb, "protein"))
  if (!is.na(chain_filter) && nzchar(chain_filter)) {
    requested_chains <- strsplit(chain_filter, ",", fixed = TRUE)[[1]]
    requested_chains[requested_chains == "_"] <- ""
    atom_chains <- as.character(pdb$atom$chain)
    atom_chains[is.na(atom_chains)] <- ""
    selected_atoms <- which(atom_chains %in% requested_chains)
    if (length(selected_atoms) == 0) {
      stop(paste(sample_id, "has no atoms for chain filter", chain_filter))
    }
    selected_xyz <- as.vector(rbind(
      3L * selected_atoms - 2L,
      3L * selected_atoms - 1L,
      3L * selected_atoms
    ))
    pdb <- trim.pdb(
      pdb,
      inds = list(atom = selected_atoms, xyz = selected_xyz),
      sse = FALSE
    )
  }
  atom_count <- nrow(pdb$atom)
  if (atom_count > 10000) {
    stop(paste(
      sample_id,
      "has more than 10,000 protein heavy atoms; this batch runner deliberately",
      "does not substitute a different by-parts approximation"
    ))
  }
  if (atom_count < 2) {
    stop(paste(sample_id, "contains fewer than two protein heavy atoms"))
  }

  prediction_started <- proc.time()[["elapsed"]]
  xyz <- matrix(pdb$xyz, ncol = 3, byrow = TRUE)
  distances <- as.matrix(dist(xyz))
  edge_indices <- which(distances < 7.0 & upper.tri(distances), arr.ind = TRUE)
  if (nrow(edge_indices) == 0) {
    stop(paste(sample_id, "has no graph edges at the 7 Angstrom cutoff"))
  }
  edge_matrix <- matrix(as.integer(edge_indices), ncol = 2)
  degrees <- tabulate(as.vector(edge_matrix), nbins = atom_count)
  if (any(degrees == 0)) {
    stop(paste(sample_id, "contains isolated protein heavy atoms"))
  }
  graphlet_counts <- orca::count4(edge_matrix)
  if (nrow(graphlet_counts) != atom_count || ncol(graphlet_counts) != 15) {
    stop(paste(sample_id, "returned an unexpected ORCA matrix shape"))
  }

  covalent_adjacency <- distances <= 2.1
  diag(covalent_adjacency) <- FALSE
  smoothed_counts <- graphlet_counts * 0
  for (atom_index in seq_len(atom_count)) {
    neighbours <- which(covalent_adjacency[atom_index, ])
    if (length(neighbours) == 0) {
      smoothed_counts[atom_index, ] <- graphlet_counts[atom_index, ]
    } else if (length(neighbours) == 1) {
      smoothed_counts[atom_index, ] <- 0.5 * (
        graphlet_counts[atom_index, ] + graphlet_counts[neighbours, ]
      )
    } else {
      smoothed_counts[atom_index, ] <- 0.5 * (
        graphlet_counts[atom_index, ] + colMeans(graphlet_counts[neighbours, , drop = FALSE])
      )
    }
  }

  features <- scale(log(smoothed_counts + 1))
  if (any(!is.finite(features))) {
    stop(paste(sample_id, "has a zero-variance or non-finite graphlet feature"))
  }
  raw_prediction <- as.numeric(features %*% beta)
  atom_scores <- as.numeric(scale(raw_prediction))
  prediction_seconds <- proc.time()[["elapsed"]] - prediction_started

  ca_indices <- which(pdb$calpha)
  ca_atoms <- pdb$atom[ca_indices, , drop = FALSE]
  insertion <- if ("insert" %in% colnames(ca_atoms)) ca_atoms$insert else ""
  insertion[is.na(insertion)] <- ""
  chain <- ca_atoms$chain
  chain[is.na(chain) | chain == ""] <- "_"
  predictions <- data.frame(
    sample_id = rep(sample_id, length(ca_indices)),
    residue_order = seq_along(ca_indices) - 1L,
    chain = chain,
    residue_number = as.character(ca_atoms$resno),
    insertion_code = as.character(insertion),
    amino_acid_3 = as.character(ca_atoms$resid),
    fastprotflex_score = atom_scores[ca_indices],
    stringsAsFactors = FALSE
  )
  timing <- data.frame(
    sample_id = sample_id,
    atom_count = atom_count,
    residue_count = length(ca_indices),
    edge_count = nrow(edge_matrix),
    prediction_seconds = prediction_seconds,
    end_to_end_seconds = proc.time()[["elapsed"]] - total_started,
    status = "ok",
    error = "",
    stringsAsFactors = FALSE
  )
  list(predictions = predictions, timing = timing)
}

prediction_tables <- list()
timing_tables <- list()
for (row_index in seq_len(nrow(manifest))) {
  sample_id <- as.character(manifest$sample_id[[row_index]])
  pdb_path <- as.character(manifest$pdb_path[[row_index]])
  chain_filter <- if ("chain_filter" %in% colnames(manifest)) {
    as.character(manifest$chain_filter[[row_index]])
  } else {
    ""
  }
  result <- tryCatch(
    predict_one(sample_id, pdb_path, chain_filter),
    error = function(condition) {
      list(
        predictions = NULL,
        timing = data.frame(
          sample_id = sample_id,
          atom_count = NA_integer_,
          residue_count = NA_integer_,
          edge_count = NA_integer_,
          prediction_seconds = NA_real_,
          end_to_end_seconds = NA_real_,
          status = "error",
          error = conditionMessage(condition),
          stringsAsFactors = FALSE
        )
      )
    }
  )
  if (!is.null(result$predictions)) {
    prediction_tables[[length(prediction_tables) + 1L]] <- result$predictions
  }
  timing_tables[[length(timing_tables) + 1L]] <- result$timing
  cat(sprintf(
    "\rFastProtFlex %d/%d: %s (%s)",
    row_index, nrow(manifest), sample_id, result$timing$status[[1]]
  ))
  flush.console()
}
cat("\n")

if (length(prediction_tables) > 0) {
  predictions <- do.call(rbind, prediction_tables)
} else {
  predictions <- data.frame(
    sample_id = character(), residue_order = integer(), chain = character(),
    residue_number = character(), insertion_code = character(),
    amino_acid_3 = character(), fastprotflex_score = double()
  )
}
timings <- do.call(rbind, timing_tables)
write.csv(predictions, prediction_path, row.names = FALSE, na = "")
write.csv(timings, timing_path, row.names = FALSE, na = "")

if (any(timings$status != "ok")) {
  warning(sprintf(
    "%d of %d structures failed; see the timing table",
    sum(timings$status != "ok"), nrow(timings)
  ))
}
