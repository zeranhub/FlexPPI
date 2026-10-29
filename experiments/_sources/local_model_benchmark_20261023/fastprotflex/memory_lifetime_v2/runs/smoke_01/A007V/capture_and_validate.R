# Called only by the explicitly requested future validation run.
# The original wrapper and original predict.r are byte-identical isolated copies.
source("capture_full_precision.R")

# predict.r clears the global workspace, so obtain arguments after it returns.
validation_arguments <- commandArgs(trailingOnly=TRUE)
stopifnot(length(validation_arguments) == 1)
stopifnot(file.exists(validation_arguments[1]))
reference <- readRDS(validation_arguments[1])
candidate <- readRDS("atom_predictions_exact.rds")
stopifnot(is.data.frame(reference$atom), is.data.frame(candidate$atom))

same_columns <- identical(names(reference$atom), names(candidate$atom))
same_count <- nrow(reference$atom) == nrow(candidate$atom) &&
              length(reference$prediction) == length(candidate$prediction)
identity_columns <- setdiff(names(reference$atom), c("x", "y", "z", "b"))
all_atom_identities_identical <- same_columns && identical(
    reference$atom[,identity_columns,drop=FALSE],
    candidate$atom[,identity_columns,drop=FALSE], num.eq=FALSE)
all_atom_coordinates_identical <- same_columns && identical(
    reference$atom[,c("x","y","z"),drop=FALSE],
    candidate$atom[,c("x","y","z"),drop=FALSE], num.eq=FALSE)
all_atom_table_identical <- identical(reference$atom, candidate$atom, num.eq=FALSE)
all_predictions_finite <- all(is.finite(reference$prediction)) && all(is.finite(candidate$prediction))
prediction_vector_identical <- identical(reference$prediction, candidate$prediction, num.eq=FALSE)
maximum_absolute_difference <- if (same_count && all_predictions_finite) {
    max(abs(reference$prediction - candidate$prediction))
} else { NA_real_ }
passed <- same_columns && same_count && all_atom_identities_identical &&
          all_atom_coordinates_identical && all_atom_table_identical &&
          all_predictions_finite && prediction_vector_identical &&
          is.finite(maximum_absolute_difference) && maximum_absolute_difference == 0
validation_boolean <- function(value) if (isTRUE(value)) "true" else "false"
writeLines(c(
    "{",
    sprintf('  "passed": %s,', validation_boolean(passed)),
    sprintf('  "same_atom_count": %s,', validation_boolean(same_count)),
    sprintf('  "reference_atom_count": %d,', nrow(reference$atom)),
    sprintf('  "candidate_atom_count": %d,', nrow(candidate$atom)),
    sprintf('  "same_columns": %s,', validation_boolean(same_columns)),
    sprintf('  "all_atom_identities_identical": %s,', validation_boolean(all_atom_identities_identical)),
    sprintf('  "all_atom_coordinates_identical": %s,', validation_boolean(all_atom_coordinates_identical)),
    sprintf('  "all_atom_table_identical": %s,', validation_boolean(all_atom_table_identical)),
    sprintf('  "all_predictions_finite": %s,', validation_boolean(all_predictions_finite)),
    sprintf('  "all_atom_prediction_vector_identical_num_eq_FALSE": %s,', validation_boolean(prediction_vector_identical)),
    sprintf('  "maximum_absolute_prediction_difference": %s', if (is.finite(maximum_absolute_difference)) sprintf("%.17g",maximum_absolute_difference) else "null"),
    "}"
), "comparison.json", useBytes=TRUE)
stopifnot(isTRUE(passed))
