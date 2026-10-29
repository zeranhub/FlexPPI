# Observational wrapper: run the exact upstream script, then save unrounded output.
.libPaths(c("E:/Project/protein_flexibility_models/R_library",
            "E:/Project/protein_flexibility_models/r_libs", .libPaths()))
options(digits=17)
source("predict.r")
stopifnot(nrow(pdb$atom) == length(as.numeric(RMSFPR)))
stopifnot(all(is.finite(as.numeric(RMSFPR))))
if (nrow(pdb$atom) > 10000) {
  # Official chunking internally calls clean.pdb. Require atom order to stay exact.
  cleaned_for_validation <- suppressWarnings(clean.pdb(pdb, consecutive=TRUE, force.renumber=TRUE))
  stopifnot(nrow(cleaned_for_validation$atom) == nrow(pdb$atom))
  stopifnot(identical(as.numeric(cleaned_for_validation$xyz), as.numeric(pdb$xyz)))
  stopifnot(identical(as.character(cleaned_for_validation$atom$elety), as.character(pdb$atom$elety)))
  stopifnot(identical(as.character(cleaned_for_validation$atom$resid), as.character(pdb$atom$resid)))
}
saveRDS(list(atom=pdb$atom, prediction=as.numeric(RMSFPR)), "atom_predictions_exact.rds")
write.csv(data.frame(pdb$atom, prediction=as.numeric(RMSFPR)),
          "atom_predictions_full_precision.csv", row.names=FALSE, na="")
writeLines(capture.output(sessionInfo()), "R_session_info.txt")
