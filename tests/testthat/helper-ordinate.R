# Build a small microbData object for ordinate()/ordination.coords() tests.
# Uses only the top `n.features` most abundant ASVs so that NMDS (metaMDS)
# and the other ordination methods run quickly in tests.
make_small_mD <- function(n.features = 30, with.distances = TRUE, metrics = "Bray-Curtis") {
  data("metadata_dt", package = "microbData", envir = environment())
  data("asv_mat", package = "microbData", envir = environment())
  data("taxonomy_dt", package = "microbData", envir = environment())
  data.table::setkey(metadata.dt, Sample)
  data.table::setkey(taxonomy.dt, ASV)

  top.features <- names(sort(colSums(asv.mat), decreasing = TRUE))[seq_len(n.features)]
  small.mat <- asv.mat[, top.features]

  mD <- microbData(
    metadata = metadata.dt,
    abundances = small.mat,
    assignments = taxonomy.dt[ASV %in% top.features]
  )
  if (with.distances) {
    # Subsetting to the top features for test speed can leave a sample with
    # zero total abundance among just those features; vegdist() warns about
    # this but it is expected and harmless for these tests.
    mD <- suppressWarnings(beta.diversity(mD, metrics = metrics, update.mD = TRUE))
  }
  mD
}
