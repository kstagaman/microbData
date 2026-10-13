test_that("ordinate() stores results unnamed and builds a matching Ordination.log", {
  mD <- make_small_mD()

  mD <- ordinate(mD, method = "PCoA")
  expect_length(mD@Ordinations, 1)
  expect_null(names(mD@Ordinations))
  expect_s3_class(mD@Ordinations[[1]], "pcoa")

  log <- mD@Other.data[["Ordination.log"]]
  expect_s3_class(log, "data.table")
  expect_equal(nrow(log), 1)
  expect_equal(log$Index, 1L)
  expect_equal(log$Method, "PCoA")
  expect_equal(log$Distance.metric, "Bray-Curtis")
  expect_true(is.na(log$Formula))
  expect_false(log$Optimized)
})

test_that("repeated ordinate() calls append to Ordinations and Ordination.log rather than overwrite", {
  mD <- make_small_mD()

  mD <- ordinate(mD, method = "PCoA")
  mD <- ordinate(mD, method = "PCA")

  expect_length(mD@Ordinations, 2)
  expect_null(names(mD@Ordinations))
  expect_s3_class(mD@Ordinations[[1]], "pcoa")
  expect_s3_class(mD@Ordinations[[2]], "prcomp")

  log <- mD@Other.data[["Ordination.log"]]
  expect_equal(nrow(log), 2)
  expect_equal(log$Index, 1:2)
  expect_equal(log$Method, c("PCoA", "PCA"))
  expect_equal(log$Distance.metric, c("Bray-Curtis", NA_character_))
  expect_true(all(is.na(log$Formula)))
})

test_that("ordinate() distinguishes multiple dbRDA formulas run on the same distance matrix", {
  mD <- make_small_mD()

  mD <- ordinate(mD, method = "dbRDA", formula = ~ Exposure)
  mD <- ordinate(mD, method = "dbRDA", formula = ~ 1)

  expect_length(mD@Ordinations, 2)
  expect_null(names(mD@Ordinations))

  log <- mD@Other.data[["Ordination.log"]]
  expect_equal(nrow(log), 2)
  expect_equal(log$Method, c("dbRDA", "dbRDA"))
  expect_equal(log$Distance.metric, c("Bray-Curtis", "Bray-Curtis"))
  expect_equal(log$Formula, c("~Exposure", "~1"))
  # Same method + same distance metric, but distinguishable by formula
  expect_false(log$Formula[1] == log$Formula[2])
})

test_that("ordinate() validates inputs and reports informative errors/warnings", {
  mD <- make_small_mD()

  expect_error(ordinate(mD, method = "dbRDA"), "formula")
  expect_error(ordinate(mD, method = "dbRDA", formula = "not a formula"), "formula")
  expect_error(ordinate(1:3, method = "PCoA"), "microbData")
  expect_error(ordinate(mD, method = "not.a.method"))

  expect_warning(ordinate(mD, method = "PCoA", optimize = TRUE), "does not support")
  expect_warning(ordinate(mD, method = "PCA", only.for = "Bray-Curtis"), "ignored")

  mD.no.dist <- make_small_mD(with.distances = FALSE)
  expect_error(ordinate(mD.no.dist, method = "PCoA"), "has none")
})

test_that("ordinate() with update.mD = FALSE returns results without modifying mD", {
  mD <- make_small_mD()

  result <- ordinate(mD, method = "PCoA", update.mD = FALSE)
  expect_type(result, "list")
  expect_named(result, c("Ordinations", "Log"))
  expect_length(result$Ordinations, 1)
  expect_null(names(result$Ordinations))
  expect_equal(nrow(result$Log), 1)
  expect_equal(result$Log$Method, "PCoA")

  expect_null(mD@Ordinations)
  expect_null(mD@Other.data[["Ordination.log"]])
})

test_that("ordinate() supports all four registered methods and records them correctly", {
  mD <- make_small_mD()

  mD <- ordinate(mD, method = "PCoA")
  # The small test dataset can yield a (near-)perfect NMDS fit; the resulting
  # "stress is (nearly) zero" warning is expected here and not a real problem.
  mD <- suppressWarnings(ordinate(mD, method = "NMDS", trymax = 2, trace = FALSE))
  mD <- ordinate(mD, method = "dbRDA", formula = ~ Exposure)
  mD <- ordinate(mD, method = "PCA")

  expect_length(mD@Ordinations, 4)
  expect_equal(
    vapply(mD@Ordinations, function(x) class(x)[1], character(1)),
    c("pcoa", "metaMDS", "capscale", "prcomp")
  )

  log <- mD@Other.data[["Ordination.log"]]
  expect_equal(log$Method, c("PCoA", "NMDS", "dbRDA", "PCA"))
  expect_equal(log$Distance.metric, c("Bray-Curtis", "Bray-Curtis", "Bray-Curtis", NA_character_))
  expect_equal(log$Formula, c(NA_character_, NA_character_, "~Exposure", NA_character_))
})


test_that("ordination.coords() returns combined, correctly-labeled coordinates for all four methods", {
  mD <- make_small_mD()
  mD <- ordinate(mD, method = "PCoA")
  mD <- suppressWarnings(ordinate(mD, method = "NMDS", trymax = 2, trace = FALSE))
  mD <- ordinate(mD, method = "dbRDA", formula = ~ Exposure)
  mD <- ordinate(mD, method = "PCA")

  coords <- ordination.coords(mD, combine = TRUE)
  expect_true(all(c("Axis.labs", "Samples") %in% names(coords)))
  expect_true(all(c("Axis1", "Axis2", "Ord.method", "Beta.metric") %in% names(coords$Samples)))
  expect_setequal(coords$Samples$Ord.method, c("PCoA", "NMDS", "dbRDA", "PCA"))

  n.samples <- nrow(get.metadata(mD))
  expect_equal(nrow(coords$Samples), n.samples * 4)

  # PCA does not come from a distance matrix
  pca.rows <- coords$Samples[Ord.method == "PCA"]
  expect_true(all(is.na(pca.rows$Beta.metric)))

  # distance-based methods are correctly tagged with their metric
  dist.based <- coords$Samples[Ord.method %in% c("PCoA", "NMDS", "dbRDA")]
  expect_true(all(dist.based$Beta.metric == "Bray-Curtis"))
})

test_that("ordination.coords() combine = FALSE returns one set of coordinates per ordination", {
  mD <- make_small_mD()
  mD <- ordinate(mD, method = "PCoA")
  mD <- ordinate(mD, method = "PCA")

  result <- ordination.coords(mD, combine = FALSE)
  expect_length(result, 2)
  expect_null(names(result))
  for (r in result) {
    expect_true("Samples" %in% names(r))
    expect_true(all(c("Axis1", "Axis2") %in% names(r$Samples)))
  }
})

test_that("ordination.coords() infers ordination type from class when no mD/log is supplied", {
  mD <- make_small_mD()
  pcoa.result <- ordinate(mD, method = "PCoA", update.mD = FALSE)$Ordinations[[1]]
  pca.result <- ordinate(mD, method = "PCA", update.mD = FALSE)$Ordinations[[1]]

  coords <- ordination.coords(
    ord.list = list(pcoa.result, pca.result),
    metadata = get.metadata(mD),
    combine = TRUE
  )
  expect_setequal(coords$Samples$Ord.method, c("PCoA", "PCA"))
  expect_true(all(is.na(coords$Samples[Ord.method == "PCA"]$Beta.metric)))
})

test_that("ordination.coords() errors on an unsupported ordination class", {
  expect_error(
    ordination.coords(
      ord.list = list(structure(list(), class = "lm")),
      metadata = data.table::data.table(Sample = "A01")
    ),
    "not supported"
  )
})
