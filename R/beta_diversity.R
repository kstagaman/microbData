#' @name beta.diversity
#' @title Estimate Beta-diversity
#' @description This function estimates beta-diversity (see `metrics`) from the abundances table in a given microbData object. It returns a distance matrix (or list of distance matrices).
#' @param mD required; microbData object from which data will be read.
#' @param metrics character; a vector of beta-diveristy metrics to be estimated from the Abundances table in the microbData object. Supported metrics include "Manhattan", "Euclidean", "Canberra", "Clark", "Bray-Curtis", "Kulczynski", "Jaccard", "Gower", "Alt Gower", "Morisita", "Horn", "Mountford", "Raup", "Binomial", "Chao", "Cao", "Mahalanobis", "Chisq", "Chord", "Aitchison", "Robust Aitchison", "Sorensen" (binary Bray-Curtis), and any UniFrac distance with an alpha (abundance weighting parameter) between 0 and 1, inclusive. Alpha can be specified by putting in the metric name with a space before "UniFrac", e.g., "0.2 UniFrac". "W UniFrac" (weighted UniFrac) is equivalent to "1 UniFrac" and "U UniFrac" (unweighted UniFrac) is equivalent to "0 UniFrac". The default metrics are "Sorensen", "Canberra", "Bray-Curtis", "W UniFrac", "0.5 UniFrac", and "U UniFrac".
#' @param presence.absence logical; should abundance data be calculated on presence/absence data rather than abundance counts? This only applies to metrics supported by \code{\link[vegan]{vegdist}}. `metric = "Bray-Curtis"` with `presence.absence = TRUE` is equivalent to `metric = "Sorensen"` (which will ignore the `presence.absence` argument). Default is FALSE.
#' @param ncores integer; number of cores to use. 1 will run beta-diversity estimation serially. >1 will run beta-diversity estimation in parallel; this will ultimately only use a number of cores equal to the number of diversity metrics supplied. Default is 1.
#' @param update.mD logical; should this function return a new microbData object with the beta-diversity distance matrices and `metrics` added to Other.data (TRUE) or just the list of the beta-diversity distance matrices (FALSE)? Default is TRUE.
#' @details The function utilizes \code{\link[vegan]{vegdist}} and/or \code{\link[GUniFrac]{GUniFrac}} to estimate pairwise beta-diversity between each sample in a `microbData` object.
#' @returns Either a `microbData` object with the results add to the Distance.matrices slot and the names of the metrics recorded in the Other.data slot, OR a named list of the distance matrices.
#' @seealso \code{\link[vegan]{vegdist}}, \code{\link[GUniFrac]{GUniFrac}}
#' #' @examples
#' data("mD_rar")
#' mD.new <- beta.diversity(mD.rar, metrics = c("Bray-Curtis", "Canberra"))
#' print(mD.new)
#' @export

beta.diversity <- function(
    mD,
    metrics = c("Sorensen", "Canberra", "Bray-Curtis", "W UniFrac", "0.5 UniFrac", "U UniFrac"),
    presence.absence = FALSE,
    ncores = 1,
    update.mD = TRUE
) {
  vegdist.metric.names <- data.table(
    My.name = c(
      "Manhattan", "Euclidean", "Canberra", "Clark", "Bray-Curtis", "Kulczynski", "Jaccard", "Gower", "Alt Gower", "Morisita", "Horn", "Mountford", "Raup", "Binomial", "Chao", "Cao", "Mahalanobis", "Chisq", "Chord", "Aitchison", "Robust Aitchison", "Sorensen"
    ),
    Vegdist.name = c(
      "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "aitchison", "robust.aitchison", "bray"
    )
  ) %>% setkey(My.name)
  if (any(!str_detect(metrics, "UniFrac"))) {
    if (!all(str_subset(metrics, "UniFrac", negate = T) %in% vegdist.metric.names$My.name)) {
      rlang::abort(
        paste(
          "Arguments supplied to `metrics' must be of the following:",
          paste(vegdist.metric.names$My.name, collapse = ", "),
          "\nOR must include 'UniFrac' preceeded by U, W or a number between 0 and 1, inclusive"
        )
      )
    }
  }
  if (ncores < 1) {
    rlang::abort("Argument `ncore' must be an integer ≥ 1")
  }
  ncores <- as.integer(ncores)
  if (ncores == 1) {
    dist.list <- NULL
    for (metric in metrics) {
      if (metric %in% vegdist.metric.names$My.name) {
        orig.pa <- presence.absence
        if (metric == "Sorensen") {
          presence.absence <- TRUE
        }
        dist.list[[metric]] <- vegan::vegdist(
          mD@Abundances,
          method = vegdist.metric.names[metric]$Vegdist.name,
          binary = presence.absence
        )
        if (metric == "Sorensen") {
          presence.absence <- orig.pa
        }
      } else if (str_detect(metric, "UniFrac")) {
        alpha <- str_split(metric, " ") %>%
          sapply(head, 1) %>%
          sapply(function(a) {
            ifelse(a == "W", 1, ifelse(a == "U", 0, as.numeric(a)))
          })
        dist.list[[metric]] <- GUniFrac::GUniFrac(
          otu.tab = mD@Abundances,
          tree = mD@Phylogeny,
          alpha = alpha
        )$unifracs[, , 1] %>%
          as.dist()
      }
    }
  } else {
    require(foreach)
    require(doParallel)
    usable.cores <- min(c(ncores, length(metrics)))
    cl <- makeCluster(usable.cores, type = "FORK")
    registerDoParallel(cl, usable.cores)
    dist.list <- foreach::foreach(
      metric = metrics,
      .final = function(x) setNames(x, metrics)
    ) %dopar% {
      if (metric %in% vegdist.metric.names$My.name) {
        orig.pa <- presence.absence
        if (metric == "Sorensen") {
          presence.absence <- TRUE
        }
        dist.mat <- vegan::vegdist(
          mD@Abundances,
          method = vegdist.metric.names[metric]$Vegdist.name,
          binary = presence.absence
        )
        if (metric == "Sorensen") {
          presence.absence <- orig.pa
        }
      } else if (str_detect(metric, "UniFrac")) {
        alpha <- str_split(metric, " ") %>%
          sapply(head, 1) %>%
          sapply(function(a) {
            ifelse(a == "W", 1, ifelse(a == "U", 0, as.numeric(a)))
          })
        dist.mat <- GUniFrac::GUniFrac(
          otu.tab = mD@Abundances,
          tree = mD@Phylogeny,
          alpha = alpha
        )$unifracs[, , 1] %>%
          as.dist()
      }
      return(dist.mat)
    } %>%
      try(silent = TRUE)
    stopCluster(cl)
    if ("try-error" %in% class(dist.list)) {
      cat(dist.list, sep = "\n")
      rlang::abort("Beta-diversity estimation unsuccessful")
    }
  }
  if (update.mD) {
    cat(names(dist.list), sep = "\n")
    mD@Other.data[["Beta.metrics"]] <- metrics
    mD@Distance.matrices <- dist.list
    return(mD)
  } else {
    return(dist.list)
  }
}
