#' @name ordination.coords
#' @title Get Plotting Coordinates from (an) Ordination(s)
#' @description Get the appropriate coordinates to plot the results of one or more ordinations. This also returns axis labels with percentage variance explained for relevant ordination methods (i.e., for all supported methods but NMDS).
#' @param mD microbData; a \code{microbData} object that includes ordination results (typically labelled "Ordinations" in the Other.data slot). Required if \code{ord.list} and \code{metadata} are NULL. Default is NULL.
#' @param ord.list.name character; required if \code{mD} has an list of ordinations stores under something other than Other.data$Ordinations. Default is "Ordinations".
#' @param ord.list list; a list of ordinations, required if \code{mD} is NULL. This list **must** be named if \code{combine} is TRUE. Default is NULL.
#' @param metadata array; a data.frame/data.table/matrix with metadata for the samples used in the ordination(s). Default is NULL.
#' @param feature.coords logical; if FALSE will only get coordinates for samples. If TRUE, will also get coordinates for features (e.g. ASVs, KOs, ...). Default is FALSE.
#' #' @param feature.tbl array; a data.frame/data.table/matrix with data for the features used in the ordination(s). This is only necessary if there is no \code{mD} and \code{feature.coords} is TRUE. Default is NULL.
#' @param constraint.coords logical; only relevant for dbRDA ordinations. If FALSE, will only get coordinates for samples. If TRUE, will also get coordinates for constraining variables to plot as vectors. Default is FALSE.
#' @param combine logical; if TRUE, will combine related results into data.tables and add columns of ordination names so results can all be plotted together, e.g. by using the \code{\link[ggplot2]{facet_wrap}} or \code{\link[ggplot2]{facet_grid}} functions with \code{\link[ggplot2]{ggplot}}. If FALSE, will return a list of data.tables with coordinates. Default is TRUE.
#' @seealso \code{\link[vegan]{scores}}, \code{\link[vegan]{eigenvals}}
#' @export

ordination.coords <- function(
    mD = NULL,
    ord.list.name = "Ordinations",
    ord.list = NULL,
    metadata = NULL,
    feature.coords = FALSE,
    feature.tbl = NULL,
    constraint.coords = FALSE,
    combine = TRUE,
    axes = 1:2,
    axis.digits = 2
) {
  if (is.null(mD) & is.null(ord.list)) {
    rlang::abort(
      "Either `mD' or `ord.list' must be supplied to get ordination coordinates."
    )
  }
  if (!is.null(mD)) {
    ord.list <- mD@Other.data[[ord.list.name]]
    metadata <- mD@Metadata
    smpl.col.name <- mD@Sample.col
  }
  if (!is.null(metadata)) {
    if (!{"data.table" %in% class(metadata)}) {
      metadata <- as.data.table(metadata)
    }
  }
  if (feature.coords) {
    if (!is.null(mD)) {
      feature.tbl <- mD@Features
    }
    if (!is.null(feature.tbl)) {
      if (!{"data.table" %in% class(feature.tbl)}) {
        feature.tbl <- as.data.table(feature.tbl)
      }
    }
  }

  ord.dts <- lapply(seq_along(ord.list), function(i) {
    sections <- NULL
    if (combine) {
      if (str_detect(names(ord.list)[i], "NMDS$|PCoA$|dbRDA$")) {
        dist.name <- str_split(names(ord.list)[i], "\\_") %>% unlist() %>% head(1)
        ord.type <- str_split(names(ord.list)[i], "\\_") %>% unlist() %>% tail(1)
      } else {
        ord.class <- class(ord.list[[i]])
        ord.type <- ifelse(
          "metaMDS" %in% ord.class, "NMDS",
          ifelse(
            "capscale" %in% ord.class, "dbRDA",
            ifelse("pcoa" == ord.class, "pcoa", NA)
          )
        )
        if (is.na(ord.type)) {
          rlang::abort(paste("The ordination supplied at element", i, "is not supported."))
        }
      }
    }
    ord <- ord.list[[i]]
    if (ord.type == "PCoA") {
      smpl.mat <- ord$vectors[, axes]
    } else {
      smpl.mat <- vegan::scores(ord, display = "sites", choices = axes)
    }
    orig.axis.names <- colnames(smpl.mat)[1:2]
    if (ord.type == "NMDS") {
      sections[["Axis.labs"]] <- setNames(orig.axis.names, orig.axis.names)
    } else {
      if (ord.type == "PCoA") {
        pct.explained <- setNames(ord$values$Relative_eig[axes], orig.axis.names)
      } else {
        pct.explained <- { vegan::eigenvals(ord) / sum(vegan::eigenvals(ord)) }[1:2]
      }
      sections[["Axis.labs"]] <- sapply(orig.axis.names, function(n) {
        paste0(n, " (", round(pct.explained[n], axis.digits), "%)")
      })
    }
    if (combine) {
      sections[["Axis.labs"]] <- as.data.table(list(Label = sections[["Axis.labs"]])) %>%
        .[, `:=`(Axis = c("Axis1", "Axis2"), Ord.method = ord.type, Beta.metric = dist.name)]
    }
    colnames(smpl.mat)[1:2] <- c("Axis1", "Axis2")
    if (!is.null(metadata)) {
      if ("sorted" %in% names(attributes(metadata))) {
        smpl.col.name <- attributes(metadata)$sorted
      } else {
        smpl.col.name <- names(metadata)[
          apply(metadata, MARGIN = 2, FUN = function(x) {all(row.names(smpl.mat) %in% x)})
        ]
        if (length(smpl.col.name) == 0) {
          rlang::abort(
            paste(
              "The supplied metadata does not contain all samples in the ordination at element", i, ", please review data and re-try."
            )
          )
        } else {
          setkeyv(metadata, smpl.col.name)
        }
      }
      sections[["Samples"]] <- as.data.table(
        smpl.mat,
        keep.rownames = smpl.col.name
      ) %>%
        setkeyv(smpl.col.name) %>%
        merge(metadata)
    } else {
      sections[["Samples"]] <- as.data.table(smpl.mat, keep.rownames = "Sample") %>%
        setkey(Sample)
    }

    if (combine) {
      sections[["Samples"]][, `:=`(Ord.method = ord.type, Beta.metric = dist.name)]
      label.pad <- 1.1
      sections[["Axis.labs"]][
        , `:=`(
          Axis1 = c(max(sections[["Samples"]]$Axis1) * label.pad, min(sections[["Samples"]]$Axis1) * label.pad),
          Axis2 = c(min(sections[["Samples"]]$Axis2) * label.pad, max(sections[["Samples"]]$Axis2) * label.pad),
          Angle = ifelse(Axis == "Axis1", 0, 90),
          Vjust = ifelse(Axis == "Axis1", 1, 0)
        )
      ]
    }

    if (feature.coords) {
      feat.mat <- vegan::scores(ord, display = "species")
      colnames(feat.mat)[1:2] <- c("Axis1", "Axis2")
      if (!is.null(feature.tbl)) {
        if ("sorted" %in% names(attributes(feature.tbl))) {
          feat.col.name <- attributes(feature.tbl)$sorted
        } else {
          feat.col.name <- names(feature.tbl)[
            apply(feature.tbl, MARGIN = 2, FUN = function(x) {
              all(row.names(feat.mat) %in% x)
            })
          ]
          if (length(feat.col.name) == 0) {
            rlang::abort(
              paste(
                "The supplied feature.tbl does not contain all samples in the ordination at element", i, ", please review data and re-try."
              )
            )
          } else {
            setkeyv(feature.tbl, feat.col.name)
          }
        }
        sections[["Features"]] <- as.data.table(feat.mat, keep.rownames = feat.col.name) %>%
          setkeyv(feat.col.name) %>%
          merge(feature.tbl)
      } else {
        sections[["Features"]] <- as.data.table(feat.mat, keep.rownames = "Feature") %>%
          setkey(Feature)
      }
      if (combine) {
        sections[["Features"]][, `:=`(Ord.method = ord.type, Beta.metric = dist.name)]
      }
    }
    if (ord.type == "dbRDA" & constraint.coords) {
      sections[["Vectors"]] <- scores(ord, display = "bp") %>%
        as.data.table(keep.rownames = "Old")
      sections[["Vectors"]][
        , Variable := str_extract(Old, paste(names(metadata), collapse = "|"))
      ]
      setcolorder(sections[["Vectors"]], ncol(sections[["Vectors"]]))
      replace.cols <- str_which(names(sections[["Vectors"]]), paste(orig.axis.names, collapse = "|"))
      names(sections[["Vectors"]])[replace.cols] <- c("Axis1", "Axis2")
      if (combine) {
        sections[["Vectors"]][, `:=`(Ord.method = ord.type, Beta.metric = dist.name)]
      }

      sections[["Centroids"]] <- scores(ord, display = "cn") %>%
        as.data.table(keep.rownames = "Class")
      sections[["Centroids"]][
        , Variable := str_extract(Class, paste(names(metadata), collapse = "|"))
      ]
      sections[["Centroids"]][, Class := str_remove(Class, paste0("^", Variable))]
      setcolorder(sections[["Centroids"]], ncol(sections[["Centroids"]]))
      replace.cols <- str_which(names(sections[["Centroids"]]), paste(orig.axis.names, collapse = "|"))
      names(sections[["Centroids"]])[replace.cols] <- c("Axis1", "Axis2")
      if (combine) {
        sections[["Centroids"]][, `:=`(Ord.method = ord.type, Beta.metric = dist.name)]
      }
    }
    return(sections)
  })
  if (combine) {
    section.names <- lapply(ord.dts, names) %>% unlist() %>% unique() %>% set_names()
    lapply(section.names, function(sn) {
      lapply(ord.dts, function(x) { x[[sn]] }) %>% rbindlist()
    }) %>% return()
  } else {
    names(ord.dts) <- names(ord.list)
    return(ord.dts)
  }
}
