## This script used to cover siganc object into anndata object
## 2024.4.30


# Function used to covert signac object into anndata object
#' @title SignacToAnndata
#' @name Function SignacToAnndata
#' @param object Siganc object
#' @param outpath output dir
#' @param assays vector of assay names to export
#' @param assay.types vector of the type of assays, mainly peak, gene, motif
#' @param markersDF a named list of marker df vars
#' @param groups vector of cell metadata colnames to export
#' @param reductions vector of reduction names to export
#' @param export_links whether to expport coaccess link informs
#' @param export_pwm whether to export pwm informas
#' @importFrom Signac GetMotifData Links
#' @importFrom Seurat GetAssayData DefaultAssay Cells Embeddings
#' @importFrom jsonlite toJSON
#' @importFrom anndata AnnData write_h5ad
#' @export
#'

SignacToAnndata <- function(object,
                            outpath,
                            assays,
                            assay.types,
                            markersDF = NULL,
                            groups = NULL,
                            reductions = NULL,
                            export_links = FALSE,
                            export_pwm = FALSE) {
  # Deal the meta and coords inform
  cell_meta <- object[[]]
  if (!is.null(groups) && groups %in% colnames(object[[]])) {
    cell_meta <- cell_meta[, which(colnames(cell_meta) %in% groups)]
  }

  ## Add embedding information
  dr <- object@`reductions`
  coords_list <- list()
  if (!is.null(reductions)) {
    ReducNames <- intersect(reductions, names(dr))
    if (length(ReducNames) == 0) {
      stop("the reduction name provided not in the object")
    } else {
      ReducNames <- names(dr)
      message("Using all embeddings contained in the object: ", paste(ReducNames, collapse = ", "))
    }
  }

  for (embedding in ReducNames) {
    emb <- Seurat::Embeddings(object = object, embedding)
    if (ncol(emb) > 2) {
      emb <- emb[, 1:2]
    }
    colnames(emb) <- c(sprintf("X_%s", embedding), sprintf("Y_%s", embedding))
    coords_list[[embedding]] <- as.matrix(emb)
  }


  # Export assays
  for (i in seq(length(assays))) {
    Seurat::DefaultAssay(object) <- assays[i]
    feature_meta <- data.frame(feature = rownames(object))
    rownames(feature_meta) <- feature_meta[["feature"]]

    ## Attension the exp data fetch method between seurat v3 ~ seurat v5
    if (!is.null(object@assays[[assays[i]]]@`data`)) {
      exp_adata <- anndata::AnnData(
        X = Matrix::t(Seurat::GetAssayData(object = object, assay = assays[i], slot = "data")),
        obs = cell_meta,
        var = feature_meta
      )
    } else {
      exp_adata <- anndata::AnnData(
        X = Matrix::t(Seurat::GetAssayData(object = object, assay = assays[i], slot = "counts")),
        obs = cell_meta,
        var = feature_meta
      )
    }


    # Add motif pwm information
    if (export_pwm && assay.types[i] == "motif") {
      peak_index <- which(assay.types == "peak")
      pwm_data <- Signac::GetMotifData(object = object, slot = "pwm", assay = assays[peak_index])
      motif_name <- Signac::GetMotifData(object = object, slot = "motif.names", assay = assays[peak_index])
      new_pwm_name <- as.vector(paste(names(pwm_data), motif_name, sep = "_"))
      names(pwm_data) <- new_pwm_name
      pwm_json <- jsonlite::toJSON(pwm_data, auto_unbox = TRUE)
      exp_adata$uns[["motif_logo"]] <- pwm_json
    }


    # Add peak gene links
    if (export_links && assay.types[i] == "peak") {
      co_data <- as.data.frame(Signac::Links(object = object[[assays[i]]]))
      exp_adata$uns[["coacc_link"]] <- co_data
    }


    # Add marker information
    if (!is.null(markersDF) && !is.null(markersDF[[assays[i]]])) {
      exp_adata$uns[["markers"]] <- markersDF[[assays[i]]]
    } else {
      message("no marker provided")
    }



    # Add coords informations
    for (eb in ReducNames) {
      exp_adata$obsm[[eb]] <- coords_list[[eb]]
    }

    # Create output dir
    if (!dir.exists(outpath)) {
      dir.create(outpath)
    }


    # Write h5ad
    file_path <- file.path(outpath, sprintf("%s.h5ad", assays[i]))
    exp_adata$write_h5ad(file_path)
  }
}
