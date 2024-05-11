## This script used to covert ArchR object into anndata object
## 2024.4.30

# Function used to covert archr object into anndata object
#' @title ArchRToAnndata
#' @param object ArchR object
#' @param outpath output dir
#' @param assays vector of assay names to export
#' @param assay.types vector of the type of assays, mainly peak, gene, motif
#' @param markersDF a named list of marker df vars
#' @param groups vector of cell metadata colnames to export
#' @param reductions vector of reduction names to export
#' @param export_links whether to expport coaccess link informs
#' @param export_pwm whether to export pwm informas
#' @importFrom ArchR getCellColData getMatrixFromProject getPeakAnnotation getCoAccessibility getEmbedding
#' @importFrom SummarizedExperiment rowRanges rowData colData assay seqnames start end assays
#' @importFrom GenomicRanges GRanges
#' @importFrom Matrix t
#' @importFrom anndata AnnData write_h5ad
#' @importFrom jsonlite toJSON
#' @export
#'
ArchrToAnndata <- function(object,
                           outpath,
                           assays,
                           assay.types,
                           markersDF = NULL,
                           groups = NULL,
                           reductions = NULL,
                           export_links = FALSE,
                           export_pwm = FALSE) {
  # Deal with metadata and coordinates information
  cell_meta <- as.data.frame(ArchR::getCellColData(ArchRProj = object))
  if (!is.null(groups) && groups %in% colnames(cell_meta)) {
    cell_meta <- cell_meta[, groups]
  }


  # Gain embedding
  coords_list <- list()
  dr <- object@`embeddings`
  if (!is.null(reductions)) {
    ReducNames <- intersect(reductions, names(dr))
    if (length(ReducNames) == 0) {
      stop("the reduction name provided not in the object")
    }
  } else {
    ReducNames <- names(dr)
    message("Using all embeddings contained in the object: ", paste(ReducNames, collapse = ", "))
  }

  for (embedding in ReducNames) {
    emb <- ArchR::getEmbedding(ArchRProj = object, embedding = embedding, returnDF = TRUE)
    colnames(emb) <- c(sprintf("X_%s", embedding), sprintf("Y_%s", embedding))
    coords_list[[embedding]] <- as.matrix(emb)
  }


  # Export assays
  for (i in seq_along(assays)) {
    exp_sce <- ArchR::getMatrixFromProject(ArchRProj = object, useMatrix = assays[i], binarize = TRUE)
    assay_names <- names(SummarizedExperiment::assays(exp_sce))
    if (length(assay_names) > 1) {
      exp_counts <- SummarizedExperiment::assay(exp_sce, assay_names[1])
    } else {
      exp_counts <- SummarizedExperiment::assay(exp_sce, assay_names)
    }


    # Get feature names
    # judge the sce object type, mainly:"RangedSummarizedExperiment" and "SummarizedExperiment"
    if (inherits(exp_sce, "RangedSummarizedExperiment")) {
      rowdata <- SummarizedExperiment::rowRanges(exp_sce)
      feature_names <- paste(SummarizedExperiment::seqnames(rowdata), SummarizedExperiment::start(rowdata), SummarizedExperiment::end(rowdata), sep = "-")
      feature_meta <- as.data.frame(rowdata)
      feature_meta$`feature` <- feature_names
      feature_meta <- feature_meta[, c("feature", "seqnames", "start", "end", "width", "strand", "idx")]
      # rownames(feature_meta) <- feature_names
    } else if (inherits(exp_sce, "SummarizedExperiment")) {
      feature_names <- SummarizedExperiment::rowData(exp_sce)$`name`
      feature_meta <- data.frame(feature = feature_names)
      # rownames(feature_meta) <- feature_names
    } else {
      stop("The matrix SummarizedExperiment object is not valid.")
    }
    rownames(feature_meta) <- feature_meta$`feature`



    # Create AnnData object
    # exp_adata <- sprintf("%s_adata", assays[i])
    exp_adata <- anndata::AnnData(
      X = Matrix::t(exp_counts),
      obs = cell_meta,
      var = feature_meta
    )


    # Add motif pwm information
    if (export_pwm && assay.types[i] == "motif") {
      motif_obj <- ArchR::getPeakAnnotation(object, name = "Motif")
      motif_pwm_obj <- motif_obj[["motifs"]]
      motif_pwm_list <- lapply(names(motif_pwm_obj), function(x) {
        motif_pwm_obj[[x]]@`profileMatrix`
      })
      names(motif_pwm_list) <- names(motif_pwm_obj)
      motif_pwm_json <- jsonlite::toJSON(motif_pwm_list, auto_unbox = TRUE)
      exp_adata$uns[["motif_logo"]] <- motif_pwm_json
    }


    # Add peak gene links
    if (export_links && assay.types[i] == "peak") {
      co_acc <- data.frame(ArchR::getCoAccessibility(object))
      co_data <- data.frame(
        seqnames = co_acc$`seqnames`,
        start = co_acc$`start`,
        end = co_acc$`end`,
        width = co_acc$`width`,
        strand = co_acc$`strand`,
        score = co_acc$`value`,
        group = co_acc$`group`
      )
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
    #exp_adata$write_h5ad(file_path)
    anndata::write_h5ad(exp_adata, file_path)
  }
}
