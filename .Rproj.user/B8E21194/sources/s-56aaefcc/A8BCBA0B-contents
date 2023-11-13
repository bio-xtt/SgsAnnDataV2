## This script used to cover siganc object into anndata object
# 对于siganc 对象其motif logo信息只能存在peak adata对象里
## 2023.8.1
# require(Signac)
# require(Seurat)
# require(jsonlite)
# require(anndata)
# require(Matrix)
# require(parallel)
# require(SeuratObject)
# options(Seurat.object.assay.version = "v5")




# Function used to covert signac object into anndata object
#' @title SignacToAnndata
#' @name Function SignacToAnndata
#' @param object siganc object
#' @param outpath anndata path of the output dir
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
# @useDynLib sgsAnnData2
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

  ############ deal the meta and coords inform
  if (!is.null(groups) && groups %in% colnames(object[[]])) {
    cell_meta <- object[[]]
    cell_meta <- cell_meta[, which(colnames(cell_meta) %in% groups)]
  } else {
    cell_meta <- object[[]]
  }

  ##### add embedding information
  dr <- object@`reductions`
  if (!is.null(reductions)) {
    if (all(reductions %in% names(dr))) {
      ReducNames <- reductions
    } else {
      stop("the reduction name provided not in the object")
    }
  } else {
    if (length(dr) >= 1) {
      ReducNames <- names(dr)
      message("Using all embeddings contained in the object: ", ReducNames)
    } else {
      message("please check the reduction imformation and run again")
    }
  }

  coords_list <- list()
  for (embedding in ReducNames) {
    emb <- Seurat::Embeddings(object = object, embedding)
    if (ncol(emb) > 2) {
      emb <- emb[, 1:2]
    }
    colnames(emb) <- c(sprintf("X_%s", embedding), sprintf("Y_%s", embedding))
    coords_list[[embedding]] <- as.matrix(emb)
  }



  ######## build init anndata object(mutlit exp)
  ## fetech the exp data
  # for (i in seq_len(assays)) {
  for (i in seq(length(assays))) {

    Seurat::DefaultAssay(object) <- assays[i]
    feature_meta <- data.frame(feature = rownames(object))
    rownames(feature_meta) <- feature_meta[["feature"]]

    ## attension the exp data fetch method between seurat v3 ~ seurat v5
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



    ### add motif pwm data
    # if (export_pwm && length(GetMotifData(object = object, slot = "pwm", assay = assays[i])) > 0 && assay.types[i] == "motif") {
    if (export_pwm && assay.types[i] == "motif") {
      peak_index <- which(assay.types == "peak")
      pwm_data <- Signac::GetMotifData(object = object, slot = "pwm", assay = assays[peak_index])
      motif_name <- Signac::GetMotifData(object = object, slot = "motif.names", assay = assays[peak_index])
      new_pwm_name <- as.vector(paste(names(pwm_data), motif_name, sep = "_"))
      names(pwm_data) <- new_pwm_name
      pwm_json <- jsonlite::toJSON(pwm_data, auto_unbox = TRUE)
      exp_adata$uns[["motif_logo"]] <- pwm_json
    }



    ##### add peak gene links
    if (export_links && assay.types[i] == "peak") {
    # if (export_links && length(Links(object = object[[assays[i]]])) > 0) {
      co_data <- as.data.frame(Signac::Links(object = object[[assays[i]]]))
      exp_adata$uns[["coacc_link"]] <- co_data
    }



    ######## add marker data
    # if (!is.null(markersDF)) {
    #   for (i in seq(length(markersDF))) {
    #     df_names <- names(markersDF)
    #     exp_adata$uns[["markers"]] <- markersDF[[i]]
    #   }
    # } else {
    #   message("no marker provided")
    # }

    if (!is.null(markersDF) && !is.null(markersDF[[assays[i]]])) {
      exp_adata$uns[["markers"]] <- markersDF[[assays[i]]]
    } else {
      message("no marker provided")
    }



    #### add coords informations
    for (eb in ReducNames) {
      exp_adata$obsm[[eb]] <- coords_list[[eb]]
    }

    ### judge the dir
    if (!dir.exists(outpath)) {
      dir.create(outpath)
    }


    #### write h5ad
    file_path <- file.path(outpath, sprintf("%s.h5ad", assays[i]))
    exp_adata$write_h5ad(file_path)

  }
}
