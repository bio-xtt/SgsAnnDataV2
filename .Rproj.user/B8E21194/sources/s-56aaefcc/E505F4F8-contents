## This script used to cover seurat object into anndata object
## 2023.8.3
# require(Seurat)
# require(anndata)
# require(Matrix)


###### Function used to parser spatial informs of different ST-seq
#' @param object seurat object
#' @param img.name name of the image
#' @param adata adata object
#' @return a list
#' @importFrom Seurat GetTissueCoordinates ScaleFactors
#'

Gain_spatial <- function(object ,img.name, adata) {
  ## "SlideSeq","VisiumV1","STARmap"(seqbased);"FOV"(image-based)
  ## class(object@`images`[[img.name]]) == "VisiumV1"
  if (class(object[[img.name]]) == "VisiumV1") {
    ## this was changed by xtt
    #spatial_coords <- GetTissueCoordinates(object = object, image = img.name)[, c("imagecol", "imagerow")]

    spatial_coords <- object@`images`[[img.name]]@`coordinates`[ ,c("imagecol", "imagerow")]
    scale <- as.list(Seurat::ScaleFactors(object[[img.name]]))
    scale_list <- list(
      "tissue_hires_scalef" = scale$`hires`,
      "tissue_lowres_scalef" = scale$`lowres`,
      "fiducial_diameter_fullres" = scale$`fiducial`,
      "spot_diameter_fullres" = scale$`spot`)

    img_arrary <- object@`images`[[img.name]]@`image`
    adata$uns[["spatial"]][[img.name]][["images"]][["lowres"]] <- img_arrary * 255
    adata$uns[["spatial"]][[img.name]][["scalefactors"]] <- scale_list

  } else if (class(object[[img.name]]) == "SlideSeq"){

    # SlideSeq没有图片，故也不需要存储相应的缩放因子
    spatial_coords <- Seurat::GetTissueCoordinates(object = object, image = img.name)[, c("y", "x")]

  } else if (class(object[[img.name]]) == "STARmap") {
    # 目前暂不确定STARmap是否与slideseq存储数据一致，故分别存储
    spatial_coords <- Seurat::GetTissueCoordinates(object = object, image = img.name)[, c("y", "x")]

  } else if (class(object@`images`[[img.name]]) == "FOV") {
    # 目前这个主要用于存储基于成像的空间转录组数据
    spatial_coords <- Seurat::GetTissueCoordinates(object = object[[img.name]], which = "centroids")
    rownames(spatial_coords) <- spatial_coords$`cell`
    spatial_coords <- spatial_coords[, c("y", "x")]

    # 处理细胞分割信息
    if ("segmentation" %in% names(object[[img.name]])) {
      Segmentation_coords <- Seurat::GetTissueCoordinates(object = object[[img.name]], which = "segmentation")
      Segmentation_coords <- as.matrix(Segmentation_coords[, c("x", "y")])
      ### add the segamnet data into uns
      adata$uns[["spatial"]][[img.name]][["images"]][["Segmentation_coords"]] <- Segmentation_coords
    }

  } else {
    message("spatial information fetech failed")
    stop()
  }

  #### 将lirary id添加到细胞meta表
  library_id <- data.frame(library_id = rep(img.name, nrow(spatial_coords)), cell_id = rownames(spatial_coords))
  rownames(library_id) <- library_id$`cell_id`

  result <- list(adata = adata, sp_coords = spatial_coords, library_id = library_id)
  return(result)

}




# Function used to ccovert seurat object into anndata(multi assay)
#' @title SeuratToAnndata
#' @name Function SeuratToAnndata
#' @param object seurat object
#' @param outpath outputpath of the anndata object
#' @param assays a vector of assay names to export in adata
#' @param groups vector of groups to export,if null export all
#' @param reductions vector of reduction names to export,if null export all
#' @param markersDF a dataframe of the marker files,need the colnames in c(gene, cluster, ...) order
#' @importFrom Seurat GetAssayData GetTissueCoordinates ScaleFactors  DefaultAssay GetAssayData Images Embeddings
#' @importFrom anndata AnnData write_h5ad
#' @importFrom utils packageVersion
#' @export
# @useDynLib sgsAnnData2
#'

SeuratToAnndata <- function(object,
                            outpath,
                            assays,
                            groups = NULL,
                            reductions = NULL,
                            markersDF = NULL) {

  #### this ccovert function mainly for seurat V3-seurat V5
  message("Seurat Version installed: ", packageVersion("Seurat"))
  message("Object was created with Seurat version ", object@`version`)


  #### get the cell meta data
  cell_meta <- object[[]]
  # # cell_meta$cell_id <- as.vector(seq(rownames(cell_meta)) )
  cell_meta$`cell_idex` <- as.vector(seq(rownames(cell_meta)))
  if (!is.null(groups) && groups %in% colnames(object[[]])) {
    cell_meta <- object[[]]
    groups <- append(groups, "cell_idex")
    # cell_meta <- cell_meta[, which(colnames(cell_meta) %in% groups)]
  } else {
    cell_meta <- object[[]]
  }
  ### add cell id column
  cell_meta$cell_id <- rownames(cell_meta)


  ### get cell embedding
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



  ### fetch the exp data
  for (i in seq(length(assays))) {


    ### fetch the feature meta
    Seurat::DefaultAssay(object) <- assays[i]
    feature_meta <- data.frame(feature = rownames(object))
    rownames(feature_meta) <- feature_meta[["feature"]]


    ## attension the exp data fetch method between seurat v3 ~ seurat v5
    if (class(object@assays[[i]]) == "Assay") {
      if (!is.null(object@assays[[i]]@`data`)) {
        adata <- anndata::AnnData(
          X = Matrix::t(Seurat::GetAssayData(object = object, assay = assays[i], slot = "data")),
          obs = cell_meta,
          var = feature_meta
        )
      } else {
        adata <- anndata::AnnData(
          X = Matrix::t(Seurat::GetAssayData(object = object, assay = assays[i], slot = "counts")),
          obs = cell_meta,
          var = feature_meta
        )
      }
    } else if (class(object@assays[[i]]) == "Assay5") {
      options(Seurat.object.assay.version = "v5")

      if (!is.null(object@assays[[i]]@`layers`[["data"]])) {
        adata <- anndata::AnnData(
          X = Matrix::t(Seurat::GetAssayData(object = object, assay = assays[i], layer = "data")),
          obs = cell_meta,
          var = feature_meta
        )
      } else {
        adata <- anndata::AnnData(
          X = Matrix::t(Seurat::GetAssayData(object = object, assay = assays[i], layer = "counts")),
          obs = cell_meta,
          var = feature_meta
        )
      }
    }


    ##### add coords informations
    for (eb in ReducNames) {
      adata$obsm[[eb]] <- coords_list[[eb]]
    }


    ##### add spatial informations if exist
    if (length(Images(object = object, assay = assays[i])) > 0) {
      images <- Seurat::Images(object = object, assay = assays[i])
    } else {
      images <- Seurat::Images(object = object)
    }
    imageNames <- images

    if (length(imageNames) >= 1 && class(object[[imageNames[1]]]) == "VisiumV1") {

      ## changed to merged all spatial coords together
      all_coords <- data.frame()
      library_ids <- data.frame()
      for (img in imageNames) {
        result <- Gain_spatial(object = object, img.name = img, adata = adata)

        spatial_coords <- result[["sp_coords"]]
        all_coords <- rbind(spatial_coords, all_coords)

        library_id <- result[["library_id"]]
        library_ids <- rbind(library_id, library_ids)
      }

      all_coords <- all_coords[rownames(object[[]]), ]
      library_ids <- library_ids[rownames(object[[]]), ]
      adata <- result[["adata"]]
      adata$obsm[["spatial"]] <- as.matrix(all_coords)
      adata$obs[["library_id"]] <- library_ids$`library_id`
    } else if (length(imageNames) >= 1) {
      for (img_names in imageNames) {
        result <- Gain_spatial(object = object, img.name = img_names, adata = adata)
      }
      spatial_coords <- result[["sp_coords"]]
      library_id <- result[["library_id"]]
      adata <- result[["adata"]]
      adata$obsm[["spatial"]] <- as.matrix(spatial_coords)
      adata$obs[["library_id"]] <- library_id$`library_id`
    } else {
      message("Object has no image informations")
    }


    ######## add marker data
    if (!is.null(markersDF) && !is.null(markersDF[[assays[i]]])) {

      adata$uns[["markers"]] <- markersDF[[assays[i]]]
    } else {
      message("no marker provided")
    }


    ### judge the dir
    if (!dir.exists(outpath)) {
      dir.create(outpath)
    }


    #### write h5ad
    file_path <- file.path(outpath, sprintf("%s.h5ad", assays[i]))
    adata$write_h5ad(file_path)
  }
}






























































