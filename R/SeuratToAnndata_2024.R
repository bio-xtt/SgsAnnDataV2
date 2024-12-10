## This script used to cover seurat object into anndata object
## 2024.12.10

#' Function used to parser spatial informs of different ST-seq
#' @param object Seurat object
#' @param img.name name of the image
#' @param adata  AnnData object
#' @return a list
#' @importFrom Seurat GetTissueCoordinates ScaleFactors
#'

Gain_spatial <- function(object, img.name, adata) {
  ## "SlideSeq","VisiumV1","STARmap"(seqbased);"FOV"(image-based)
  if (class(object[[img.name]]) == "VisiumV1") {
    spatial_coords <- object@`images`[[img.name]]@`coordinates`[, c("imagecol", "imagerow")]
    scale <- as.list(Seurat::ScaleFactors(object[[img.name]]))
    scale_list <- list(
      "tissue_hires_scalef" = scale$`hires`,
      "tissue_lowres_scalef" = scale$`lowres`,
      "fiducial_diameter_fullres" = scale$`fiducial`,
      "spot_diameter_fullres" = scale$`spot`
    )
    img_arrary <- object@`images`[[img.name]]@`image`
    adata$uns[["spatial"]][[img.name]][["images"]][["lowres"]] <- img_arrary * 255
    adata$uns[["spatial"]][[img.name]][["scalefactors"]] <- scale_list  
    
  } else if (class(object[[img.name]]) == "VisiumV2"){
    # spatial_coords <- GetTissueCoordinates(object = object@images[[img.name]])[,c("x","y")]
    spatial_coords <- GetTissueCoordinates(object = object@images[[img.name]])[,c("y","x")]
    scale <-  as.list(Seurat::ScaleFactors(object[[img.name]]))
    scale_list <- list(
      "tissue_hires_scalef" = scale$`hires`,
      "tissue_lowres_scalef" = scale$`lowres`,
      "fiducial_diameter_fullres" = scale$`fiducial`,
      "spot_diameter_fullres" = scale$`spot`)
    img_arrary <- object@`images`[[img.name]]@`image`
    adata$uns[["spatial"]][[img.name]][["images"]][["lowres"]] <- img_arrary * 255
    adata$uns[["spatial"]][[img.name]][["scalefactors"]] <- scale_list
      
  } else if (class(object[[img.name]]) == "SlideSeq") {
    spatial_coords <- Seurat::GetTissueCoordinates(object = object, image = img.name)[, c("y", "x")]
      
  } else if (class(object[[img.name]]) == "STARmap") {
    spatial_coords <- Seurat::GetTissueCoordinates(object = object, image = img.name)[, c("y", "x")]
      
  } else if (class(object@`images`[[img.name]]) == "FOV") {
    spatial_coords <- Seurat::GetTissueCoordinates(object = object[[img.name]], which = "centroids")
    rownames(spatial_coords) <- spatial_coords$`cell`
    spatial_coords <- spatial_coords[, c("y", "x")]
    if ("segmentation" %in% names(object[[img.name]])) {
      Segmentation_coords <- Seurat::GetTissueCoordinates(object = object[[img.name]], which = "segmentation")
      Segmentation_coords <- as.matrix(Segmentation_coords[, c("x", "y")])
      ### add the segamnet data into uns
      adata$uns[["spatial"]][[img.name]][["images"]][["Segmentation_coords"]] <- Segmentation_coords
    }
      
  } else {
    stop("Spatial information fetch failed")
  }
  library_id <- data.frame(library_id = rep(img.name, nrow(spatial_coords)), cell_id = rownames(spatial_coords))
  rownames(library_id) <- library_id$`cell_id`
  result <- list(adata = adata, sp_coords = spatial_coords, library_id = library_id)
  return(result)
}




# Function used to ccovert seurat object into anndata(multi assay)
#' @title SeuratToAnndata
#' @name Function SeuratToAnndata
#' @param object Seurat object
#' @param outpath output dir
#' @param assays vector of assay names to export
#' @param groups vector of groups to export, if null export all
#' @param reductions vector of reduction names to export, if null export all
#' @param markersDF a named list of marker df vars,need the colnames in c(gene, cluster, ...) order
#' @importFrom Seurat GetAssayData GetTissueCoordinates ScaleFactors  DefaultAssay GetAssayData Images Embeddings
#' @importFrom anndata AnnData write_h5ad
#' @importFrom utils packageVersion
#' @export
#'
SeuratToAnndata <- function(object,
                            outpath,
                            assays,
                            groups = NULL,
                            reductions = NULL,
                            markersDF = NULL) {
  # This function mainly designed to convert Seurat V3 to V5 objects
  message("Seurat Version installed: ", packageVersion("Seurat"))
  message("Object was created with Seurat version ", object@`version`)


  # Get cell meta data
  cell_meta <- object[[]]
  cell_meta$`cell_idex` <- as.vector(seq(rownames(cell_meta)))
  # remove NA column
  cell_meta =  cell_meta[, apply(cell_meta, 2, function(y) any(!is.na(y)))]
  if (!is.null(groups) && groups %in% colnames(cell_meta)) {
    groups <- append(groups, "cell_idex")
    cell_meta <- cell_meta[, groups]
  }
  # Add cell id column
  cell_meta$`cell_id` <- rownames(cell_meta)


                                 
  # Export assays
  for (i in seq_along(assays)) {
    Seurat::DefaultAssay(object) <- assays[i]
    feature_meta <- data.frame(feature = rownames(object))
    rownames(feature_meta) <- feature_meta[["feature"]]
    ## Attension the exp data fetch method between seurat v3 ~ seurat v5
    if (class(object@assays[[i]]) == "Assay") {
      if (!is.null(object@assays[[i]]@`data`)) {
        exp <- Seurat::GetAssayData(object = object, assay = assays[i], slot = "data")
        adata <- anndata::AnnData(
          X = Matrix::t(exp),
          obs = cell_meta,
          var = feature_meta
        )
      } else {
        exp <- Seurat::GetAssayData(object = object, assay = assays[i], slot = "counts")
        adata <- anndata::AnnData(
          X = Matrix::t(exp),
          obs = cell_meta,
          var = feature_meta
        )
      }
    } else if (class(object@assays[[i]]) == "Assay5") {
      options(Seurat.object.assay.version = "v5")

      if (!is.null(object@assays[[i]]@`layers`[["data"]])) {
        exp <- Seurat::GetAssayData(object = object, assay = assays[i], layer = "data")
        adata <- anndata::AnnData(
          X = Matrix::t(exp),
          obs = cell_meta,
          var = feature_meta
        )
      } else {
        exp <- Seurat::GetAssayData(object = object, assay = assays[i], layer = "counts")
        adata <- anndata::AnnData(
          X = Matrix::t(exp),
          obs = cell_meta,
          var = feature_meta
        )
      }
    }


    # Add coords informations
    dr <- object@`reductions`
    !is.null(reductions){
      ReducNames <- intersect(reductions, names(dr))
      if (length(ReducNames) == 0) {
      stop("the reduction name provided not in the object")
    }
    } else {
      ReducNames <- names(dr)
      message("Using all embeddings contained in the object: ", paste(ReducNames, collapse = ", "))
    }
    message("Using all embeddings contained in the object: ", paste(ReducNames, collapse = ", "))
    for (embedding in ReducNames) {
    if (assays[i] == object@reductions[[embedding]]@`assay.used`) {
        emb <- Seurat::Embeddings(object = object, embedding)
        if (ncol(emb) > 2) {
          emb <- emb[, 1:2]
        }
        colnames(emb) <- c(sprintf("X_%s", embedding), sprintf("Y_%s", embedding))
        # coords_list[[embedding]] <- as.matrix(emb)
        adata$obsm[[embedding]] <- as.matrix(emb)
      } else if (nrow(Seurat::Embeddings(object = object, embedding)) == ncol(object@`assays`[[assays[i]]])) {  
        emb <- Seurat::Embeddings(object = object, embedding)
        if (ncol(emb) > 2) {
          emb <- emb[, 1:2]
        }
        colnames(emb) <- c(sprintf("X_%s", embedding), sprintf("Y_%s", embedding))
        # coords_list[[embedding]] <- as.matrix(emb)
        adata$obsm[[embedding]] <- as.matrix(emb)
            
      } else {
        message("The dimensionality reduction coordinates are not computed from this matrix")
        break
      }
    }
        

    

    # Add spatial informations if exist
    if (!is.null(Seurat::Images(object = object, assay = assays[i])) && length(Images(object = object, assay = assays[i])) > 0) {
      images <- Seurat::Images(object = object, assay = assays[i])
    } else {
      images <- Seurat::Images(object = object)
    }
    imageNames <- images

    if (length(imageNames) >= 1 && class(object[[imageNames[1]]]) == "VisiumV1") {
      # Merged all spatial coords
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
    } else if (length(imageNames) >= 1 && class(object[[imageNames[1]]]) == "VisiumV2") {
        # "VisiumV2"
        all_coords <- data.frame()
        library_ids <- data.frame()
        for (img in imageNames) {
          result <- Gain_spatial(object = object, img.name = img, adata = adata)
          spatial_coords <- result[["sp_coords"]]
          all_coords <- rbind(spatial_coords, all_coords)
          library_id <- result[["library_id"]]
          library_ids <- rbind(library_id, library_ids)
        }
        adata <- result[["adata"]]
        ## filter the object
        all_coords <- all_coords[rownames(adata$`obs`), ]
        library_ids <- library_ids[rownames(adata$`obs`), ]
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


    # Add marker data
    if (!is.null(markersDF) && !is.null(markersDF[[assays[i]]])) {
      adata$uns[["markers"]] <- markersDF[[assays[i]]]
    } else {
      message("no marker provided")
    }


    # Create output dir
    if (!dir.exists(outpath)) {
      dir.create(outpath)
    }

    # Write h5ad
    file_path <- file.path(outpath, sprintf("%s.h5ad", assays[i]))
    #adata$write_h5ad(file_path)
    anndata::write_h5ad(adata, file_path)
  }
}
