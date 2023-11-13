# require(anndata)
# require(Giotto)
# require(magick)
# require(imager)
# require(png)
# require(jsonlite)
# require(dplyr)

#' @title list_expression_names
#' @name list_expression_names
#' @param gobject gobject
#' @param spat_unit spatial unit ‘cell’
#' @param feat_type feature type ‘rna’
#' @description lists the available matrices names for a given spatial unit and feature type
#' @return vector with names of available matrices
list_expression_names = function(gobject,
                                 spat_unit = NULL,
                                 feat_type = NULL) {

  if(is.null(spat_unit)) stop('spat_unit must be given\n')
  if(is.null(feat_type)) stop('feat_type must be given\n')

  expression_names = names(gobject@expression[[spat_unit]][[feat_type]])

  return(expression_names)
}



#' @title set_default_feat_type
#' @name set_default_feat_type
#' @param gobject gobject
#' @param feat_type feature type ‘rna’
#' @param spat_unit spatial unit ‘cell’
set_default_feat_type = function(gobject,
                                 feat_type = NULL,
                                 spat_unit) {

  # set spatial unit
  if(is.null(feat_type)) {
    feat_type = getOption('giotto.feat_type')
    if(is.null(feat_type)) {
      if(!is.null(gobject@expression) & length(gobject@expression) > 0) {
        feat_type = names(gobject@expression[[spat_unit]])[[1]]
        if(is.null(feat_type)) stop('valid spat_unit input needed \n')
      } else if(!is.null(gobject@feat_info)){
        feat_type = names(gobject@feat_info)[[1]]
      } else {
        warning('No default for feat_type could be set \n')
      }
    }

  }

  return(feat_type)
}

#' @title set_default_spat_unit
#' @name set_default_spat_unit
#' @param gobject gobject
#' @param spat_unit spatial unit ‘cell’
#' @keywords internal
set_default_spat_unit = function(gobject,
                                 spat_unit = NULL) {


  # set spatial unit
  if(is.null(spat_unit)) {
    spat_unit = getOption('giotto.spat_unit')
    if(is.null(spat_unit)) {
      if(!is.null(gobject@expression) & length(gobject@expression) > 0) {
        spat_unit = names(gobject@expression)[[1]]
      } else if(!is.null(gobject@spatial_info)){
        spat_unit = names(gobject@spatial_info)[[1]]
      } else {
        warning('No default for spat_unit could be set \n')
      }
    }

  }
  return(spat_unit)
}


## expression values slot ####
#' @title  Get expression values
#' @name  get_expression_values
#' @description Function to get expression values from giotto object
#' @param gobject image object
#' @param spat_unit spatial unit ‘cell’
#' @param feat_type feature type ‘rna’
#' @param values expression values to extract (e.g. "raw", "normalized", "scaled")
#' @param output what object type to retrieve the expression as. Currently either
#' 'matrix' for the matrix object contained in the exprObj or 'exprObj' (default) for
#' the exprObj itself are allowed.
#' @importFrom methods new slot
#' @return expression matrix
#' @family expression accessor functions
#' @family functions to get data from giotto object
#' @export
get_expression_values = function(gobject,
                                 values = NULL,
                                 spat_unit = NULL,
                                 feat_type = NULL,
                                 output = c('exprObj', 'matrix')) {


  output = match.arg(output, choices = c('exprObj', 'matrix'))

  # 1. Set feat_type and spat_unit

  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # 2. Find object

  potential_values = list_expression_names(gobject = gobject,
                                           spat_unit = spat_unit,
                                           feat_type = feat_type)

  if(is.null(values)) values = potential_values[[1]]

  ## special cases for giotto standard pipeline
  if(values == 'scaled' & is.null(gobject@expression[[spat_unit]][[feat_type]][[values]])) {
    stop('run first scaling (& normalization) step')
  } else if(values == 'normalized' & is.null(gobject@expression[[spat_unit]][[feat_type]][[values]])) {
    stop('run first normalization step')
  } else if(values == 'custom' & is.null(gobject@expression[[spat_unit]][[feat_type]][[values]])) {
    stop('first add custom expression matrix')
  }

  if(!values %in% potential_values) stop("The spatial unit ", spat_unit ," for expression matrix ", feat_type, " and with name ","'", values, "'"," can not be found \n")

  # 3. Get object in desired format

  expr_values = gobject@expression[[spat_unit]][[feat_type]][[values]]

  if(output == 'exprObj') {
    if(!inherits(expr_values, 'exprObj')) {
      expr_values = methods::new('exprObj',
                        name = values,
                        exprMat = expr_values,
                        spat_unit = spat_unit,
                        feat_type = feat_type,
                        provenance = spat_unit, # assumed
                        misc = NULL)
    }

    if(!inherits(expr_values, 'exprObj')) stop('Cannot convert to exprObj')

    # return exprObj
    return(expr_values)

  } else if(output == 'matrix') {

    if(inherits(expr_values, 'exprObj')) expr_values = methods::slot(expr_values, 'exprMat')

    # return 'matrix'
    return(expr_values)
  }

}




#' @title getNoshiftSpatialLocs
#' @name Function getNoshiftSpatialLocs
#' @description Coveret shiftData To NoShiftData
#' @param object giotto object
#' @param cell_meta_anno cell meta table
#' @param spatial_locs spatial location
#' @importFrom dplyr mutate arrange
getNoshiftSpatialLocs <- function(object,cell_meta_anno = NULL,spatial_locs = NULL) {
  add_to_x_dict <- list()
  slice_names <- names(object@images)
  for (i in seq_along(slice_names)) {
    slice <- object@images[[i]]

    xmax_sloc <- slice@minmax[[1]]
    xmin_sloc <- slice@minmax[[2]]

    xmax_adj <- slice@boundaries[[1]]
    xmin_adj <- slice@boundaries[[2]]

    add_to_x <- (i-1) * (xmax_sloc - xmin_sloc + xmax_adj + xmin_adj) + (i-1) * 1000

    add_to_x_dict[[slice@name]] <- add_to_x
  }

  slice_names <- sub("^(.*?)\\-.*", "\\1", spatial_locs$cell_ID)

  slice_names <- paste0(slice_names, "-image")

  add_to_x_values <- unlist(add_to_x_dict[slice_names])


  processed_spatial_locs <- spatial_locs %>%
    dplyr::mutate(sdimx = sdimx - add_to_x_values,
           sdimy = abs(sdimy)) %>%
    dplyr::arrange(match(cell_ID, cell_meta_anno$list_ID))

  return(processed_spatial_locs)
}


#' @title giottoToAnnData
#' @name Function giottoToAnnData
#' @description Coveret Giotto To Anndata
#' @param object a giotto object
#' @param outpath the path to write h5ad
#' @param markerDF marker table to exprt
#' @importFrom anndata AnnData
#' @importFrom anndata write_h5ad
#' @importFrom Giotto get_spatial_locations
#' @importFrom Giotto get_dimReduction fDataDT pDataDT
#' @importFrom png readPNG
#' @importFrom jsonlite fromJSON
#' @importFrom utils read.table
#' @export
giottoToAnnData <- function(object = NULL,
                            outpath = NULL,
                            markerDF = NULL){

  # object <- visium_brain
  # Check gobject
  invalid_obj = !("giotto" %in% class(object))
  if (is.null(object) || invalid_obj) {
    stop(wrap_msg("Please provide a valid Giotto Object for conversion."))
  }

  # Check outpath directory, make it if it doesn't exist
  if (is.null(outpath)) outpath = paste0(getwd(),"/GiottoToAnndata.h5ad")

  raw_x = get_expression_values(gobject = object,
                                values = 'raw',
                                spat_unit = 'cell',
                                feat_type = 'rna',
                                output = "matrix")

  # cell_meta_anno <- data.frame(visium_brain@cell_metadata$cell$rna@metaDT)
  cell_meta_anno <- Giotto::pDataDT(object, spat_unit = NULL, feat_type = NULL)

  cell_meta_anno[['library_id']] <- cell_meta_anno$list_ID
  rownames(cell_meta_anno) <- cell_meta_anno$cell_ID

  gene_metadata <- Giotto::fDataDT(object, spat_unit = NULL, feat_type = NULL)
  rownames(gene_metadata) <- gene_metadata$feat_ID

  if (!is.null(raw_x) ) {
    adata <- anndata::AnnData(
      X=t(raw_x),
      obs = cell_meta_anno,
      var = gene_metadata)
  } else {
    message("raw_x information fetech failed, please make sure it is not null! ")
  }


  ## visium_brain@spatial_locs$cell$raw@coordinates
  # spatial_coordinates <- visium_brain@spatial_locs$cell$raw@coordinates
  spatial_locs <- Giotto::get_spatial_locations(gobject = object,
                                        output = "data.table",
                                        spat_unit = 'cell')


  if(!is.null(spatial_locs)) {
    # img_spatial <- sprintf("%s_spatial",visium_brain@images$image@name)
    # dim(obsm$spatial); obsm$spatial  [1]    2 2698
    # my_assay_int <- apply(my_assay, c(1, 2), as.integer)
    if (length(names(object@images)) > 1){
      spatial_locs <- getNoshiftSpatialLocs(object, cell_meta_anno, spatial_locs)
    }

    spatial_locs <- as.matrix(abs(spatial_locs[, c("sdimx", "sdimy")]))
    adata$obsm[['spatial']] <- spatial_locs
  } else {
    message("spatial_locs information fetech failed, please make sure it is not null! ")
  }


  if (!is.null(object@dimension_reduction)) {
    dr <- object@dimension_reduction$cells$cell$rna
    if (length(dr) >= 1){
      ReducNames <- list(names(dr))
      # ReducNames[[1]]
      for (embedding in ReducNames[[1]]) {
        print(embedding)
        adata$obsm[[embedding]] <- Giotto::get_dimReduction(
          object,
          spat_unit = NULL,
          feat_type = NULL,
          reduction = c("cells", "feats"),
          reduction_method = embedding,
          name = embedding,
          output = "data.table",
          set_defaults = TRUE
        )
        message("Dimension reduction information ", embedding, " is successfully stored in anndata object.")
      }
    }
  } else {
    message("dimension_reduction is null!")
  }


  if (!is.null(object@images)) {
    dr <- object@images
    if (length(dr) >= 1){
      ImagesNames <- list(names(dr))
      # ReducNames[[1]]
      for (image_name in ImagesNames[[1]]) {
        image_path <- object@images[[image_name]]@file_path
        image_path <- gsub("[\\/]+", "/", image_path)
        image_path <- normalizePath(image_path)
        print(image_path)
        img <- png::readPNG(image_path)
        img_arrary <- as.array(img)
        adata_image_name <- gsub("-image", "", image_name)
        adata$uns[["spatial"]][[adata_image_name]][["images"]][["lowres"]] <- img_arrary*255.0

        scale_factor_file_path <-  gsub("tissue_lowres_image.png", "scalefactors_json.json", image_path)
        file_content <- jsonlite::fromJSON(txt = scale_factor_file_path)
        spot_diameter_fullres <- file_content$spot_diameter_fullres
        tissue_hires_scalef <- file_content$tissue_hires_scalef
        tissue_lowres_scalef <- file_content$tissue_lowres_scalef
        fiducial_diameter_fullres <- file_content$fiducial_diameter_fullres
        scalefactors <- list()
        scalefactors$spot_diameter_fullres <- spot_diameter_fullres
        scalefactors$tissue_hires_scalef <- tissue_hires_scalef
        scalefactors$tissue_lowres_scalef <- tissue_lowres_scalef
        scalefactors$fiducial_diameter_fullres <- fiducial_diameter_fullres
        adata$uns[["spatial"]][[adata_image_name]][["scalefactors"]] <- scalefactors
      }
    }
  } else {
    message("images is null!")
  }

  if (is.character(markerDF)) {
    if (file.exists(markerDF)) {
      markerDF <- gsub("[\\/]+", "/", markerDF)
      absolutePath_markerDF <- normalizePath(markerDF)
      marker_data <- utils::read.table(absolutePath_markerDF, header = TRUE)
      adata$uns[["marker"]] <- marker_data
    } else {
      stop("invalid path")
    }
  } else if (exists("markerDF")) {

    marker_data <- as.data.frame(markerDF)
    adata$uns[["marker"]] <- marker_data
  } else {
    stop("invalid var")
  }

  ### write the adata into h5ad
  # adata$write_h5ad(file.path(outdir,sprintf("test.h5ad")))
  # outpath <- './test.h5ad'
  adata$write_h5ad(outpath)

  adata

}







