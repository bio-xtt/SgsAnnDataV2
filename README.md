# SgsAnnDataV2
A package used to convert Seurat, Giotto, Signac, ArchR analysis object into AnnData format 

#### Example1
Convert scRNA Seurat object into AnnData, this automatically generates RNA.h5ad and SCT.h5ad respectively in the output folder
(```)
SeuratToAnndata(object=scRNA,
                outpath="/test_adata",
                assays=c("RNA", "SCT"),
                groups = NULL,
                reductions = NULL,
                markersDF = list("RNA"=rna_marker_df, "SCT"=sct_marker_df)) 
(```)

#### Example2
Convert Spatial Seurat object into AnnData, this automatically generates RNA.h5ad and Spatial.h5ad respectively in the output folder
(```)
SeuratToAnndata(object=scRNA,
                outpath="/test_adata",
                assays=c("RNA", "Spatial"),
                groups = NULL,
                reductions = NULL,
                markersDF = list("RNA"=rna_marker_df, "Spatial"=st_marker_df)) 
(```)

#### Example3
Convert Signac object into AnnData, this automatically generates RNA.h5ad and Spatial.h5ad respectively in the output folder
(```)
SeuratToAnndata(object=scATAC,
                outpath = "/test_adata",
                assays=c("RNA", "Peaks", "Motif"),
                assay.types=("gene","peak","motif"),
                markersDF = list("RNA"=gene_marker_df, "Peaks"=marker_peak_df, "Motif"=marker_motif_df),
                groups = NULL,
                reductions = NULL,
                export_links = TRUE,  # to export coaccess link 
                export_pwm = TRUE)   # to export motif pwm 
(```)

#### Example4
Convert ArchR object into AnnData, this automatically generates RNA.h5ad, Peaks.h5ad and Motif.h5ad respectively in the output folder
(```)
ArchrToAnndata(object=project5,
               outpath="/test_adata",
               assays=c("RNA", "Peaks", "Motif"),
               assay.types=("gene","peak","motif"),
               markersDF=list("RNA"=gene_marker_df, "Peaks"=marker_peak_df, "Motif"=marker_motif_df),
               groups = NULL,
               reductions = NULL,
               export_links = FALSE,
               export_pwm =TRUE) 
(```)

#### Example5
Convert Giotto object into AnnData, this automatically generates giotto.h5ad in the output folder
(```)
giottoToAnnData <- function(object = giotto,
                            outpath = "/test_adata",
                            markerDF = NULL)
(```)












