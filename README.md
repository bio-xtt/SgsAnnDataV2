# SgsAnnDataV2
SgsAnnDataV2 is an R package that facilitates the seamless conversion of single-cell analysis object from popular tools such as Seurat, Giotto, Signac, and ArchR into the AnnData. This format can be directly visualized in SGS which is an interactive browser for single-cell and spatial multimodal datasets.

> ## Warning
> Please note that the usage of SGS requires the prior installation of R packages such as Seurat, Giotto, Signac, and anndata. It is important to ensure that these packages are installed beforehand for smooth execution of SGS workflows.

## Installation
Currently, SgsAnnDataV2 can only be installed through GitHub. However, we are actively working on providing a BiocManager installation method in the near future. 
```
devtools::install_github("bio-xtt/SgsAnnDataV2")
```

## Usage
### Seurat to AnnData
### Single-Cell Transcriptome Object Conversion
When Seurat contains multiple assays, users can provide multi assay names they wish to export. This will automatically generates the h5ad for each assay type, like RNA.h5ad,integration.h5ad etc. 
```
library(Seurat)
library(anndata)
SeuratToAnndata(object=scRNA,
                outpath="/test_adata",
                assays=c("RNA", "SCT"),
                groups = NULL,
                reductions = c("tsne","umap"),
                markersDF = list("RNA"=rna_marker_df, "SCT"=sct_marker_df)) 
```
### Spatial Transcriptome Object Conversion
Convert Spatial Seurat object into AnnData, this automatically generates RNA.h5ad and Spatial.h5ad respectively in the output folder

When converting spatial transcriptomic data, we automatically determine the spatial data type based on the object's structure. We return the corresponding **h5ad** file based on the identified type. The spatial coordinates are stored within **adata.obsm['spatial']**. If the analyzed object includes spatially organized slice information, we store it in **adata.uns['spatial']**. Multiple slices are differentiated using different **"library_id"**, ensuring clear distinction between them.
```
SeuratToAnndata(object=scRNA,
                outpath="/test_adata",
                assays=c("RNA", "Spatial"),
                groups = c("seurat_clusters", "region"),
                reductions = "tsne",
                markersDF = list("RNA"=rna_marker_df, "Spatial"=st_marker_df)) 
```

#### Example3
Convert Signac object into AnnData, this automatically generates RNA.h5ad and Spatial.h5ad respectively in the output folder
```
SeuratToAnndata(object=scATAC,
                outpath = "/test_adata",
                assays=c("RNA", "Peaks", "Motif"),
                assay.types=("gene","peak","motif"),
                markersDF = list("RNA"=gene_marker_df, "Peaks"=marker_peak_df, "Motif"=marker_motif_df),
                groups = NULL,
                reductions = NULL,
                export_links = TRUE,  # to export coaccess link 
                export_pwm = TRUE)   # to export motif pwm 
```
### Siganac to Anndata




### ArchR to adata
Convert ArchR object into AnnData, this automatically generates RNA.h5ad, Peaks.h5ad and Motif.h5ad respectively in the output folder
```
ArchrToAnndata(object=project5,
               outpath="/test_adata",
               assays=c("RNA", "Peaks", "Motif"),
               assay.types=("gene","peak","motif"),
               markersDF=list("RNA"=gene_marker_df, "Peaks"=marker_peak_df, "Motif"=marker_motif_df),
               groups = NULL,
               reductions = NULL,
               export_links = FALSE,
               export_pwm =TRUE) 
```

#### Example5
Convert Giotto object into AnnData, this automatically generates giotto.h5ad in the output folder
```
giottoToAnnData <- function(object = giotto,
                            outpath = "/test_adata",
                            markerDF = NULL)
```












