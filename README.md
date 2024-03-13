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
When Seurat contains multiple assays, users can provide multi assay names they wish to export. This will automatically generates the h5ad for each assay type, like RNA.h5ad,integration.h5ad etc. We also offer the flexibility for users to output marker dataframes as a list, with the corresponding information stored in **adata.uns['markers']**.
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
When converting spatial transcriptomic data, we automatically determine the spatial data type based on the object's structure. We return the corresponding **h5ad** file based on the identified type. The spatial coordinates are stored within **adata.obsm['spatial']**. If the analyzed object includes spatially organized slice information, we store it in **adata.uns['spatial']**. Multiple slices are differentiated using different **"library_id"**, ensuring clear distinction between them.
```
SeuratToAnndata(object=scRNA,
                outpath="/test_adata",
                assays=c("RNA", "Spatial"),
                groups = c("seurat_clusters", "region"),
                reductions = "tsne",
                markersDF = list("RNA"=rna_marker_df, "Spatial"=st_marker_df)) 
```

### Signac to Anndata
We recommend exporting **peak matrix**, **genescore matrix**, and **motif matrix** from single-cell ATAC data analysis as h5ad files. If desired, these files can be further converted into [**Mudata**](https://mudata.readthedocs.io/en/latest/), which offers a comprehensive and annotated multimodal dataset structure. Additionally, users can use **export_pwm** to output results from co-accessibility analysis and motif enrichment analysis, respectively. Users can access motif-related information via **adata.uns.["motifs"]**.
```
SignacToAnndata(object=scATAC,
                outpath = "/test_adata",
                assays=c("RNA", "Peaks", "Motif"),
                assay.types=("gene","peak","motif"),
                markersDF = list("RNA"=gene_marker_df, "Peaks"=marker_peak_df, "Motif"=marker_motif_df),
                groups = c("clusters","sample", "age"),
                reductions = c("umap","tsne"),
                export_pwm = TRUE)   # to export motif pwm 
```

### ArchR to AnnData
Similar to the Signac object conversion, users can convert ArchR objects and export relevant information based on their needs.
```
ArchrToAnndata(object=project5,
               outpath="/test_adata",
               assays=c("RNA", "Peaks", "Motif"),
               assay.types=("gene","peak","motif"),
               markersDF=list("RNA"=gene_marker_df, "Peaks"=marker_peak_df, "Motif"=marker_motif_df),
               groups = NULL, ## export all cell meta
               reductions = NULL, ## export all reduction coords
               export_pwm =TRUE) 
```

### Giotto to AnnData
Similar to Seurat's conversion of spatial transcriptomic analysis results, we seamlessly incorporate tissue slices and spatial coordinates from Giotto object into the output h5ad object by default.
```
giottoToAnnData(object = giotto,
                outpath = "/test_adata",
                markerDF = NULL)
```












