# ScAT: Stereo-seq and scRNA-seq Analysis Toolkit

## Overview

```bash
> ${scat_dir}/bin/ScAT -h

usage: ScAT [-h] [-v]
            {QcControl,Clustering,Trajectory,Annotation,AnnotscCATCH,Enrichment,CellCellInteraction,TFregulon,SpatialCluster,DeconvSpot,SpatialPatternGene}
            ...

Stereo-seq and scRNA-seq Analysis Toolkit (ScAT), a user-friendly Python
package that provides multiple modules for scRNA-seq and Stereo-seq data
analysis.

positional arguments:
  {QcControl,Clustering,Trajectory,Annotation,AnnotscCATCH,Enrichment,CellCellInteraction,TFregulon,SpatialCluster,DeconvSpot,SpatialPatternGene}
    QcControl           Perform QC control by Seurat
    Clustering          Perform cell cluster by Seurat
    Trajectory          Perform trajectory analysis by monocle3
    Annotation          Perform cell annotation by SingleR
    AnnotscCATCH        Perform cell annotation by scCATCH
    Enrichment          Perform go enrichment by clusterProfiler
    CellCellInteraction
                        Perform cell-cell communication analysis by CellChat
    TFregulon           Perform single-cell regulatory network inference.
    SpatialCluster      Perform spatial cluster by BayeSpaces
    DeconvSpot          Perform deconvolution of cell types by SPOTlight
    SpatialPatternGene  Statistical analysis of spatial expression pattern by
                        SPARK

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         print version info.

Use ScAT {command} -h to get help on individual commands

```

## Installation & Environment
```bash
## Download ScAT source code
git clone https://github.com/liaoshangfeng/ScAT.git
## Create environment
conda create -n ScAT
wget https://github.com/liaoshangfeng/ScAT/blob/main/env/install_env.sh
bash install_env.sh
```

## Usage

### ScAT QcControl

```bash
> ${scat_dir}/bin/ScAT QcControl -h
usage: ScAT QcControl [-h] --Rscript RSCRIPT --input INPUT [--binsize BINSIZE]
                      [--geneColumn GENECOLUMN] [--fileType {gem,matrix}]
                      --out OUT [--prefix PREFIX] [--ptsize PTSIZE]
                      [--mtPattern MTPATTERN] [--rbPattern RBPATTERN]
                      [--hbPattern HBPATTERN] [--platPattern PLATPATTERN]
                      [--mFilter] [--minCount MINCOUNT] [--maxCount MAXCOUNT]
                      [--minFeature MINFEATURE] [--maxFeature MAXFEATURE]
                      [--minCell MINCELL] [--mtRatio MTRATIO]
                      [--rbRatio RBRATIO] [--hbRatio HBRATIO]
                      [--platRatio PLATRATIO]

optional arguments:
  -h, --help            show this help message and exit

input arguments:
  --Rscript RSCRIPT     Rscript path (default: None)
  --input INPUT         Input Stereo-seq *.gem; 10X single cell directory;
                        raw_matrix.txt.gz; raw_seuratObj.rds filename
  --binsize BINSIZE     Bin size to binning, your input should be in binSize1
                        if you set this (default 1)
  --geneColumn GENECOLUMN
                        Specify which column of genes.tsv or features.tsv to
                        use for gene names (default 2)
  --fileType {gem,matrix}
                        Set as matrix to support *.matrix.txt.gz, otherwise,
                        no need to change (default gem)

output arguments:
  --out OUT             Directory to save file
  --prefix PREFIX       sample ID, will be used as output prefix and seurat
                        object ident (default sample)
  --ptsize PTSIZE       Point size for spatial dim plot (default 1)

quality control arguments:
  --mtPattern MTPATTERN
                        Mitochondria gene pattern (default ^MT-)
  --rbPattern RBPATTERN
                        Ribosome gene pattern (default ^RP[SL])
  --hbPattern HBPATTERN
                        Hemoglobin gene pattern (default ^HB[^(P)])
  --platPattern PLATPATTERN
                        platelet gene pattern (default PECAM1|PF4)

filter arguements:
  --mFilter             Add "--mTop_markers" to allow filteration (default:
                        False)
  --minCount MINCOUNT   Minimum UMI number per cell/BIN [default 0]
  --maxCount MAXCOUNT   Maximum UMI number per cell/BIN [default 1000000]
  --minFeature MINFEATURE
                        Minimum feature number per cell/BIN [default 0]
  --maxFeature MAXFEATURE
                        Maximum feature number per cell/BIN [default 100000]
  --minCell MINCELL     Minimum cell or BIN number express the gene (for genes
                        filterationn) (default 3)
  --mtRatio MTRATIO     Maximum mitochondria gene UMI ratio per cell/BIN
                        (default 100)
  --rbRatio RBRATIO     Maximum ribosome gene UMI ratio per cell/BIN (default
                        100)
  --hbRatio HBRATIO     Maximum hemoglobin gene UMI ratio per cell/BIN
                        (default 100)
  --platRatio PLATRATIO
                        Maximum platelet gene UMI ratio per cell/BIN (default
                        100)


```

**Note** geneColumn = 1 when using C4 data

### ScAT Clustering

```bash
> ${scat_dir}/bin/ScAT Clustering -h
/data/work/previous/ScAT_github/core/bin/ScAT.py
usage: ScAT Clustering [-h] --Rscript RSCRIPT --input INPUT [--pc PC]
                       [--resolution RESOLUTION] [--ptsize PTSIZE]
                       [--topN {1,2,3,4,5,6,7,8,9}] --out OUT
                       [--prefix PREFIX] [--findSpatialVar]

optional arguments:
  -h, --help            show this help message and exit

input arguments:
  --Rscript RSCRIPT     Rscript path (default: None)
  --input INPUT         Input *_seuratObj.rds with SCT assay
  --pc PC               PCA dimensionality for clustering. (default: 30)
  --resolution RESOLUTION
                        Resolution for clustering, use a value above (below)
                        1.0 if you want obtain a larger (smaller) number of
                        communities. (default: 0.8)
  --ptsize PTSIZE       Point size for spatil dim plot
  --topN {1,2,3,4,5,6,7,8,9}
                        Select top N makers in each cluster for further
                        visualization (default: 3)

output arguments:
  --out OUT             Directory to save file
  --prefix PREFIX       sample ID, will be used as output file prefix
                        (default: sample)

boolean arguments:
  --findSpatialVar      Add "--findSpatialVar" to find spatial variable
                        features by Seurat (default: False)

```

### ScAT Annotation

```bash
> ${scat_dir}/bin/ScAT Annotation -h
usage: ScAT Annotation [-h] --Rscript RSCRIPT --input INPUT
                       [--ref_name {BlueprintEncodeData,DatabaseImmuneCellExpressionData,HumanPrimaryCellAtlasData,ImmGenData,MonacoImmuneData,MouseRNAseqData,NovershternHematopoieticData}]
                       [--cluster CLUSTER] --out OUT [--prefix PREFIX]

optional arguments:
  -h, --help            show this help message and exit

input arguments:
  --Rscript RSCRIPT     Rscript path (default: None)
  --input INPUT         Input normalized data (*txt.gz)
  --ref_name {BlueprintEncodeData,DatabaseImmuneCellExpressionData,HumanPrimaryCellAtlasData,ImmGenData,MonacoImmuneData,MouseRNAseqData,NovershternHematopoieticData}
                        Ref database name (default: HumanPrimaryCellAtlasData)
  --cluster CLUSTER     cluster file generated by Seurat workflow (default:
                        None)

output arguments:
  --out OUT             Directory to save file
  --prefix PREFIX       sample ID, will be used as output file prefix
                        (default: sample)
        
```

### ScAT AnnotscCATCH

```bash
> ${scat_dir}/bin/ScAT AnnotscCATCH -h
usage: ScAT AnnotscCATCH [-h] --Rscript RSCRIPT --input INPUT --species
                         SPECIES --tissue TISSUE [--cancer CANCER] --out OUT
                         [--prefix PREFIX]

optional arguments:
  -h, --help         show this help message and exit

input arguments:
  --Rscript RSCRIPT  Rscript path (default: None)
  --input INPUT      Rds file of gene expression
  --species SPECIES  Select species, support Human or Mouse (default: None)
  --tissue TISSUE    tiisue type: https://github.com/ZJUFanLab/scCATCH/wiki
                     (default: None)
  --cancer CANCER    cancer type:https://github.com/ZJUFanLab/scCATCH/wiki
                     (default: None)

output arguments:
  --out OUT          Directory to save file
  --prefix PREFIX    sample ID, will be used as output file prefix (default:
                     sample)

```

### ScAT Enrichment

```bash
> ${scat_dir}/bin/ScAT Enrichment -h
usage: ScAT Enrichment [-h] --Rscript RSCRIPT --input INPUT --species
                       {Human,Mouse} --cluster CLUSTER --out OUT

optional arguments:
  -h, --help            show this help message and exit

input arguments:
  --Rscript RSCRIPT     Rscript path (default: None)
  --input INPUT         Input markerGene filename which generated by Seurat
                        workflow.
  --species {Human,Mouse}
                        Species name
  --cluster CLUSTER     Cluser file name to compare go enrichment (one column)

output arguments:
  --out OUT             Directory to save file

```

### ScAT Trajectory
```bash
> ${scat_dir}/bin/ScAT Trajectory -h
usage: ScAT Trajectory [-h] --Rscript RSCRIPT --input INPUT [--gene GENE]
                       [--mDEG] [--mGCN] [--mTop_markers] --out OUT
                       [--prefix PREFIX]

optional arguments:
  -h, --help         show this help message and exit

input arguments:
  --Rscript RSCRIPT  Rscript path (default: None)
  --input INPUT      Input RDS file with seurat object, e.g.
                     *flt_norm_cluster_SeuratObject.
  --gene GENE        Gene name to define root of trajectory (default: None)
  --mDEG             Add "--mTop_markers" to allow filteration (default:
                     False)
  --mGCN             Add "--mGCN" to allow gene co-expression analysis
                     (default: False)
  --mTop_markers     Add "--mTop_markers" to allow top markers analysis
                     (default: False)

output arguments:
  --out OUT          Directory to save file
  --prefix PREFIX    Sample ID, will be used as output prefix and seuratObj
                     ident

```

### ScAT CellCellInteraction
```bash
> ${scat_dir}/bin/ScAT CellCellInteraction -h
usage: ScAT CellCellInteraction [-h] --Rscript RSCRIPT --input INPUT
                                [--species {Human,Mouse}] --meta META --out
                                OUT [--prefix PREFIX]

optional arguments:
  -h, --help            show this help message and exit

input arguments:
  --Rscript RSCRIPT     Rscript path (default: None)
  --input INPUT         Input raw data file (*.txt.gz)
  --species {Human,Mouse}
                        Species name (default: Human)
  --meta META           Cell labels (Recommand the cell types results
                        generated by Annotation module (default: None)

output arguments:
  --out OUT             Directory to save file
  --prefix PREFIX       Sample ID, will be used as output prefix and seuratObj
                        ident

```

### ScAT TFregulon
```bash
> ${scat_dir}/bin/ScAT TFregulon -h
usage: ScAT TFregulon [-h] --GeneMatrix GENEMATRIX --species SPECIES
                      [--annot ANNOT] [--TFDataset TFDATASET] [--TFName TF]
                      [--NetFile NET] [--SCENICloom SCENICLOOM] [--h5ad H5AD]
                      [--TFreglist TFREGLIST] [--cpu CPU] --outdir OUTDIR

optional arguments:
  -h, --help            show this help message and exit

input arguments:
  --GeneMatrix GENEMATRIX
                        Gene expression matrix of single cell.
  --species SPECIES     Human or mouse.
  --annot ANNOT         Cell type annotation file.
  --TFDataset TFDATASET
                        TF dataset.
  --TFName TF           Provide a TF name.
  --NetFile NET         Co-expression network inferred by pyscenic named
                        '02.Network.Adj.tsv' in the outdir.
  --SCENICloom SCENICLOOM
                        Loom file outputed by SCENIC named
                        'SCENIC.Regulon.AUC.loom '.
  --h5ad H5AD           h5ad file including AdataExpression and regulonAUC
                        value info named 'Adata.Expression.Regulon.AUC.h5ad'.
  --TFreglist TFREGLIST
                        TF regulons list file;One column.
  --cpu CPU             Cpu numbers.

output arguments:
  --outdir OUTDIR       Output result directory

```

### ScAT SpatialCluster
```bash
> ${scat_dir}/bin/ScAT SpatialCluster -h
usage: ScAT SpatialCluster [-h] --Rscript RSCRIPT --Stdata STDATA
                           [--clusterNumber CLUSTERNUMBER] --iterations
                           ITERATIONS --out OUT

optional arguments:
  -h, --help            show this help message and exit

input arguments:
  --Rscript RSCRIPT     Rscript path (default: None)
  --Stdata STDATA       Seurat rds before clustering...
  --clusterNumber CLUSTERNUMBER
                        Set cluster number of BayesSpace spatial cluster
  --iterations ITERATIONS
                        iterations number of MCMC

output arguments:
  --out OUT             Directory to save file
  
```

### ScAT SpatialPatternGene
```bash
> ${scat_dir}/bin/ScAT SpatialPatternGene -h
usage: ScAT SpatialPatternGene [-h] --Rscript RSCRIPT --input INPUT --out OUT
                               [--prefix PREFIX]

optional arguments:
  -h, --help         show this help message and exit

input arguments:
  --Rscript RSCRIPT  Rscript path.(default: None)
  --input INPUT      Seurat rds.

output arguments:
  --out OUT          Directory to save file
  --prefix PREFIX    sample ID, will be used as output file prefix (default:
                     sample)

```

### ScAT DeconvSpot
```bash
> ${scat_dir}/bin/ScAT DeconvSpot -h
usage: ScAT DeconvSpot [-h] --Rscript RSCRIPT --sc_rds SC_RDS --st_rds ST_RDS
                       [--piescale PIESCALE] [--ClusterFile CLUSTERFILE] --out
                       OUT

optional arguments:
  -h, --help            show this help message and exit

input arguments:
  --Rscript RSCRIPT     Rscript path (default: None)
  --sc_rds SC_RDS       Input seurate rds single-cell data after cell cluster
                        annotation.
  --st_rds ST_RDS       Input seurate rds spatial transcriptome data.
  --piescale PIESCALE   Pie spot size.
  --ClusterFile CLUSTERFILE
                        View a set of specific cluster locations.

output arguments:
  --out OUT             Directory to save file

```


```python

```
