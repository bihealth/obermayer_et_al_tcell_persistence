# Code repository for Obermayer et al. "Identification of transcriptional attributes associated with persistence of T-cells in allogeneic hematopoietic stem cell transplantation"

## data access

- processed data is available upon request or (pending data privacy review) from NCBI GEO
- `pbmc_multimodal.h5seurat` is available [here](https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat) 
- MSigDB is available [here](http://www.gsea-msigdb.org/gsea/downloads.jsp)

input data structure should looks like so:

```
data
├── bulk_TCR
│   └── DE07NGSUKBD128950_ct.tsv.gz 
├── cellranger
│   ├── D10postPBMC_filtered_contigs.csv 
│   ├── D10postPBMC_filtered_counts.h5 
│   ├── D10prePBMC_filtered_contigs.csv 
│   ├── D10prePBMC_filtered_counts.h5 
│   ├── ...
├── DE
├── scDiffCom
├── seurat
│   └── pbmc_multimodal.h5seurat
├── VDJdb
│   ├── input
│   ├── TRA
│   └── TRB
└── vireo
    └── D31postR24d30PBMC_donor_ids.tsv
```

## Seurat processing

run `R/processing.Rmd`

## scDiffCom

use `R/scDiffCom.R` to run [scDiffCom](https://github.com/CyrilLagger/scDiffCom)

## VDJdb

VDJmatch is available from [here](https://github.com/antigenomics/vdjmatch) and should be run like so

```
Rscript R/prepare_vdjdb.R
bash scripts/run_vdjdb.sh
```

## effectorness (=pseudotime) score

use `R/pseudotime.Rmd`

## paper figures

all code for paper figures is in `R/paper_figures.Rmd`
