# Code repository for Obermayer et al. "Identification of transcriptional attributes associated with persistence of T-cells in allogeneic hematopoietic stem cell transplantation"

## data access

- processed data is available from NCBI GEO under accession [GSE222633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222633)
- `pbmc_multimodal.h5seurat` is available [here](https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat) 
- MSigDB is available [here](http://www.gsea-msigdb.org/gsea/downloads.jsp)

input data structure should looks like so (after stripping `GSE222633_` prefixes from files):

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
├── scDiffCom
├── seurat
│   └── pbmc_multimodal.h5seurat
├── tables
└── VDJdb
    ├── input
    ├── TRA
    └── TRB
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

## session info

**R version 4.0.3 (2020-10-10)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=en_US.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=en_US.UTF-8_, _LC_COLLATE=en_US.UTF-8_, _LC_MONETARY=en_US.UTF-8_, _LC_MESSAGES=en_US.UTF-8_, _LC_PAPER=en_US.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=en_US.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 

* grid 
* parallel 
* stats4 
* stats 
* graphics 
* grDevices 
* utils 
* datasets 
* methods 
* base 


**other attached packages:** 

* pander(v.0.6.5) 
* org.Hs.eg.db(v.3.12.0) 
* scDiffCom(v.0.1.0) 
* forcats(v.0.5.1) 
* ggfortify(v.0.4.11) 
* ggExtra(v.0.9) 
* caret(v.6.0-90) 
* lattice(v.0.20-41) 
* randomForest(v.4.6-14) 
* ggpubr(v.0.4.0) 
* pheatmap(v.1.0.12) 
* sctransform(v.0.3.3) 
* DT(v.0.17) 
* lme4(v.1.1-27.1) 
* gtools(v.3.9.2) 
* AnnotationDbi(v.1.52.0) 
* RColorBrewer(v.1.1-2) 
* Matrix(v.1.3-4) 
* igraph(v.1.2.11) 
* stringr(v.1.4.0) 
* ggforce(v.0.3.2.9000) 
* ggraph(v.2.0.5) 
* tmod(v.0.46.2) 
* Binarize(v.1.3) 
* diptest(v.0.76-0) 
* dendextend(v.1.15.1) 
* ComplexHeatmap(v.2.6.2) 
* DESeq2(v.1.30.1) 
* SummarizedExperiment(v.1.20.0) 
* Biobase(v.2.50.0) 
* MatrixGenerics(v.1.2.1) 
* matrixStats(v.0.61.0) 
* GenomicRanges(v.1.42.0) 
* GenomeInfoDb(v.1.26.7) 
* IRanges(v.2.24.1) 
* S4Vectors(v.0.28.1) 
* BiocGenerics(v.0.36.1) 
* plyr(v.1.8.6) 
* ggalluvial(v.0.12.3) 
* fitdistrplus(v.1.1-8) 
* survival(v.3.2-7) 
* MASS(v.7.3-53) 
* cluster(v.2.1.0) 
* tidyr(v.1.2.0) 
* scRepertoire(v.1.1.2) 
* dplyr(v.1.0.8) 
* ggrepel(v.0.9.1) 
* cowplot(v.1.1.1) 
* ggplot2(v.3.3.5) 
* SeuratDisk(v.0.0.0.9018) 
* SeuratObject(v.4.0.4) 
* Seurat(v.4.1.0) 


**loaded via a namespace (and not attached):** 

* SparseM(v.1.81) 
* scattermore(v.0.8) 
* ModelMetrics(v.1.2.2.2) 
* evmix(v.2.12) 
* bit64(v.4.0.5) 
* knitr(v.1.31) 
* irlba(v.2.3.5) 
* DelayedArray(v.0.16.3) 
* data.table(v.1.14.2) 
* rpart(v.4.1-15) 
* RCurl(v.1.98-1.3) 
* doParallel(v.1.0.16) 
* generics(v.0.1.2) 
* RSQLite(v.2.2.6) 
* RANN(v.2.6.1) 
* VGAM(v.1.1-5) 
* future(v.1.24.0) 
* bit(v.4.0.4) 
* lubridate(v.1.7.10) 
* spatstat.data(v.2.1-2) 
* httpuv(v.1.6.5) 
* isoband(v.0.2.5) 
* assertthat(v.0.2.1) 
* viridis(v.0.5.1) 
* gower(v.0.2.2) 
* xfun(v.0.22) 
* hms(v.1.0.0) 
* plotwidgets(v.0.4) 
* promises(v.1.2.0.1) 
* fansi(v.1.0.3) 
* readxl(v.1.3.1) 
* caTools(v.1.18.2) 
* DBI(v.1.1.1) 
* geneplotter(v.1.68.0) 
* htmlwidgets(v.1.5.4) 
* powerTCR(v.1.10.3) 
* spatstat.geom(v.2.3-2) 
* purrr(v.0.3.4) 
* ellipsis(v.0.3.2) 
* backports(v.1.2.1) 
* permute(v.0.9-5) 
* annotate(v.1.68.0) 
* deldir(v.1.0-6) 
* vctrs(v.0.3.8) 
* Cairo(v.1.5-12.2) 
* ROCR(v.1.0-11) 
* abind(v.1.4-5) 
* cachem(v.1.0.6) 
* withr(v.2.5.0) 
* vegan(v.2.5-7) 
* goftest(v.1.2-3) 
* gsl(v.2.1-6) 
* lazyeval(v.0.2.2) 
* crayon(v.1.5.0) 
* genefilter(v.1.72.1) 
* hdf5r(v.1.3.3) 
* labeling(v.0.4.2) 
* recipes(v.0.1.17) 
* pkgconfig(v.2.0.3) 
* tweenr(v.1.0.2) 
* nlme(v.3.1-151) 
* nnet(v.7.3-14) 
* rlang(v.1.0.2) 
* globals(v.0.14.0) 
* lifecycle(v.1.0.1) 
* miniUI(v.0.1.1.1) 
* cellranger(v.1.1.0) 
* polyclip(v.1.10-0) 
* lmtest(v.0.9-40) 
* carData(v.3.0-4) 
* boot(v.1.3-25) 
* zoo(v.1.8-9) 
* beeswarm(v.0.4.0) 
* ggridges(v.0.5.3) 
* GlobalOptions(v.0.1.2) 
* png(v.0.1-7) 
* viridisLite(v.0.4.0) 
* rjson(v.0.2.20) 
* bitops(v.1.0-7) 
* pROC(v.1.18.0) 
* KernSmooth(v.2.23-18) 
* Biostrings(v.2.58.0) 
* blob(v.1.2.1) 
* shape(v.1.4.5) 
* parallelly(v.1.30.0) 
* spatstat.random(v.2.1-0) 
* rstatix(v.0.7.0) 
* ggsignif(v.0.6.2) 
* scales(v.1.2.0) 
* memoise(v.2.0.0) 
* magrittr(v.2.0.2) 
* ica(v.1.0-2) 
* gplots(v.3.1.1) 
* zlibbioc(v.1.36.0) 
* compiler(v.4.0.3) 
* clue(v.0.3-58) 
* cli(v.3.2.0) 
* XVector(v.0.30.0) 
* listenv(v.0.8.0) 
* patchwork(v.1.1.1) 
* pbapply(v.1.5-0) 
* mgcv(v.1.8-33) 
* tidyselect(v.1.1.2) 
* stringi(v.1.7.6) 
* locfit(v.1.5-9.4) 
* tools(v.4.0.3) 
* future.apply(v.1.8.1) 
* rio(v.0.5.27) 
* circlize(v.0.4.12) 
* rstudioapi(v.0.13) 
* foreach(v.1.5.1) 
* foreign(v.0.8-81) 
* tagcloud(v.0.6) 
* gridExtra(v.2.3) 
* cubature(v.2.0.4.1) 
* prodlim(v.2019.11.13) 
* farver(v.2.1.0) 
* Rtsne(v.0.15) 
* digest(v.0.6.29) 
* lava(v.1.6.10) 
* shiny(v.1.7.1) 
* Rcpp(v.1.0.8.3) 
* car(v.3.0-11) 
* broom(v.0.7.9) 
* later(v.1.3.0) 
* RcppAnnoy(v.0.0.19) 
* httr(v.1.4.2) 
* colorspace(v.2.0-3) 
* XML(v.3.99-0.6) 
* tensor(v.1.5) 
* reticulate(v.1.24) 
* splines(v.4.0.3) 
* uwot(v.0.1.11) 
* spatstat.utils(v.2.3-0) 
* graphlayouts(v.0.7.1) 
* plotly(v.4.10.0) 
* xtable(v.1.8-4) 
* jsonlite(v.1.8.0) 
* nloptr(v.1.2.2.2) 
* truncdist(v.1.0-2) 
* tidygraph(v.1.2.0) 
* timeDate(v.3043.102) 
* ipred(v.0.9-12) 
* R6(v.2.5.1) 
* pillar(v.1.7.0) 
* htmltools(v.0.5.2) 
* mime(v.0.12) 
* glue(v.1.6.2) 
* fastmap(v.1.1.0) 
* minqa(v.1.2.4) 
* colorDF(v.0.1.4) 
* BiocParallel(v.1.24.1) 
* class(v.7.3-17) 
* codetools(v.0.2-18) 
* utf8(v.1.2.2) 
* spatstat.sparse(v.2.1-0) 
* tibble(v.3.1.6) 
* evd(v.2.3-3) 
* curl(v.4.3.2) 
* leiden(v.0.3.9) 
* zip(v.2.1.1) 
* openxlsx(v.4.2.3) 
* munsell(v.0.5.0) 
* GetoptLong(v.1.0.5) 
* GenomeInfoDbData(v.1.2.4) 
* iterators(v.1.0.13) 
* haven(v.2.4.3) 
* reshape2(v.1.4.4) 
* gtable(v.0.3.0) 
* spatstat.core(v.2.4-0) 
