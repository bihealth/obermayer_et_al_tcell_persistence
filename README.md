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

## session info

**R version 4.0.3 (2020-10-10)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=en_US.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=en_US.UTF-8_, _LC_COLLATE=en_US.UTF-8_, _LC_MONETARY=en_US.UTF-8_, _LC_MESSAGES=en_US.UTF-8_, _LC_PAPER=en_US.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=en_US.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 
_grid_, _parallel_, _stats4_, _stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
_pander(v.0.6.5)_, _org.Hs.eg.db(v.3.12.0)_, _scDiffCom(v.0.1.0)_, _forcats(v.0.5.1)_, _ggfortify(v.0.4.11)_, _ggExtra(v.0.9)_, _caret(v.6.0-90)_, _lattice(v.0.20-41)_, _randomForest(v.4.6-14)_, _ggpubr(v.0.4.0)_, _pheatmap(v.1.0.12)_, _sctransform(v.0.3.3)_, _DT(v.0.17)_, _lme4(v.1.1-27.1)_, _gtools(v.3.9.2)_, _AnnotationDbi(v.1.52.0)_, _RColorBrewer(v.1.1-2)_, _Matrix(v.1.3-4)_, _igraph(v.1.2.11)_, _stringr(v.1.4.0)_, _ggforce(v.0.3.2.9000)_, _ggraph(v.2.0.5)_, _tmod(v.0.46.2)_, _Binarize(v.1.3)_, _diptest(v.0.76-0)_, _dendextend(v.1.15.1)_, _ComplexHeatmap(v.2.6.2)_, _DESeq2(v.1.30.1)_, _SummarizedExperiment(v.1.20.0)_, _Biobase(v.2.50.0)_, _MatrixGenerics(v.1.2.1)_, _matrixStats(v.0.61.0)_, _GenomicRanges(v.1.42.0)_, _GenomeInfoDb(v.1.26.7)_, _IRanges(v.2.24.1)_, _S4Vectors(v.0.28.1)_, _BiocGenerics(v.0.36.1)_, _plyr(v.1.8.6)_, _ggalluvial(v.0.12.3)_, _fitdistrplus(v.1.1-8)_, _survival(v.3.2-7)_, _MASS(v.7.3-53)_, _cluster(v.2.1.0)_, _tidyr(v.1.2.0)_, _scRepertoire(v.1.1.2)_, _dplyr(v.1.0.8)_, _ggrepel(v.0.9.1)_, _cowplot(v.1.1.1)_, _ggplot2(v.3.3.5)_, _SeuratDisk(v.0.0.0.9018)_, _SeuratObject(v.4.0.4)_ and _Seurat(v.4.1.0)_

**loaded via a namespace (and not attached):** 
_SparseM(v.1.81)_, _scattermore(v.0.8)_, _ModelMetrics(v.1.2.2.2)_, _evmix(v.2.12)_, _bit64(v.4.0.5)_, _knitr(v.1.31)_, _irlba(v.2.3.5)_, _DelayedArray(v.0.16.3)_, _data.table(v.1.14.2)_, _rpart(v.4.1-15)_, _RCurl(v.1.98-1.3)_, _doParallel(v.1.0.16)_, _generics(v.0.1.2)_, _RSQLite(v.2.2.6)_, _RANN(v.2.6.1)_, _VGAM(v.1.1-5)_, _future(v.1.24.0)_, _bit(v.4.0.4)_, _lubridate(v.1.7.10)_, _spatstat.data(v.2.1-2)_, _httpuv(v.1.6.5)_, _isoband(v.0.2.5)_, _assertthat(v.0.2.1)_, _viridis(v.0.5.1)_, _gower(v.0.2.2)_, _xfun(v.0.22)_, _hms(v.1.0.0)_, _plotwidgets(v.0.4)_, _promises(v.1.2.0.1)_, _fansi(v.1.0.3)_, _readxl(v.1.3.1)_, _caTools(v.1.18.2)_, _DBI(v.1.1.1)_, _geneplotter(v.1.68.0)_, _htmlwidgets(v.1.5.4)_, _powerTCR(v.1.10.3)_, _spatstat.geom(v.2.3-2)_, _purrr(v.0.3.4)_, _ellipsis(v.0.3.2)_, _backports(v.1.2.1)_, _permute(v.0.9-5)_, _annotate(v.1.68.0)_, _deldir(v.1.0-6)_, _vctrs(v.0.3.8)_, _Cairo(v.1.5-12.2)_, _ROCR(v.1.0-11)_, _abind(v.1.4-5)_, _cachem(v.1.0.6)_, _withr(v.2.5.0)_, _vegan(v.2.5-7)_, _goftest(v.1.2-3)_, _gsl(v.2.1-6)_, _lazyeval(v.0.2.2)_, _crayon(v.1.5.0)_, _genefilter(v.1.72.1)_, _hdf5r(v.1.3.3)_, _labeling(v.0.4.2)_, _recipes(v.0.1.17)_, _pkgconfig(v.2.0.3)_, _tweenr(v.1.0.2)_, _nlme(v.3.1-151)_, _nnet(v.7.3-14)_, _rlang(v.1.0.2)_, _globals(v.0.14.0)_, _lifecycle(v.1.0.1)_, _miniUI(v.0.1.1.1)_, _cellranger(v.1.1.0)_, _polyclip(v.1.10-0)_, _lmtest(v.0.9-40)_, _carData(v.3.0-4)_, _boot(v.1.3-25)_, _zoo(v.1.8-9)_, _beeswarm(v.0.4.0)_, _ggridges(v.0.5.3)_, _GlobalOptions(v.0.1.2)_, _png(v.0.1-7)_, _viridisLite(v.0.4.0)_, _rjson(v.0.2.20)_, _bitops(v.1.0-7)_, _pROC(v.1.18.0)_, _KernSmooth(v.2.23-18)_, _Biostrings(v.2.58.0)_, _blob(v.1.2.1)_, _shape(v.1.4.5)_, _parallelly(v.1.30.0)_, _spatstat.random(v.2.1-0)_, _rstatix(v.0.7.0)_, _ggsignif(v.0.6.2)_, _scales(v.1.2.0)_, _memoise(v.2.0.0)_, _magrittr(v.2.0.2)_, _ica(v.1.0-2)_, _gplots(v.3.1.1)_, _zlibbioc(v.1.36.0)_, _compiler(v.4.0.3)_, _clue(v.0.3-58)_, _cli(v.3.2.0)_, _XVector(v.0.30.0)_, _listenv(v.0.8.0)_, _patchwork(v.1.1.1)_, _pbapply(v.1.5-0)_, _mgcv(v.1.8-33)_, _tidyselect(v.1.1.2)_, _stringi(v.1.7.6)_, _locfit(v.1.5-9.4)_, _tools(v.4.0.3)_, _future.apply(v.1.8.1)_, _rio(v.0.5.27)_, _circlize(v.0.4.12)_, _rstudioapi(v.0.13)_, _foreach(v.1.5.1)_, _foreign(v.0.8-81)_, _tagcloud(v.0.6)_, _gridExtra(v.2.3)_, _cubature(v.2.0.4.1)_, _prodlim(v.2019.11.13)_, _farver(v.2.1.0)_, _Rtsne(v.0.15)_, _digest(v.0.6.29)_, _lava(v.1.6.10)_, _shiny(v.1.7.1)_, _Rcpp(v.1.0.8.3)_, _car(v.3.0-11)_, _broom(v.0.7.9)_, _later(v.1.3.0)_, _RcppAnnoy(v.0.0.19)_, _httr(v.1.4.2)_, _colorspace(v.2.0-3)_, _XML(v.3.99-0.6)_, _tensor(v.1.5)_, _reticulate(v.1.24)_, _splines(v.4.0.3)_, _uwot(v.0.1.11)_, _spatstat.utils(v.2.3-0)_, _graphlayouts(v.0.7.1)_, _plotly(v.4.10.0)_, _xtable(v.1.8-4)_, _jsonlite(v.1.8.0)_, _nloptr(v.1.2.2.2)_, _truncdist(v.1.0-2)_, _tidygraph(v.1.2.0)_, _timeDate(v.3043.102)_, _ipred(v.0.9-12)_, _R6(v.2.5.1)_, _pillar(v.1.7.0)_, _htmltools(v.0.5.2)_, _mime(v.0.12)_, _glue(v.1.6.2)_, _fastmap(v.1.1.0)_, _minqa(v.1.2.4)_, _colorDF(v.0.1.4)_, _BiocParallel(v.1.24.1)_, _class(v.7.3-17)_, _codetools(v.0.2-18)_, _utf8(v.1.2.2)_, _spatstat.sparse(v.2.1-0)_, _tibble(v.3.1.6)_, _evd(v.2.3-3)_, _curl(v.4.3.2)_, _leiden(v.0.3.9)_, _zip(v.2.1.1)_, _openxlsx(v.4.2.3)_, _munsell(v.0.5.0)_, _GetoptLong(v.1.0.5)_, _GenomeInfoDbData(v.1.2.4)_, _iterators(v.1.0.13)_, _haven(v.2.4.3)_, _reshape2(v.1.4.4)_, _gtable(v.0.3.0)_ and _spatstat.core(v.2.4-0)_
