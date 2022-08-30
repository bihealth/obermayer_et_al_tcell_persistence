library(Seurat)
library(dplyr)

pbmc.all <- readRDS(file.path('data','seurat','pbmc_all.rds'))

for (lib in unique(pbmc.all$orig.ident)) {
  infile <- file.path('data','cellranger',paste0(lib,'_filtered_contigs.csv'))
  for (sample in unique(pbmc.all$sample[pbmc.all$orig.ident==lib])) {
    bcs <- gsub('.*_','',Cells(pbmc.all)[(pbmc.all$sample==sample) & (pbmc.all$predicted.celltype.l1 %in% c('CD4 T','CD8 T','other T'))])
    if (file.exists(infile)) {
      tmp <- read.csv(infile,
                      stringsAsFactors=FALSE) %>%
        dplyr::mutate(is_cell=is_cell=='true',
                      high_confidence=high_confidence=='true',
                      full_length=full_length=='true',
                      productive=productive=='true') %>%
        dplyr::filter(barcode %in% bcs) %>%
        dplyr::mutate(count=umis,
                      freq=umis/sum(umis),
                      cdr3nt=cdr3_nt,
                      cdr3aa=cdr3,
                      v=v_gene,
                      d=d_gene,
                      j=j_gene) %>%
        dplyr::mutate(d=ifelse(d=='',paste0(chain,'D'),d)) %>%
        dplyr::group_by(cdr3nt,cdr3aa) %>%
        dplyr::summarise(count=sum(count),
                         freq=sum(freq),
                         cdr3nt=dplyr::first(cdr3nt),
                         cdr3aa=dplyr::first(cdr3aa),
                         v=dplyr::first(v),
                         d=dplyr::first(d),
                         j=dplyr::first(j),
                         chain=dplyr::first(chain)) %>%
        dplyr::select(count,freq,cdr3nt,cdr3aa,v,d,j,chain) 
      write.table(tmp %>% dplyr::filter(chain=='TRA') %>% dplyr::select(-chain),
                  file.path('data','VDJdb','input',paste0(sample,'_TRA.txt')),sep='\t',row.names=FALSE,quote=FALSE)
      write.table(tmp %>% dplyr::filter(chain=='TRB') %>% dplyr::select(-chain),
                  file.path('data','VDJdb','input',paste0(sample,'_TRB.txt')),sep='\t',row.names=FALSE,quote=FALSE)
    }
  }
}
