rnaseq.dirs <- list.files('data')
print(rnaseq.dirs)
rnaseq.dirs[1]

list.files paste0(list.files(paste0('data/', rnaseq.dirs[1]), pattern='.tsv'))
#raw_counts <- paste0(list.files(paste0('data/', rnaseq.dirs[1]), pattern='.tsv'))
                     
raw_counts <- read.delim(paste0('./data/', rnaseq.dirs[1],'/', list.files(paste0('data/', rnaseq.dirs[1]), pattern='.tsv')), skip=1)[-c(1:4), ]
raw_counts <-raw_counts [, c(1,4)]
head(raw_counts)
colnames(raw_counts)[2] <-rnaseq.dirs[1]

c=3
for (i in rnased.dirs [-1]){
    tmp <- read.delim(paste0('./data/rnaseq', i, '/', list.files(paste0('data/rnaseq', i)
    pattern= '.tsv')), skip=1) [-c(1:4), ]
    tmp <-tmp[,4]
    raw_counts <-cbind(raw_counts, tmp)
    colnames(raw_count)[c] <- i
    c= c+1
    print(c)
}
rm(tmp)
head(raw_counts)
dim(raw_counts)
rownnames(raw_counts) <- raw_counts$gene_idraw_counts <-as.matrix(raw_counts[,-1])
class(raw_counts)
is.numeric(raw_counts)
raw_counts[1:4, 1:4] #salvar em um diretório (mkdir data/processed)
saveRDS(raw_counts, file= 'data/processed_data/TCGA_COAD_RNASEQ.rds')

#para fazer a matriz deixar só valores numéricos
#unstranded, gene_id, gene_name 

for (i in rnased.dirs [-1]){
    tmp <- read.delim(paste0('./data/rnaseq', i, '/', list.files(paste0('data/rnaseq', i)
    pattern= '.tsv')), skip=1) [-c(1:4), ]
    tmp <-tmp[,4]
    raw_counts <-cbind(raw_counts, tmp)
    colnames(raw_count)[c] <- i
    c= c+1
    print(c)
}