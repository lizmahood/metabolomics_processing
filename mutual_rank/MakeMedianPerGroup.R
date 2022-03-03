
get_groups <- function(df, typ){
  if (typ == 'trans'){
    ct <- colnames(df)
    groups <- substr(ct, 1, nchar(ct)-2)
  }else if (typ == 'metab'){
    ct <- colnames(df[,33:ncol(df)])
    groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  }
  return(groups)
}

match_metab_with_trans <- function(tdf, mdf, outliers){
  newm <- mdf[,33:ncol(mdf)]
  newc <- unlist(lapply(strsplit(colnames(newm), '_Run'), '[', 1))
  colnames(newm) <- newc
  toremove <- which(newc %in% outliers)
  newm <- newm[,-toremove]
  
  ## now for matching between trans and metab
  notrans <- which(!(colnames(newm) %in% colnames(tdf)))
  newm <- newm[,-notrans]
  nometabs <- which(!(colnames(tdf) %in% colnames(newm)))
  if (length(nometabs) > 0){
    tdf <- tdf[,-nometabs]
  }
  ordertdf <- tdf[,order(colnames(tdf))]
  ordermdf <- newm[,order(colnames(newm))]
  out <- rbind.data.frame(ordertdf, ordermdf)
  return(out)
}

argg <- commandArgs(T)

if (length(argg) != 4){
  stop('ARGS: Transformed transcript file 2) Normalized metabolite file
       3) Outliers, separated by commas (ex Hydro.Cop.Root_2,Sym.SporeW.Root_3)
       4) Brachy metabolic genes file')
}

trans_rlog <- read.table(argg[1], sep = '\t', header = T, quote = "", fill = NA, row.names = 1)
metabs <- read.table(argg[2], sep = '\t', header = T, quote = "", fill = NA, row.names = 1)
outliers <- unlist(strsplit(argg[3], ','))
metabolic_genes <- read.table(argg[4], sep = '\t', header = T, quote = "", fill = NA)

## Expression file for only metabolic genes
trans_rlog_metabolic <- trans_rlog[which(rownames(trans_rlog) %in% metabolic_genes$genes),]

combtransmetab <- match_metab_with_trans(trans_rlog, metabs, outliers)
comb_metabolic <- match_metab_with_trans(trans_rlog_metabolic, metabs, outliers)
grps <- get_groups(combtransmetab, 'trans')

outl <- list()
outlm <- list() ## for metabolic genes only

for (grp in unique(grps)){
  tmp <- combtransmetab[,which(grepl(grp, colnames(combtransmetab)))]
  outl[[grp]] <- matrixStats::rowMedians(as.matrix(tmp))
  
  tmpm <- comb_metabolic[,which(grepl(grp, colnames(comb_metabolic)))]
  outlm[[grp]] <- matrixStats::rowMedians(as.matrix(tmpm))
}

brachy_rlogmed <- do.call(cbind.data.frame, outl)
row.names(brachy_rlogmed) <- row.names(combtransmetab)
write.table(brachy_rlogmed, file = paste0(argg[1], 'metab_trans_median.tab'),
            sep = '\t', quote = F, row.names = T)

## File for metabolic genes only
rlog_metabolic <- do.call(cbind.data.frame, outlm)
row.names(rlog_metabolic) <- row.names(comb_metabolic)
write.table(rlog_metabolic, file = paste0(argg[1], 'metab_trans_median_metabolicgenes.tab'),
            sep = '\t', quote = F, row.names = T)

print('Done!')
