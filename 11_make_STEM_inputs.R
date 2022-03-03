###Script for normalizing metabolites for STEM filtering

argg <- commandArgs(T)

if (length(argg) != 4){
  stop('ARGS: 1) PeakArea file (with only MSMS and adducts removed) 
       2) POS or NEG? 3) Leaf, Root or Both?
       4) Experiment name (for ex, All_exps)')
}

########Function##########

get_groups <- function(df){
  ct <- colnames(df)
  groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  return(groups)
}

#######Body##########

##getting zscores of values
infil <- read.table(argg[1], sep = '\t', stringsAsFactors = F, quote = "", fill = NA, header = T)

toscl <- t(infil[,33:ncol(infil)])
zscrs <- scale(toscl)
goodzscrs <- t(zscrs)

##now getting average of the groups
groups <- get_groups(goodzscrs)
out <- c()

for (grp in unique(groups)){
  colstoavg <- goodzscrs[,grepl(grp, colnames(goodzscrs), fixed = T)]
  avgrg <- rowMeans(colstoavg)
  out <- cbind(out, avgrg)
}

finalout <- cbind(infil[,1], out)
colnames(finalout) <- c('SPOT', unique(groups))

write.table(finalout, file = paste0('E:/MS_Data/BrachyMetabolites/STEM/', argg[4], '/STEM_input_',argg[2], argg[3], '.txt'),
            sep = '\t', row.names = F, quote = F)

print('Done!')