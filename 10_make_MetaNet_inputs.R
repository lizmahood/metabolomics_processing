##For making PCC network with MetaNetter
##Input is filtered peak area file


########Function##########

get_groups <- function(df){
  ct <- colnames(df)
  groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  return(groups)
}

#######Body##########
argg <- commandArgs(T)

if (length(argg) != 4){
  stop('ARGS: 1) input file (PeakArea) 2) POS or NEG? 
      3) Leaf OR Root OR Both
      4) Experiment name (for ex, All_exps)')
}

infil <- read.table(argg[1], sep = '\t', stringsAsFactors = F, header = T, quote = "", fill = NA)

##remove non fragmented metabolites
frags <- infil[which(infil$MS.MS.spectrum != ''),]

##remove adducts (for NEG)
if (argg[2] == 'NEG'){
  frags <- frags[-which(frags$Adduct.type == '[M+FA-H]-'),]
  frags <- frags[-which(frags$Adduct.type == '[M-H2O-H]-' & grepl('[M-H]-', frags$Post.curation.result, fixed = T)),]
}

##remove adducts (for POS)
if (argg[2] == 'POS'){
  frags <- frags[-which(frags$Adduct.type != '[M+H]+'),]
}

##ordering columns
#frags[,c(33:ncol(frags))] <- frags[,order(colnames(frags)[33:ncol(frags)])]

##writing this out
write.table(frags, file = paste0('E:/MS_Data/BrachyMetabolites/MetaNetterClustering/', argg[4], '/',
                                 basename(argg[1]), '_', argg[3], '_msms.tab'), row.names = F, quote = F, sep = '\t')


##now getting the averages of each group, deleting rest of columns
frags <- frags[,-c(1,2,4:32)]
groups <- get_groups(frags[,2:ncol(frags)])

out <- c()

for (grp in unique(groups)){
  print(grp)
  colstoavg <- frags[,grepl(grp, colnames(frags), fixed = T)]
  avgrg <- rowMeans(colstoavg)
  out <- cbind(out, avgrg)
}

finalout <- cbind(frags[,1], out)
colnames(finalout) <- c('mass', unique(groups))

write.table(finalout, sep = '\t', file = paste0('E:/MS_Data/BrachyMetabolites/MetaNetterClustering/', argg[4], '/',
                                                basename(argg[1]), '_', argg[3], '_INPT.txt'), col.names = T, row.names = F, quote = F)

print('Done!')
