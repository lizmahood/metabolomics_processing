### For splitting up DAM files into up vs. down accumulated
### Also for putting canopus annotations onto the DAMs

argg <- commandArgs(T)

if (length(argg) != 3){
  stop('ARGS: 1) path to DAMs 2) path to canopus_summary filtered
       3) desired output file')
}

get_up_vs_down <- function(damfiles){
  olist <- list()
  for (fil in damfiles){
    hmm <- unlist(strsplit(fil, '_diff_'))[2]
    cond <- gsub('.tab', '', hmm, fixed = T)
    fildf <- read.table(fil, header = T, sep = '\t', fill = NA, quote = "")
    condup <- fildf[fildf$nfold_change > 0,]
    conddown <- fildf[fildf$nfold_change < 0,]
    olist[[cond]] <- list(condup, conddown)
  }
  return(olist)
}

merge_canopus <- function(condlist, canopus, odir){
  canopus_to_merge <- canopus[, c('name', 'class', 'superclass', 'subclass', 'level.5')]
  i <- 1                            
  for (cond in condlist){
    up <- cond[[1]]
    down <- cond[[2]]
    canopus_up <- merge(up, canopus_to_merge, all.x = T, by.x = 'Alignment.ID', by.y = 'name')
    canopus_down <- merge(down, canopus_to_merge, all.x = T, by.x = 'Alignment.ID',
                          by.y = 'name')
    canopus_up[is.na(canopus_up)] <- 'None'
    canopus_down[is.na(canopus_down)] <- 'None'
    write.table(canopus_up, file = paste0(odir, '_', names(condlist)[i], '_upDAM_canopus.tsv'),
                sep = '\t', quote = F, row.names = F)
    write.table(canopus_down, file = paste0(odir, '_', names(condlist)[i], '_downDAM_canopus.tsv'),
                sep = '\t', quote = F, row.names = F)
    i <- i + 1
  }
}

damfiles <- list.files(argg[1], pattern = 'FDR', full.names = T)
canopus <- read.table(argg[2], header = T, fill = NA, quote = "", sep = '\t')

condlist <- get_up_vs_down(damfiles)
merge_canopus(condlist, canopus, argg[3])
print('Done!')
