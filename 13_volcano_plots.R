##script for graphing. Makes volcano plots from tissue specific metabolites
library(dplyr)
library(ggplot2)
library(tidyverse)

mk_fil <- function(indir){
  ##indir is a string, path to directory with Leaf and Root subdirs
  ##RETURNS: list with each tissue-cond pair as a dataframe
  tissues <- c('Leaf', 'Root')
  
  filelist <- list()
  for (tis in tissues){
    readdir <- paste(indir, tis, 'MSDIAL_out', sep = '/')
    print(readdir)
    fils <- list.files(readdir, pattern = 'all_pval_fc', full.names = T)
    if (length(fils) <1){
      print('wow')
      readdir <- paste0(indir, '/MSDIAL_out')
      fils <- list.files(readdir, pattern = 'all_pval_fc', full.names = T)
    }
    print(fils)
    for (fil in fils){
      temp <- read.table(fil, sep = '\t', header = T, quote = "", fill = NA)
      thistis <- rep(tis, nrow(temp))
      ##python has spoiled me
      tcond <- strsplit(fil, '_', fixed = T)[[1]][length(strsplit(fil, '_', fixed = T)[[1]])]
      print(tcond)
      cond <- strsplit(tcond, '.', fixed = T)[[1]][1]
      thiscond <- rep(cond, nrow(temp))
      temp <- cbind(temp, thistis, thiscond)
      filelist <- append(filelist, list(temp))
    }
  }
  return(filelist)
}

argg <- commandArgs(T)
if (length(argg) != 3){
  stop('ARGS: 1) Path to input directory (with Leaf/Root subdirs)
       2) Path to desired output directory
       3) Experiment: Sym OR Hydro')
}

##setting arguments
indir <- argg[1]; odir <- argg[2]; exp <- argg[3]

##reading all in as list
filelist <- mk_fil(indir)

##coecring to df
dat_tf <- bind_rows(filelist)
print(head(dat_tf))
dat_tf$nfold_change[dat_tf$nfold_change > 15] <- 15

##adding additional cols
ndat <- dat_tf %>% mutate(logpval = -log10(adj_pval)) %>%
          mutate(threshold = if_else(nfold_change >= 2 & logpval >= 1.30103 |nfold_change <= -2 & logpval >= 1.30103,"A", "B"))%>%
          mutate(nlabel = paste(thistis, thiscond, threshold, sep = ' '))

##make experiment-specific variables
possiblecolors <- c("darkblue", "darkorchid1", "darkolivegreen3", "sandybrown", "darkturquoise",
                    "deeppink1")
totallayers <- length(unique(ndat$nlabel))/2

colors <- c(possiblecolors[1:totallayers], rep('black', totallayers))
print(colors)
labs <- c(unique(ndat$nlabel)[grepl('A$', unique(ndat$nlabel))], rep('', totallayers))

##finally plotting
pdf(paste0(odir, '/DAM_volcano_plot.pdf'))
ggplot(ndat, aes(nfold_change,logpval, colour = interaction(thistis, thiscond, threshold), size = threshold)) +
    geom_point(alpha = 0.5) + 
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5, size = 1) + 
    geom_vline(xintercept = 2, linetype = 2, alpha = 0.5, size = 1) +
    geom_vline(xintercept = -2, linetype = 2, alpha = 0.5, size = 1) +
    scale_colour_manual(values = colors, labels = labs) +
    scale_size_manual(values=c(3,1))+
    xlab("log2 fold change") + ylab("-log10 q-value")+labs(colour = '')+
    theme_bw()+guides(size = F)

dev.off()
