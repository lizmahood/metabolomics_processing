##script for graphing. Makes volcano plots from tissue specific metabolites
library(dplyr)
library(ggplot2)
library(tidyverse)
library(grid)

mk_fil <- function(indir, string){
  ##indir is a string, path to directory with Leaf and Root subdirs
  ##string is a string, pattern in files to read in
  ##RETURNS: list with each tissue-cond pair as a dataframe
  tissues <- c('Leaf', 'Root')
  fillist <- list()
  for (tis in tissues){
    readdir <- paste(indir, tis, 'MSDIAL_out', sep = '/')
    print(readdir)
    fils <- list.files(readdir, pattern = string, full.names = T)
    if (length(fils) <1){
      print('wow')
      readdir <- paste0(indir, '/MSDIAL_out')
      fils <- list.files(readdir, pattern = string, full.names = T)
    }
    for (fil in fils){
      temp <- read.table(fil, sep = '\t', header = T, quote = "", fill = NA)
      thistis <- rep(tis, nrow(temp))
      ##python has spoiled me
      tcond <- strsplit(fil, '_', fixed = T)[[1]][length(strsplit(fil, '_', fixed = T)[[1]])]
      print(tcond)
      cond <- strsplit(tcond, '.', fixed = T)[[1]][1]
      thiscond <- rep(cond, nrow(temp))
      temp <- cbind(temp, thistis, thiscond)
      ##making additional columns
      ndat <- temp %>% mutate(logpval = -log10(adj_pval)) %>%
        mutate(threshold = if_else(nfold_change >= 2 & logpval >= 1.30103 |nfold_change <= -2 & logpval >= 1.30103,"A", "B"))
      ndat$nfold_change[ndat$nfold_change >= 10] <- 10
      fillist <- c(fillist, list(ndat))
    }
  }
  return(fillist)
}

make_volcano <- function(fillist, odir){
  ##fillist is a list, odir is a string
  ##OUTPUTS: one volcano plot per tissue/condition pair
  for (fil in fillist){
    ndat <- data.frame(fil)
    tis <- ndat$thistis[1]; cond <- ndat$thiscond[1]
    grob <- grobTree(textGrob(paste0("n = ", length(which(ndat$nfold_change >= 2 & ndat$logpval >= 1.30103))), x=0.8,  y= 0.85, hjust=0,
                              gp=gpar(col="blue", fontsize=15)))
    grob2 <- grobTree(textGrob(paste0("n = ", length(which(ndat$nfold_change <= -2 & ndat$logpval >= 1.30103))), x=0.1,  y= 0.84, hjust=0,
                               gp=gpar(col="blue", fontsize=15)))
    
    if (tis == 'Leaf') color <- 'darkgreen'
    if (tis == 'Root') color <- 'saddlebrown'
    pdf(paste0(odir, '/', cond, tis, 'DAM_volcano_plot.pdf'))
    plot <- ggplot(ndat, aes(nfold_change,logpval, colour = interaction(thistis, thiscond, threshold), size = threshold)) +
      geom_point(alpha = 0.5) + 
      geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5, size = 1) + 
      geom_vline(xintercept = 2, linetype = 2, alpha = 0.5, size = 1) +
      geom_vline(xintercept = -2, linetype = 2, alpha = 0.5, size = 1) +
      scale_colour_manual(values = c(color, 'black'), labels = cond) +
      scale_size_manual(values=c(3,1))+
      xlab("log2 fold change") + ylab("-log10 q-value")+labs(colour = '')+
      annotation_custom(grob)+ annotation_custom(grob2)+
      theme_bw()+guides(size = F, colour = F) + ggtitle(paste(tis, cond, sep = ' '))
    print(plot)
  }
}

make_barplot <- function(fillist, odir){
  ##fillist is a list, output of mk_fil
  ##odir is a string
  ##OUTPUTS: barplot
  odf <- data.frame('condition' = NA, 'count' = NA, 'direction' = NA)
  for (fil in fillist){
    ndat <- data.frame(fil)
    tis <- ndat$thistis[1]; cond <- ndat$thiscond[1]
    up <- length(which(ndat$nfold_change >= 2 & ndat$logpval >= 1.30103))
    down <- length(which(ndat$nfold_change <= -2 & ndat$logpval >= 1.30103))
    odf <- rbind(odf, c(paste(tis, cond, sep = ' '), as.numeric(up), 'Up'))
    odf <- rbind(odf, c(paste(tis, cond, sep = ' '), as.numeric(down), 'Down'))
  }
  odf <- odf[complete.cases(odf),]
  odf$count <- as.numeric(odf$count); odf$condition <- as.factor(odf$condition)
  odf$direction <- as.factor(odf$direction)
  
  pdf(paste0(odir, '/barchart_diff_abund_metabs.pdf'))
  print(ggplot(odf, aes(x =condition,y=count)) + geom_bar(stat = 'identity') + 
    facet_grid(direction ~ ., scales = "free_y") +
    theme(strip.text = element_text(size = 20, colour = 'orange'), 
          axis.text.x = element_text(size = 10, color = 'orange'),
          axis.text.y = element_text(size = 10, color = 'orange'))+
    theme_bw()+ ylab('') + xlab(''))
}

argg <- commandArgs(T)
if (length(argg) != 2){
  stop('ARGS: 1) Path to input directory (with Leaf/Root subdirs)
       2) Path to desired output directory')
}

##setting arguments
indir <- argg[1]; odir <- argg[2]

##reading all in as list
allmetabs <- mk_fil(indir, 'all_pval_fc')

make_volcano(allmetabs, odir)
print('Done with volcano plots!')
make_barplot(allmetabs, odir)
print('All done!')