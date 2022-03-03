## Script to make boxplots of peak areas. 

library(ggplot2)
library(ggpubr)
library(reshape2)

argg <- commandArgs(T)

if (length(argg) != 4){
  stop('ARGS: 1) Peak Area totally raw 2) Peak Area after KNN 
       3) Peak Area after VSN, with no negative values
       4) Directory and name of output file')
}

get_good_columns <- function(df, good_colnames){
  return(df[, which(colnames(df) %in% good_colnames)])
}

get_groups <- function(vec){
  return(unlist(lapply(strsplit(vec, '_'), '[[', 1)))
}

make_good_names <- function(names){
  names <- unlist(lapply(strsplit(names, 'Run'), '[[', 1))
  names <- substr(names, 1, nchar(names)-1)
  return(names)
}

raw <- read.table(argg[1], header = T, sep ='\t', quote = "", fill = NA)
knn <- read.table(argg[2], header = T, sep ='\t', quote = "", fill = NA)
vsn <- read.table(argg[3], header = T, sep ='\t', quote = "", fill = NA)

vsn_use <- vsn[, -c(1:32)]
raw_use <- get_good_columns(raw, colnames(vsn_use))
knn_use <- get_good_columns(knn, colnames(vsn_use))

mvsn <- melt(vsn_use)
mraw <- melt(raw_use)
mknn <- melt(knn_use)

mvsn$step <- 'VSN'
mraw$step <- 'Raw'
mknn$step <- 'KNN'

mtoplot <- rbind(mraw, mknn, mvsn)
mtoplot$variable <- as.character(mtoplot$variable)
mtoplot$groups <- get_groups(mtoplot$variable)
mtoplot$variable <- make_good_names(mtoplot$variable)

colpal <- c('black', 'brown2', 'chartreuse3', 'dodgerblue2', 'cyan', 'magenta',
            'goldenrod', 'grey50', 'darkgreen', 'darkblue', 'yellow1', 'violetred4',
            'yellowgreen', 'springgreen', 'pink', 'wheat4', 'lightcoral')

## Plotting now, no x-axis tick labels. Assuming there's a figure already with them

pdf(paste0(argg[4], '_normalization_steps_boxplots.pdf'), height = 4, width = 10)
ggplot(mtoplot, aes(variable, value, fill = groups)) + geom_boxplot() + theme_bw() + 
  facet_wrap(~step, scales = 'free', nrow = 3, ncol = 1) + scale_fill_manual(values = colpal) + 
  theme(axis.text.x = element_blank())
 dev.off()
 
 