library(ggplot2)
library(ggpubr)
library(reshape2)
library(plyr)

## Script for making plots (line and/or density) for a certain group
## of transcripts or metabolites across our treatments.

get_groups <- function(df){
  ct <- colnames(df[,2:ncol(df)])
  groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  return(groups)
}

get_median_per_group <- function(pkarea, groups){
  outl <- list() ## for metabolic genes only
  for (grp in unique(grps)){
    tmp <- pkarea[,which(grepl(grp, colnames(pkarea)))]
    outl[[grp]] <- matrixStats::rowMedians(as.matrix(tmp))
  }
  grpmed <- do.call(cbind.data.frame, outl)
  grpmed <- cbind.data.frame(pkarea[,1], grpmed)
  colnames(grpmed)[1] <- 'Gene'
  return(grpmed)
}

make_density_plot <- function(goterm, expression, char, odir, typ, onlyheatcop = T){
  if (typ == 'trans'){
    oxgens <- unique(goterm[grepl(char, goterm$V2), 'V1'])
    oxpkarea <- expression[expression$Gene %in% oxgens,]
  }  else {
    oxgens <- unique(goterm[grepl(char, goterm$subclass), 'name'])
    oxpkarea <- expression[expression$Alignment.ID %in% oxgens,]
  }
  
  idv <- ifelse(typ == 'trans', 'Gene', 'Alignment.ID')
  oxdata <- melt(oxpkarea, id.vars = idv)
  if (onlyheatcop){
    oxhydro <- oxdata[which(grepl('Hydro', oxdata$variable)),]
    oxdata <- oxhydro[which(grepl('Root', oxhydro$variable)),]
    cvec <- c('mediumslateblue','goldenrod2',
               'forestgreen',  'darkorchid4')
  } else {
    cvec <- c('salmon', 'mediumslateblue','palegreen', 'goldenrod2',
              'lightblue', 'forestgreen', 'plum3', 'darkorchid4', 'pink',
              'azure4', 'darkturquoise', 'darkred', 'honeydew1', 'mistyrose1',
              'blue2', 'deeppink4', 'gray27')
  }
  oxdata$variable <- as.character(oxdata$variable)
  oxdata$grp <- unlist(lapply(strsplit(oxdata$variable, '_'), '[[', 1))

  oxmu <- ddply(oxdata, "grp", summarise, grp.mean=mean(value))
  pdf(paste0(odir, '/', char, '_expression_density_per_group.pdf'))
  print(ggplot(oxdata, aes(value, color = grp)) + geom_density(size = 1) + 
          theme_bw() + scale_color_manual(values = cvec) +
          geom_vline(data = oxmu, aes(xintercept = grp.mean, color = grp),
                     linetype="dashed", size = 0.75) + ggtitle(paste0(char, ' Expression per Group')) + 
          scale_color_manual(values = cvec))
    
  dev.off()
}

argg <- commandArgs(T)

if (length(argg) != 7){
  stop('ARGS: 1) file that says what gene/metab belongs to what goterm/structural class
       2) Expression file of genes/metabs 3) Do you want a line plot? y OR n
       4) Do you want a density plot? y OR n 5) output directory 6) trans OR metab
       7) What class/goterm do you want to look for?')
}

if (argg[6] == 'trans'){
  goterm <- read.table(argg[1], sep = '\t', fill = NA, quote = "", header = F)
}else goterm <- read.table(argg[1], sep = '\t', fill = NA, quote = "", header = T)

expression <- read.table(argg[2], sep = '\t', fill = NA, quote = "", header = T)

if (argg[6] == 'metab'){
  expression <- expression[,-c(2:32)]
}

## Line plots 
if (argg[3] == 'y'){  
  if (argg[6] == 'trans'){
    oxgens <- unique(goterm[grepl(argg[7], goterm$V2), 'V1'])
    pkarea <- expression[expression$Gene %in% oxgens,]
  } else {
    oxgens <- unique(goterm[grepl(argg[7], goterm$subclass), 'name'])
    pkarea <- expression[expression$Alignment.ID %in% oxgens,]
  }
  
  grps <- get_groups(pkarea)
  
  medianpkarea <- get_median_per_group(pkarea, grps)
  medianpkarea$higher <- ifelse(medianpkarea$Hydro.HeatCop.Root > medianpkarea$Hydro.Heat.Root, 1, 0)
  
  datatoplot <- melt(medianpkarea, id.vars = c('Gene', 'higher'))
  pdf(paste0(argg[5], '/', argg[7], '_line_plots.pdf'))
  print(ggplot(datatoplot, aes(variable, value, group = Gene)) + geom_line(aes(color = higher)) + geom_point(color = 'red') +
          theme(axis.text.x = element_text(angle = 45, size = 9, hjust = 1), legend.position = "none") +
          ggtitle(paste0(argg[7], ' Normalized Expression')) + labs(x = '', y = 'Normalized Expression'))
  dev.off()
}


## Now making a density plot
if (argg[4] == 'y'){
  make_density_plot(goterm, expression, argg[7], argg[5], argg[6], onlyheatcop = T)
}




