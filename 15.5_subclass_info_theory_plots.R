## For an additional information theory plot
## That plots subclasses on average pk area change per stress plots

library(reshape2)
library(ggplot2)
library(tidyverse)

argg <- commandArgs(T)

if (length(argg) != 4){
  stop('ARGS: 1) canopus file 2) pk area change per metabolite file
       3) which class/stress to plot subclasses of? Enter as:
       Class1-Stress.Name.Tissue1,Class2-Stress.Name.Tissue2
       Put spaces as _
       4) output directory and prefix')
}

parse_classes <- function(class_str){
  
  class_str <- gsub('_', ' ', class_str)
  
  if (grepl(',', class_str)){
    class_vec <- unlist(strsplit(class_str, ','))
  }else class_vec <- class_str
  
  class_df <- as.data.frame(class_vec)
  
  outdf <- class_df %>% separate(class_vec, c('Class', 'Stress'), sep = '([-])')
  return(outdf)
}

plot_subclasses <- function(canopus_sub, pc, pkareachange, ofil){
  pdf(paste0(ofil, '_subclass_avg_pk_area_change_plots.pdf'))
  for (rw in 1:nrow(pc)){
    pkareachange_tmp <- pkareachange[pkareachange$Stress == pc[rw, 2] &
                                       pkareachange$Class == pc[rw, 1],]
    toplot <- merge(pkareachange_tmp, canopus_sub, 
                    by.x = 'Alignment.ID', by.y = 'name')
    toplot[toplot$subclass == '', 'subclass'] <- 'None'
    
    print(ggplot(toplot, aes(Stress, AvgPkAreaChange, color = subclass)) + 
      geom_jitter(size = 2) + theme_minimal() + 
      scale_color_brewer(palette = 'Paired')) #theme(plot.margin = unit(c(1,1,1,1),"cm"))
    
  }
  dev.off()
}

canopus <- read.table(argg[1], sep = '\t', fill = NA, quote = "", header = T)
canopus_sub <- canopus[, c('name', 'subclass')]
pkareachange <- read.table(argg[2], sep = '\t', fill = NA, quote = "", header = T)
pc <- parse_classes(argg[3])

plot_subclasses(canopus_sub, pc, pkareachange, argg[4])
