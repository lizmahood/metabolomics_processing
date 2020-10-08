##script for plotting shannon entropy
library(ggplot2)
library(ggrepel)

argg <- commandArgs(T)
if (length(argg) != 2){
  stop('ARGS: 1) Path to shannon entropy tab file 2) output directory')
}

shannon <- read.table(argg[1], sep = '\t', stringsAsFactors = F, quote = "", header = T)
shannon$Condition <- gsub('_', ' ', shannon$Condition, fixed = T)
shannon$Condition <- gsub('.', ' ', shannon$Condition, fixed = T)

pdf(paste0(argg[2], '_shannon_entropy_plot.pdf'))
ggplot(shannon, aes(Total_peak_num, shannon)) + geom_point(size = 2, color = 'navy') +
  theme_bw() + xlab('Total Number of Peaks') + ylab('Shannon Entropy') + 
  geom_text_repel(aes(label = Condition),size = 4.5, color = 'saddlebrown')

dev.off()