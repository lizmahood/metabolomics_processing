## Quantifying correspondance between GNPS (MSDIAL) and CANOPUS

argg <- commandArgs(T)

annots <- read.table(argg[1], header = T, sep = '\t', quote = "", fill = NA)

## Removing rows if Canopus did not annotate them at all
annots <- annots[-which(annots$superclass == ''),]

## What % of superclass, class, subclass, and level 5 do they agree on?
aggrement <- c()

classis <- annots[, c('superclass', 'class', 'subclass', 'level.5',
                      'CF_Superclass', 'CF_Class', 'CF_Subclass', 'CF_Parent.Level.1')]

colnames(classis)[5:8] <- c('CF_superclass', 'CF_class', 'CF_subclass', 'CF_level.5')

for (level in 1:4){
  this_canopus <- colnames(classis)[level]
  this_cf <- paste0('CF_', this_canopus)
  
  not_classified <- which(classis[,level] == 'None')
  both <- classis[,c(this_canopus, this_cf)]
  both <- both[-not_classified,]
  
  print(this_canopus)
  print(nrow(both))
  print((length(which(both[,1] == both[,2]))))
  print((length(which(both[,1] == both[,2]))) / nrow(both))
}

