## take in clusterone clusters.csv file, genes of interest, and canopus annotation file
## parse IDs of canopus file
## get the clusters that the genes of interest are in
## for each gene - what metabolites are it clustered with?
## What's their most specifc class?
library(readxl)

argg <- commandArgs(T)

if (length(argg) != 7){
  stop('ARGS: 1) Brachy genes of interest file
       2) Positive spearman clusterONE output
       3) Positive Pcc clusterone output
       4) Negative spearman custerONE output
       5) Negative PCC clusterONE output
       6) Canopus output positive
       7) Canopus output negative')
}

genesint <- read_excel(argg[1])
pos_spr <- read.csv(argg[2])
pos_pcc <- read.csv(argg[3])
neg_spr <- read.csv(argg[4])
neg_pcc <- read.csv(argg[5])
canopusp <- read.table(argg[6], header = T, sep = '\t', fill = NA, quote = "")
canopusn <- read.table(argg[7], header = T, sep = '\t', fill = NA, quote = "")


parse_id <- function(canopus){
  if (grepl('scan', canopus$name[1])){
    tl <- unlist(lapply(strsplit(canopus$name, '_'), '[[', 4))
    canopus$name <- substr(tl, 6, nchar(tl))
  }
  return(canopus)
}

get_mets <- function(genesint, clusters){
  
  geneclusters <- c()
  for (i in 1:nrow(genesint)){
    tmp <- paste(clusters$Cluster[which(grepl(genesint$AMS_gene[i], clusters$Members))], collapse = ',')
    geneclusters <- c(geneclusters, tmp)
  }
  
  linkedmets <- c()
  
  for (clu in geneclusters){
    mets <- c()
    if (nchar(clu) > 0){
      clustersalone <- unlist(strsplit(clu, ','))
      for (tclu in clustersalone) {
        row <- clusters[which(clusters$Cluster == tclu),]
        member <- gsub('"', '', row$Members)
        members <- unlist(strsplit(member, ' '))
        if (any(grepl("^[[:digit:]]+", members))){
          thesemets <- members[which(grepl('^[[:digit:]]+', members))]
          mets <- c(mets, thesemets)
        }
      }
    } 
    linkedmets <- c(linkedmets, paste(unique(mets), collapse = '|'))
  }
  return(linkedmets)
}

get_mets_union <- function(spr, pcc){
  out <- c()
  for (metset in 1:length(spr)){
    sprmet <- unlist(strsplit(spr[metset], '|', fixed = T))
    pccmet <- unlist(strsplit(pcc[metset], '|', fixed = T))
    union <- unique(c(sprmet, pccmet))
    if (length(union) > 0){
      outstr <- paste(union, collapse = '|')
    } else outstr <- ''
    out <- c(out, outstr)
  }
  return(out)
}

get_annots <- function(linkedmet, canopus){
  annots <- c()
  for (mets in linkedmet){
    if (nchar(mets) > 0){
      annot <- c()
      mets <- unlist(strsplit(mets, '|', fixed = T))
      for (met in mets){
        if (met %in% canopus$name){
          annot <- c(annot, canopus$most.specific.class[canopus$name == met])
        } else annot <- c(annot, 'None')
      }
      annot <- paste(annot, collapse = '|')
    } else {
      annot <- ''
    }
    annots <- c(annots, annot)
  }
  return(annots)
}

canopusp <- parse_id(canopusp)
canopusn <- parse_id(canopusn)

pospccmets <- get_mets(genesint, pos_pcc)
possprmets <- get_mets(genesint, pos_spr)
negpccmets <- get_mets(genesint, neg_pcc)
negsprmets <- get_mets(genesint, neg_spr)

genesint$Pos_Mets <- get_mets_union(possprmets, pospccmets)
genesint$Neg_Mets <- get_mets_union(negsprmets, negpccmets)
genesint$Pos_CANOPUS <- get_annots(genesint$Pos_Mets, canopusp)
genesint$Neg_CANOPUS <- get_annots(genesint$Neg_Mets, canopusn)

write.table(genesint, file = paste0(argg[1], '_linkedmetabs.tab'), sep = '\t', row.names = F, quote = F)
