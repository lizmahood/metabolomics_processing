##script for mapping STEM clusters to MetaNet network output files
##ARGS: 1) MetaNet correlation matrix, 2) STEM metabolite cluster output
##3) nodes to consider (significant nodes) separated by commas
##EX: 2,4,5,7,10 
##4) STEM input file

######Functions########

get_all_metabs <- function(sin, sout){
  ##sin is df, stem input
  ##sout is df, stem output
  ##RETURNS: df, merged sin and sout

  ##split up combined SPOTS
  bad <- which(is.na(as.numeric(sout$SPOT)))
  for (i in bad){
    sp <- sout[i, 'SPOT']
    smaldf <- sout[i,]
    ids <- unlist(strsplit(sp, ';', fixed = T))
    for (j in ids){
      smaldf$SPOT <- j
      sout <- rbind(sout, smaldf)
    }
  }
  ##merging
  if (length(which(grepl(';', sout$SPOT, fixed = T))) > 0){
    sout <- sout[-which(grepl(';', sout$SPOT, fixed = T)),]
  }
  out <- merge(sin, sout, by = 'SPOT', all = T)
  return(out)
}

get_sig_clusters <- function(clusts, sigc){
  ##clusts is a numeric vector
  ##sigc is numeric vector
  ##returns numeric vector with non-sig clusters
  ##turned to 1
  clusts[is.na(clusts)] <- 1
  clusts[which(!(clusts %in% sigc))] <- 1
  return(clusts)
}

make_edge_table <- function(corf, ids){
  ##corf is the correlation network output of 
  ##MetaNetter
  ##ids is a numeric vector, alignment ids
  ##Returns: df of highly correlated pairs of ids
  ##(This is the edge file)
  rownames(corf) <- paste0('X', ids)
  corf[,1] <- NULL
  goodfirst <- as.numeric(corf[1,])
  corf [,1] <- goodfirst
  colnames(corf) <- paste0('X', ids)
  corf[1,1] <- 0
  ind <- which(upper.tri(corf, diag = TRUE), arr.ind = TRUE)
  nn <- dimnames(corf)
  test <- data.frame(row = nn[[1]][ind[, 1]],
                     col = nn[[2]][ind[, 2]],
                     val = corf[ind])
  out <- test[which(test$val == 1),]
  out$row <- gsub('X', '', out$row)
  out$col <- gsub('X', '', out$col)
  return(out)
}

#######Body########

argg <- commandArgs(T)
if (length(argg) != 4){
  stop('ARGS: 1) MetaNet correlation matrix, 2) STEM metabolite cluster output
  3) nodes to consider (significant nodes) separated by comma,
       EX: 2,4,5,7,10   4) STEM input file')
}

cornet <- read.table(argg[1], fill = NA, header = T, quote = "")
stemout <- read.table(argg[2], header = T, fill = NA, quote = "", sep = '\t')
stemin <- read.table(argg[4], header = T, quote = "")
sigclus <- as.numeric(unlist(strsplit(argg[3], ',', fixed = F)))

allmets <- get_all_metabs(stemin, stemout)
goodclusts <- get_sig_clusters(allmets$Profile, sigclus)
allmets$Profile <- goodclusts
edge <- make_edge_table(cornet, allmets$SPOT)
node <- data.frame(allmets$SPOT, allmets$Profile)

##writing out 
write.table(edge, file = paste0(argg[1], '_edge.txt'), sep = '\t', row.names = F, quote = F)
write.table(node, file = paste0(argg[1], '_node.txt'), sep = '\t', row.names = F, quote = F)
