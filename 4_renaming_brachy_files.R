##Script for renaming brachy files

change_names <- function(path, cut){
  ##path is a string, full path to files to rename
  ##cut is a str, the characters to split the name on
  torenam <- list.files(path, pattern = '.mzML', full.names = T)
  out <- c()
  for (nam in torenam){
    bnam <- basename(nam)
    short <- strsplit(bnam, cut, fixed = T)[[1]][2]
    shorter <- gsub('Rg100to800-CE102040-Brachy-S1_', '', short, fixed=T)
    shorter_full <- paste0(path, shorter)
    out <- c(out, shorter_full)
  }
  return(out)
}

get_run_order <- function(newn){
  ##newn is a character vector of new names
  runs <- unlist(lapply(strsplit(newn, '_Run', fixed = T), '[[', 2))
  rest_of_name <- unlist(lapply(strsplit(newn, '_Run', fixed = T), '[[', 1))
  runs <- gsub('.mzML', '', runs, fixed = T)
  nruns <- as.numeric(runs)
  runs_order <- paste(runs[order(runs)], 1:length(runs), sep = '_')
  runs_order <- 
  out <- paste(rest_of_name[order(runs)], '_Run', runs_order, '.mzML', sep = '')
  return(out)
}

bigdir <- 'C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/'
dirs <- c('Copper_Heat/', 'Symbiosis/', 'Tissue_ctrl/')
cuts <- c('Hydro-', 'Sym-', 'Tis-')
types <- c('POS/', 'NEG/')

for (dir in 1:length(dirs)){
  for (typ in types){
    fildir <- paste0(bigdir, dirs[dir], typ)
    oldnames <- list.files(fildir, pattern = '.mzML', full.names = T)
    ##have to order oldnames
    smallnam <- unlist(lapply(strsplit(oldnames, '_Run', fixed = T), '[[', 2))
    nameruns <- gsub('.mzML', '', smallnam, fixed = T)
    order_oldnames <- oldnames[order(nameruns)]
    ##getting newnames
    newnames <- change_names(fildir, cuts[dir])
    order_newnames <- get_run_order(newnames)
    for (fil in 1:length(order_oldnames)){
      print(order_oldnames[fil])
      print(order_newnames[fil])
      file.rename(order_oldnames[fil], order_newnames[fil])
    }
  }
}
