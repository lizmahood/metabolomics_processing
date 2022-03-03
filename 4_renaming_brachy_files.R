## Script for renaming brachy files

change_names <- function(path, cut, all){
  ## path is a string, full path to files to rename
  ## cut is a str, the characters to split the name on
  ## all is another string -- all files together?
  torenam <- list.files(path, pattern = '.mzML', full.names = T)
  out <- c()
  for (nam in torenam){
    bnam <- basename(nam)
    short <- strsplit(bnam, cut, fixed = T)[[1]][2]
    if (all == 'yes'){
      short <- gsub('^\\d*_', '', short)
    }
    if (grepl('QC-Moghe', short)){
      shorter <- gsub('Rg100to800-CE102040-Brachy-QC-Moghe_', '', short)
    }else {
      shorter <- gsub('Rg100to800-CE102040-Brachy-S1_', '', short, fixed=T)
    }
    shorter_full <- paste0(path, shorter)
    out <- c(out, shorter_full)
  }
  return(out)
}

get_run_order <- function(newn){
  ## newn is a character vector of new names
  runs <- unlist(lapply(strsplit(newn, '_Run', fixed = T), '[[', 2))
  rest_of_name <- unlist(lapply(strsplit(newn, '_Run', fixed = T), '[[', 1))
  for (i in 1:length(runs)){
    run <- runs[i]
    if (grepl('_', run, fixed = T)){
      runs[i] <- strsplit(run, '_', fixed = T)[[1]][1]
    }
  }
  runs <- gsub('.mzML', '', runs, fixed = T)
  nruns <- as.numeric(runs)
  runs_order <- paste(runs[order(nruns)], 1:length(runs), sep = '_')
  out <- paste(rest_of_name[order(nruns)], '_Run', runs_order, '.mzML', sep = '')
  return(out)
}

argg <- commandArgs(T)

if (length(argg) != 2){
  stop('ARGS: 1) Path to directory of files to rename 
       2) Is this all files together in the Brachy expt? yes OR no')
}

cuts <- 'MSMS_'
fildir <- argg[1]
oldnames <- list.files(fildir, pattern = '.mzML', full.names = T)
## have to order oldnames
smallnam <- unlist(lapply(strsplit(oldnames, '_Run', fixed = T), '[[', 2))
nameruns <- gsub('.mzML', '', smallnam, fixed = T)
order_oldnames <- oldnames[order(as.numeric(nameruns))]
print(head(order_oldnames))
## getting newnames
newnames <- change_names(fildir, cuts, argg[2])
order_newnames <- get_run_order(newnames)

## renaming
for (fil in 1:length(order_oldnames)){
  print(order_oldnames[fil])
  print(order_newnames[fil])
  file.rename(order_oldnames[fil], order_newnames[fil])
}


