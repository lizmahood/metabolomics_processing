library(xcms)
library(RColorBrewer)
library(SummarizedExperiment)
library(pander)
library(magrittr)

##This script is for getting the TICs of all files in the three experiments

get_tic <- function(fils, datt, outp){
  ##fils is character vector of paths to files to read in
  ##datt is a phenodata df
  ##outp is string, out path of plots
  raw <- readMSData(files = fils, pdata = new('NAnnotatedDataFrame', datt), mode = 'onDisk')
  print('Done reading In')
  tics <- chromatogram(raw, aggregationFun = "sum")
  typs <- length(unique(datt$conds))
  group_c <- paste0(brewer.pal(typs, 'Set3'), '60')
  names(group_c) <- unique(datt$conds)
  ##plotting TIC
  pdf(paste0(outp, '/TIC.pdf'))
  plot(tics, col = group_c[datt$conds])
  legend("topright", legend = names(group_c), fill = unname(group_c), ncol = 2)
  dev.off()
  ##plotting zoomed in
  pdf(paste0(outp, '/zoomed_in_TIC.pdf'))
  plot(tics, col = group_c[datt$conds], ylim = c(0, 2.5e+08))
  legend("topright", legend = names(group_c), fill = unname(group_c), ncol = 2)
  dev.off()
}

make_sample_nam <- function(nam, pat){
  ##nam is character vector
  ##pat is pattern to find name after
  out <- c()
  for (i in 1:length(nam)){
    bignam <- strsplit(nam[i], pat, fixed = T)[[1]][2]
    sampn <- paste(strsplit(bignam, '_', fixed = T)[[1]][c(1,2)], collapse ='_')
    out <- c(out, sampn)
  }
  return(out)
}

make_condition <- function(snams, typ){
  ##snams is a vector of sample names
  ##typ is string -- are we doing sym, copper, or ctrl?
  conds <- c()
  if (typ == 'Sym-'){
    for (i in 1:length(snams)){
      if (grepl('Leaf', snams[i], fixed = T)) tmp_t <- 'L'
      else if (grepl('Root',snams[i], fixed = T)) tmp_t <- 'R'
      if (grepl('W',snams[i], fixed = T)) tmp_c <- 'SW'
      else if (grepl('Spore',snams[i], fixed = T)) tmp_c <- 'S'
      if (grepl('Ex', snams[i], fixed = T)) tmp_c <- 'BLK'
      else if (grepl('Ctrl', snams[i], fixed = T)) tmp_c <- 'CTRL'
      out <- paste(tmp_t, tmp_c, sep = '_')
      conds <- c(conds, out)
    }
  }
  if (typ == 'Tis-'){
    for (i in 1:length(snams)){
      print(i)
      if (grepl('Culm', snams[i], fixed = T)) tmpt <- 'Culm'
      else if (grepl('Leaf', snams[i], fixed = T)) tmpt <- 'Leaf'
      else if (grepl('Spike', snams[i], fixed = T)) tmpt <- 'Spklt'
      else if (grepl('ExCt', snams[i], fixed = T)) tmpt <- 'BLK'
      conds <- c(conds, tmpt)
    }
  }
  if (typ == 'Hydro-'){
    for (i in 1:length(snams)){
      if (grepl('HeatRoot', snams[i], fixed = T)) tmpt <- 'BLK_HR'
      else if (grepl('ExCtrl', snams[i], fixed = T)) tmpt <- 'BLK_H'
      else if (grepl('HeatCop', snams[i], fixed = T)){
        if (grepl('Leaf', snams[i], fixed = T)) tmpt <- 'HCP_L'
        else tmpt <- 'HCP_R'
      }
      else if (grepl('Ctrl', snams[i], fixed = T)){
        if (grepl('Leaf', snams[i], fixed = T)) tmpt <- 'CTRL_L'
        else tmpt <- 'CTRL_R'
      }
      else if (grepl('Cop', snams[i], fixed = T)){
        if (grepl('Leaf', snams[i], fixed = T)) tmpt <- 'CP_L'
        else tmpt <- 'CP_R'
      }
      else if (grepl('Heat', snams[i], fixed = T)){
        if (grepl('Leaf', snams[i], fixed = T)) tmpt <- 'H_L'
        else tmpt <- 'H_R'
      }
      conds <- c(conds, tmpt)
    }
  }
  return(conds)
}

make_pd <- function(snams, conds){
  return(data.frame(snams, conds, stringsAsFactors = F))
}

bigpath <- 'C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files'
conditions <- c('Copper_Heat', 'Symbiosis', 'Tissue_ctrl')
pats <- c('Hydro-', 'Sym-', 'Tis-')
modes <- c('FPS', 'NEG', 'POS')

for (cnds in 1:length(conditions)){
  for (mod in modes){
    filpth <- paste(bigpath, conditions[cnds], mod, sep = '/')
    fils <- list.files(filpth, full.names = T, pattern = 'mzML')
    snames <- make_sample_nam(fils, pats[cnds])
    condts <- make_condition(snames, pats[cnds])
    pddf <- make_pd(snames, condts)
    get_tic(fils, pddf, filpth)
    print(paste0('Done with ', mod))
  }
}

