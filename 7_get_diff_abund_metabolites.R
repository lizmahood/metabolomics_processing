###Prettier script for analyzing metabolomic data. Finds differentially
###expressed metabolites for treatments. Normalizes metabolites
###first, then finds differentially expressed ones between healthy and
###not. Does FDR correction.

library(matrixStats)
library(vsn)
library(MetaboDiff)
library(dplyr)


                  ###################      
                  #### Functions ####
                  ###################

get_groups <- function(df){
  ct <- colnames(df)
  groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  return(groups)
}

checkmissing <- function(asycl, cutoff){
  ##Should not include blank
  bad <- c()
  for (i in 1:ncol(asycl)){
    cl <- asycl[,i]
    cl[cl == 0] <- NA
    print (paste0(colnames(asycl)[i], ' ', (length(which(is.na(cl))) / length(cl))))
    if ((length(which(is.na(cl))) / length(cl)) >= cutoff) bad <- c(bad, i)
  }

  if (length(bad) > 0){
    stop(paste0('ERROR! Samples ', paste0(colnames(asycl)[bad], collapse = ', '), ' had more missing values than you allow!!'))
  }
}

checksn <- function(asycl, snm, snthrsam){
  ##Should not include blank
  snasy <- as.matrix(snm[,asycl])
  bad <- c()
  for (samp in 1:ncol(snasy)){
    print(paste0(colnames(snm)[asycl[samp]], ' ', mean(snasy[,samp])))
    if (mean(snasy[,samp]) < snthrsam) bad <- c(bad, samp)
  }
  
  if (length(bad) > 0){
    stop(paste0('ERROR! Samples ', paste0(colnames(snasy)[bad], collapse = ', '), ' had lower S/N than you allow!!'))
  }
}

make_blank_single <- function(asy){
  print(colnames(asy))
  inblank <- c()
  bcol <- ncol(asy)
  
  for (i in 1:nrow(asy)){
    if (asy[i, bcol] !=0 & (asy[i, bcol] / min(asy[i,(1:bcol-1)])) > 0.5){
      inblank <- c(inblank, i)
    } 
  }
  return(inblank)
}

make_blank_mult <- function(pkl, blkc, blkg){
  ##pkl is whole alignment table, not just asy
  ##blkc is numeric vector of blank column indices
  ##blkg is a string -- which blanks go with which columns
  
  ##parsing the blkg string into a list
  blkgroups <- unlist(strsplit(blkg, '[', fixed = T))
  blkgroups <- blkgroups[blkgroups != '']
  inblank <- c()
  
  for (grp in blkgroups){
    ##parsing the string
    grp <- gsub(']', '', grp, fixed = T)
    grplst <- strsplit(grp, '-', fixed = T)[[1]]
    blkcol <- as.numeric(grplst[1])
    asycols <- as.numeric(unlist(strsplit(grplst[2], ',')))
    for (i in 1:nrow(pkl)){
      if (pkl[i, blkcol] != 0 & (pkl[i, blkcol] / min(pkl[i, asycols])) > 0.5){
        inblank <- c(inblank, i)
      }
    }  
  }
  return(unique(inblank))
}

imp_per_group <- function(df, groups, mdat){
  ##Output dataframe will have sample columns
  ##at the end of the dataframe. Blank should
  ##NOT be included at this step
  
  tdf <- as.data.frame(t(df))
  tdf$groups <- groups
  
  ##of non-zero/non-noise rows per group
  gil <- list()
  
  ##list of imputed values
  implist <- list()
  
  ##imputing for each group separately
  for (grp in unique(tdf$groups)){
    timp <- tdf[which(tdf$groups == grp),]
    imp <- t(timp)
    imp <- as.matrix(imp[-nrow(imp),])
    imp[imp == 0] <- NA
    noise <- which(rowSums(is.na(imp)) >= (0.5 * ncol(imp)))
    good <- which(rowSums(is.na(imp)) < (0.5 * ncol(imp)))
    goodimp <- imp[-noise,]
    knnimp <- impute::impute.knn(goodimp)$data
    gil[[grp]] <- good
    #names(gil)[length(names(gil))] <- grp
    implist[[grp]] <- knnimp
  }
  
  ##making copy of df to remove rows from
  outdf <- df
  
  ##removing rows if they are zero/noise across ALL conditions
  for (rw in 1:nrow(df)){
    inn <- F
    conds <- c()
    for (vec in 1:length(gil)){
      if (rw %in% gil[[vec]]) {
        inn <- T
        conds <- c(conds, vec)
      }
    } 
    
    if (inn == F) {
      outdf <- outdf[-rw,]
      mdat <- mdat[-rw,]
    }
    ##if they are not zero/noise across ALL conditions,
    ##combine the zero/noise conditions with imputed conditions
    else {
      ##for each group it was good in
      for (cnd in conds){
        ##get the columns of that group  
        cls <- which(tdf$group == names(gil)[cnd])
        ##get the imputed values for those groups
        impvals <- unname(implist[[cnd]][unname(which(gil[[cnd]] == rw)),])
        ##replace original values with these values
        for(cl in 1:length(cls)){
          outdf[rw, cls[cl]] <- impvals[cl]
        }
      }
    }
  }
  return(list(mdat, outdf))
}

vs_norm <- function(pl, idcol, ccols, tcols, blkcol, sweights, sncol, 
                    snthresh, cutoff, imp, telid, vsnn, isnm, snmatt, snsampt, blkgrp){
  
  #'pl: peak list file
  #'c/tcols: columns of pl containing samples values of metabolites in ctrl/treat conditions
  #'blkcol: column index of the blank
  #'sweights: vector of sample weights in the same order as samples.
  #'cutoff: for imputing a line with missing values more than this
  #'
  #'First: S/N threshold. 2) noise threshold. 3: Blank filter. 4: sample weight 
  #'normalization. 5: Detect and remove replicate outliers via correlation 
  #'6: impute (if wanted) 7: vsn (if wanted)
  
  if (tcols != 0){
    asycols <- c(ccols, tcols)
  }else asycols <- ccols
  #print(colnames(pl)[asycols])
  print(paste0('Total metabolites: ', nrow(pl)))
  print(paste0('Total fragmented metabolites: ', length(which(pl$MS.MS.spectrum != ''))))
  
  ##which rows have fragmented metabolites
  frags <- which(pl$MS.MS.spectrum != '')
  
  ### CHECKS for %0 values and S/N values
  checkmissing(pl[,asycols], cutoff)
  checksn(asycols, snmatt, snsampt)
  
  ### 1) doing sn thresholding
  pl[,sncol] <- as.numeric(pl[,sncol])
  plr <- nrow(pl)
  ##assuming internal std has survived SN threshold
  pl <- pl[which(pl[,sncol] > snthresh),]
  snfrags <- which(pl$MS.MS.spectrum != '')
  
  ##getting internal standard row
  print(length(which(pl[,idcol] == telid)))
  telrow <- which(pl[,idcol] == telid)
  if (length(telrow) > 1) telrow <- telrow[1]
  #print(telrow)
  #print(pl[telrow,])
  
  ##getting sample values
  asy <- pl[,c(asycols, blkcol)]
  #asy <- pl[,asycols]
  asy <- apply(asy, 2, as.numeric)
  #print('Correlation of raw values after s/n threshold')
  #print(cor(asy))
  asyr <-  nrow(asy)
  print(paste0('Num of metabolites removed from SN: ', (plr - asyr)))
  print(paste0("Num of fragmented metabolies removed: ", (length(frags) - length(snfrags))))
  badrow <- c()
  
  ### 2 Apply noise threshold: any peak with max value lower than 10k removed
  lowrow <- which(rowMaxs(asy) < 10000)
  badrow <- c(lowrow)
  ##writing out correlation before blank
  cort <- cor(asy[-lowrow,])
  write.table(cort, file = paste0(args[1], '_cor_before_blank.tab'), sep = '\t', quote = F)
  print(paste0('Num of metabolites removed from low values: ', length(lowrow)))
  print(paste0('Num of fragmented metabolites removed: ', (length(which(snfrags %in% badrow)))))
  
  ### 3) Blank filter: if average metabolite value in blank is >= 50% of
  ### min value in sample -- remove
  if (length(blkcol) == 1){
    inblank <- make_blank_single(asy)
  }else if (length(blkcol) > 1){
    inblank <- make_blank_mult(pl, blkcol, blkgrp)
  }
  
  print(paste0('Num of metabolites in blank: ', length(inblank)))
  print(paste0('Num of fragmented metabolites removed: ', length(which(snfrags %in% inblank))))
  #print(inblank[which(inblank == telrow)])
  if (telrow %in% inblank) inblank <- inblank[-which(inblank == telrow)] ##need to keep internal stndard
  badrow <- c(badrow, inblank)
  badrow <- unique(badrow)
  #print('Correlation after blank filter')
  #print(cor(asy[-badrow,]))
  print(paste0('Num of metabolites removed from blank or low row: ', length(badrow)))
  print(colnames(asy))

  ### 4) Normalizing to sample weight
  for (i in 1:length(sweights)){
    asy[,i] <- asy[,i]/sweights[i]
  }
  
  ### 5) Normalizing to internal standard, if wanted
  if (isnm == 'yes'){
    intvals <- asy[telrow,]
    
    for (i in 1:length(intvals)){
      asy[,i] <- asy[,i]/intvals[i]
    }
  }
  
  print('Done with intstd norm')
  
  ##now need to add telrow to rows getting removed
  badrow <- c(badrow, telrow)
  
  ### 6) imputing missing values, ignoring blank
  if (imp == 'yes'){
    asy2 <- asy[-badrow, 1:(ncol(asy)-length(blkcol))]
    mdat <- pl[-badrow,]
    gps <- get_groups(asy2)
    knnlis <- imp_per_group(asy2, gps, mdat)
    mdat2 <- knnlis[[1]]
    knnasy <- knnlis[[2]]
    
    print(paste0('Num of metabolites removed from knn, noise or blank: ', length(badrow)))
    print(length(badrow))
    #badrow <- c(badrow, missrow)
    badrow <- unique(badrow)
    #print('Correlation after knn imputation')
    #print(cor(knnasy))
  }
  
  else{
    ##need to remove bad metabolite rows and blanks
    if (length(blkcol) == 1){
      knnasy <- asy[-badrow, -ncol(asy)]
    } else {
      knnasy <- asy[-badrow, 1:(ncol(asy)-length(blkcol))]
    }
    #print(str(knnasy))
    gps <- get_groups(knnasy)
    #knnasy <- asy[-badrow,]
    mdat2 <- pl[-badrow,]
  }

  write.table(cor(knnasy), file = paste0(args[1], '_cor_after_blank_itsd_smpwght.tab'), sep = '\t', quote = F)
  
  ##7 removing non-correlated replicates
  
  ##getting list of colnames per group
  h <- list()
  tmp <-t(knnasy)
  tmp$groups <- gps
  for(gp in 1:length(unique(gps))){
    name <- unique(gps)[gp]
    h[[unique(gps)[gp]]] <- which(grepl(paste0('^',name), tmp$groups))
  }
  
  gdc <- c()
  for(i in 1:length(h)){
    cls <- h[[i]]
    corr <- cor(knnasy[,cls])
    
    ##which columns in h[[i]] have at least 3 correlations < 0.7
    fbad <- which(unname(colSums(corr < 0.7)) > 2)
    
    if (length(fbad) >= 1) {
      fgood <- h[[i]][-fbad]
    }
    else fgood <- h[[i]]
    gdc <- c(gdc, fgood)
  }
  
  knnasy <- knnasy[,gdc]
  newh <- list(h[[1]][h[[1]] %in% gdc], h[[2]][h[[2]] %in% gdc])
  newasy <- list(asycols[newh[[1]]], asycols[newh[[2]]])
  
  write.table(cor(knnasy), file = paste0(args[1], '_cor_low_cor_reps_removed.tab'), sep = '\t', quote = F)
  
  ##calculating fold change before normalization (Because vsn changes fold change)
  fldc <- log2(rowMeans(knnasy[,newh[[2]]]) / rowMeans(knnasy[,newh[[1]]]))
  
  ### 8) VSN (if wanted)
  ##ignoring rows that did not make it
  if (vsnn == 'yes'){
    vsnasy <- justvsn(knnasy)
    metadat <- mdat2[, -c(asycols)]
    print('Correlation after vsn')
    print(cor(vsnasy))
    
    return(list(newasy, data.frame(metadat[,1:32], vsnasy, fldc)))
  }
  
  else {
    metadat <- mdat2[, -c(asycols)]
    return(list(newasy, data.frame(metadat[,1:32], knnasy, fldc)))
  }
}

get_diff_exp <- function(pl, ctrl, stress) {
  #'Ctrl and bac are both data frames with metabolite peak areas
  #'Ctrl has peak areas for control and bac has areas for infected
  #'pl is said peaklist
  #'This function finds parametric and non-param p-values for
  #'each metabolite, then does FDR with Benjamani-Hochberg, with 
  #'adjusted p-value of 0.05
  #'
  #'Outputs: vector of rows containing differentially expressed
  #'metabolites for this species.
  
  type <- data.frame('group' = c(rep('Control', ncol(ctrl)), rep('Stress', ncol(stress))))
  wpvals <- c()
  tpvals <- c()
  fc <- c()
  
  ##getting Wilcoxon and t-test p-values for each metabolite
  for (metab in 1:nrow(ctrl)){
    ctm <- mean(as.numeric(ctrl[metab,]))
    if (ctm == 0) {ctm <- 1e-30}
    stm <- mean(as.numeric(stress[metab,]))
    ifelse(stm == 0, logfc <- 0, logfc <- log2(stm/ctm))
    if (logfc == Inf) {print(stm, ctm)}
    fc <- c(fc, logfc)
    vals <- as.numeric(unlist(c(ctrl[metab,], stress[metab,])))
    ndf <- cbind(vals, type)
    
    ##getting wilcoxon p-values for this metab 
    wpvals <- c(wpvals, wilcox.test(vals ~ group, data = ndf)[[3]])
    tpvals <- c(tpvals, t.test(vals ~ group, data = ndf)[[3]])
  }
  
  ##finding Benjamani-Hochberg fdr corrected metabs
  diffmetab <- data.frame(pl, pval = tpvals, wpval = wpvals,
                          adj_pval = p.adjust(tpvals, method = 'fdr'), 
                          nfold_change = fc, ofc = pl[,ncol(pl)])
  
  ##which metabolites are Diff abund?
  low <- which((diffmetab$pval <= 0.05) &
                 abs(diffmetab$ofc) >= 2 & 
                 (diffmetab$wpval <= 0.1 | diffmetab$wpval == min(diffmetab$wpval)))
  
  ##if there are any diff expressed metabolites before fdr, are there any after fdr?
  if (length(low) > 0){
    print('There are diff abund metabs before FDR!')
    
    ##setting up pvalue for adj p-values. If there are no/very few metabolites
    ##diff expressed with p-value at or lower than pct, increase it by 0.05
    pct <- 0.05
    adjlow <- which((diffmetab$adj_pval <= pct) & abs(diffmetab$ofc) >= 2 
                    & diffmetab$wpval <= 0.1)
    
    if (length(adjlow) >= 1) {
      print('Got some FDR Diff abund metabolites! ')
      print(length(adjlow))
      out <- data.frame(cbind(pl[adjlow,], diffmetab[adjlow,]))
      return(list(diffmetab, out))
    }else {
      print('No FDR Diff abund metabolites :(')
      out <- data.frame(cbind(pl[low,], diffmetab[low,]))
      return(list(diffmetab, out))
    }
    
  }else {print('No diff abund metabolites at all! :(((')}
  
  return(list(diffmetab, data.frame(pl, diffmetab)))
}

                #####################
                #### Coding time ####
                #####################


##parsing args
args <- commandArgs(T)
print(args)
if (length(args) != 16){
  stop('ARGS: 1) input peak area file, tsv 2) index of metab ID 3) indices of 
       treated columns, separated by a comma 4) indices of control columns
       5) index of column with blank(s) if multiple sepatate with comma
       6) vector of sample weights, separated by comma
       7) column with s/n information 8) s/n cutoff 
       9) cutoff for removal of rows with missing values 10) impute? yes OR no
       11) ID of internal standard 12) Do you want vsn? yes OR no
       13) Do you want to normalize to an internal standard? yes OR no
       14) S/N matrix file 15) threshold for S/N value per individual sample
       16) blank groups. If multiple blanks, input as:
       [blank1-asycol1,asycol2][blank2-asycol4,asycol5]')
}


peakarea <- read.table(args[1], sep = '\t', header = T, stringsAsFactors = F, quote = "")
idcol <- as.numeric(args[2])
if (grepl(',', args[3], fixed = T)){
  treatcol <- as.numeric(unlist(strsplit(args[3], ',')))
}else treatcol <- 0
ctrlcol <- as.numeric(unlist(strsplit(args[4], ',')))
##are there multiple blanks?
if (grepl(',', args[5], fixed = T)) {
  blnkcol <- as.numeric(unlist(strsplit(args[5], ',', fixed = T)))
} else blnkcol <- as.numeric(args[5])

swt <- as.numeric(unlist(strsplit(args[6], ',')))
snc <- as.numeric(args[7])
sncut <- as.numeric(args[8])
knncut <- as.numeric(args[9])
immp <- as.character(args[10])
tlid <- as.numeric(args[11])
varsn <- as.character(args[12])
isd <- as.character(args[13])
snmat <- read.table(args[14], sep = '\t', header = T, stringsAsFactors = F, quote = "")
snt <- as.numeric(args[15])
blankgroups <- args[16]

##normalizing
npeakarea_l <- vs_norm(pl = peakarea, ccols = ctrlcol, tcols = treatcol,
                       blkcol = blnkcol, sweights = swt, sncol = snc, 
                       snthresh = sncut, cutoff = knncut, imp = immp, telid = tlid,
                       vsnn = varsn, isnm = isd, snmatt = snmat, snsampt = snt, blkgrp = blankgroups)

npeakarea <- as.data.frame(npeakarea_l[[2]])

##writing out normalized peak areas, of all surviving metabolites
outname <- c(args[1])
if (varsn == 'yes'){
  outname <- c(outname, '_VSN_')
}

if (isd == 'yes'){
  outname <- c(outname, '_IS_')
}

write.table(npeakarea, file = paste(c(outname, '_filtered_normalized_metabs.tab'), collapse = ''), 
            sep = '\t', row.names = F, quote = F)

##now getting differentially abundant peaks (if there are conditions)
print(treatcol)
if (treatcol != 0){
  newgroups <- get_groups(npeakarea[,-c(1:32, ncol(npeakarea))])
  
  leaves <- newgroups[which(grepl('Leaf', newgroups))]
  roots <- newgroups[which(grepl('Root', newgroups))]
  for (tissue in list(leaves, roots)){
    if (length(tissue) >0){
      tname <- unique(tissue[grepl('Ctrl.', tissue)])
      ctrlcol <- which(grepl(tname, colnames(npeakarea)))
      treat <- tissue[-which(grepl('Ctrl', tissue))]
      uniq <- unique(treat)
      for (cond in uniq){
        print(cond)
        ttreat <- which(grepl(paste0('^', cond), colnames(npeakarea)))
        if (length(ttreat) > 1 & length(ctrlcol) > 1){
          both <- get_diff_exp(npeakarea, ctrl = npeakarea[,ctrlcol], stress = npeakarea[,ttreat])
          allpval <- both[[1]]; diff <- both[[2]] 
          write.table(diff, file = paste(c(outname, '_FDR_diff_', cond, '.tab'), collapse = ''), sep = '\t', row.names = F, quote = F)
          write.table(allpval, file = paste(c(outname, '_all_pval_fc_', cond, '.tab'), collapse = ''), sep = '\t', row.names = F, quote = F)
        }else {print('Not enough replicates to find Differentially Abundant metabolites!')}
      }
    }
  }
}
print('Done! Wow!')
