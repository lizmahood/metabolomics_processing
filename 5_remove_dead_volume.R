argg <- commandArgs(T)

if (length(argg) != 3){
  stop('ARGS: 1) PeakArea file, 2) corresponding SN_matrix, 3) PeakHeight file')
}

innam <- argg[1]
snnam <- argg[2]

#infil <- read.table(innam, sep = '\t', header = T, quote = "", fill = NA)
infil <- read.csv(argg[1], sep = "\t", header = TRUE,comment.char = "", check.names = T)
#snmat <- read.table(snnam, sep = '\t', header = T, quote = "", fill = NA)
snmat <- read.csv(argg[2], sep = "\t", header = TRUE,comment.char = "", check.names = T)
#pkheight <- read.table(argg[3], sep = '\t', header = T, quote = "", fill = NA)
pkheight <- read.csv(argg[3], sep = "\t", header = TRUE,comment.char = "", check.names = T)
print('Number of dead volume metabolites:')

deadvol <- which(infil$Average.Rt.min. <= 1.5)
print(length(deadvol))
outfil <- infil[-deadvol,]
outsn <- snmat[-deadvol,]
outheight <- pkheight[-deadvol,]
write.table(outfil, file = paste0(innam, '_deadvol_removed.tab'), row.names = F, quote = F, sep = '\t')
write.table(outsn, file = paste0(snnam, '_deadvol_removed.tab'), row.names = F, quote = F, sep = '\t')
write.table(outheight, file = paste0(argg[3], '_deadvol_remved.tab'), row.names = F, quote = F, sep = '\t')