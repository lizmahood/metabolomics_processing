library('googledrive')
library('purrr')

argg <- commandArgs(T)

if (length(argg) != 2){
  stop('ARGS: 1) directory for files to be downloaded to
       2) Pattern for files to include, in regex')
}

setwd(argg[1])

##folder of Brachy metabolite data
folderurl <- 'https://drive.google.com/drive/folders/1xVbTh1XExLfECuGAc9tDYwMb4wyPIlQb'

folder <- drive_get(as_id(folderurl))

##what mzml files are in the folder?
ptrn <- paste0(argg[2], '.*mzML')
print(ptrn)
mzml_files <- drive_ls(folder, pattern = ptrn)

##which mzml files do we already have?
localmzml <- list.files('D:/MS_Data/BrachyMetabolites/original_files/', pattern = 'mzML', full.names = F)

##which don't we have?
need <- mzml_files

##downloading them now
walk(need$id, ~ drive_download(as_id(.x)))

##Done
print('Done!')