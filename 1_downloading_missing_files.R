library('googledrive')
library('purrr')
setwd('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/')

##folder of Brachy metabolite data
folderurl <- 'https://drive.google.com/drive/folders/1xVbTh1XExLfECuGAc9tDYwMb4wyPIlQb'

folder <- drive_get(as_id(folderurl))

##what mzml files are in the folder?
mzml_files <- drive_ls(folder, pattern = "mzML")

##which mzml files do we already have?
localmzml <- list.files('D:/MS_Data/BrachyMetabolites/original_files/', pattern = 'mzML', full.names = F)

##which don't we have?
need <- mzml_files

##downloading them now
walk(need$id, ~ drive_download(as_id(.x)))

##Done