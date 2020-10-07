##Script for moving around mzml files

mzml <- list.files('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/', pattern = 'mzML', full.names = T)

sym <- mzml[which(grepl('_Sym-', mzml, fixed = T))]
hydro <- mzml[which(grepl('_Hydro-', mzml, fixed = T))]
tissue <- mzml[which(grepl('_Tis-', mzml, fixed = T))]

dir.create('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/Symbiosis/NEG')
dir.create('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/Symbiosis/POS')
dir.create('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/Symbiosis/FPS')

for (fil in sym){
  org_p <- fil
  if (grepl('FPS_MS1', fil, fixed = T)){
    ndir <- 'FPS/'
  }
  if (grepl('NEG_MSMS', fil, fixed = T)){
    ndir <- 'NEG/'
  }
  if (grepl('POS_MSMS', fil, fixed = T)){
    ndir <- 'POS/'
  }
  new_p <- paste0('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/Symbiosis/', ndir, basename(fil))
  file.copy(org_p, new_p)
  unlink(fil)
}

dir.create('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/Copper_Heat/NEG')
dir.create('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/Copper_Heat/POS')
dir.create('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/Copper_Heat/FPS')

for (fil in hydro){
  org_p <- fil
  if (grepl('FPS_MS1', fil, fixed = T)){
    ndir <- 'FPS/'
  }
  if (grepl('NEG_MSMS', fil, fixed = T)){
    ndir <- 'NEG/'
  }
  if (grepl('POS_MSMS', fil, fixed = T)){
    ndir <- 'POS/'
  }
  new_p <- paste0('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/Copper_Heat/', ndir, basename(fil))
  file.copy(org_p, new_p)
  unlink(fil)
}

dir.create('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/Tissue_ctrl/NEG')
dir.create('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/Tissue_ctrl/POS')
dir.create('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/Tissue_ctrl/FPS')

for (fil in tissue){
  org_p <- fil
  if (grepl('FPS_MS1', fil, fixed = T)){
    ndir <- 'FPS/'
  }
  if (grepl('NEG_MSMS', fil, fixed = T)){
    ndir <- 'NEG/'
  }
  if (grepl('POS_MSMS', fil, fixed = T)){
    ndir <- 'POS/'
  }
  new_p <- paste0('C:/Users/ehm79/Documents/MS_Data/BrachyMetabolites/original_files/Tissue_ctrl/', ndir, basename(fil))
  file.copy(org_p, new_p)
  unlink(fil)
}

print('Done!')

