
library(magrittr)
library(OpenImageR)

# 0. Settings

in_path <- './data/model training/label.csv'

png_dir <- './data/model training/png/'
RData_dir <- './data/model training/RData/'

# 1. Read data

ori_dat <- read.csv(in_path)

# 2. Processing

pb <- txtProgressBar(max = nrow(ori_dat), style = 3)

for (i in 1:nrow(ori_dat)) {
  
  ori_path <- paste0(png_dir, ori_dat[i,'NIH_Image_Index'])
  prc_path <- paste0(RData_dir, gsub('\\.png$', '.RData', ori_dat[i,'NIH_Image_Index']))
  
  # Resize file
  
  if (!file.exists(prc_path)) {
    
    img <- readImage(ori_path)
    
    if (dim(img)[3] %in% c(3, 4)) {
      
      img <- img[,,1:3] %>% apply(., 1:2, mean)
      
    } else {
      
      img <- array(img, dim = dim(img))
      
    }
    
    target_size <- c(256, 256)
    targer_fold <- dim(img)[1] / dim(img)[2]
    if (targer_fold > 1) {target_size[1] <- target_size[1] * targer_fold} else {target_size[2] <- target_size[2] / targer_fold}
    resize_img <- resizeImage(img, width = target_size[1], height = target_size[2], method = 'bilinear')
    
    if (!dir.exists(dirname(prc_path))) {dir.create(dirname(prc_path), recursive = TRUE)}
    
    save(resize_img, file = prc_path)
    
  }
  
  setTxtProgressBar(pb, i)
  
}

close(pb)
