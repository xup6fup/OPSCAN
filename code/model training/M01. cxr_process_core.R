library(abind)

cxr_process_core <- function(img_paths, img_crop = 224, VALID_CROP = NULL,
                             col_flip = TRUE, add_sd = 0.05, mul_sd = 0.1, pow_sd = 0.1) {
  
  batch_IMG_LIST <- list()
  
  for (i in 1:length(img_paths)) {
    
    load(img_paths[i])
    
    resize_img <- array(resize_img, dim = c(dim(resize_img), 3))
    
    row_df <- dim(resize_img)[1] - img_crop
    col_df <- dim(resize_img)[2] - img_crop
    
    if (is.null(VALID_CROP)) {
      
      random.row <- sample(0:row_df, 1)
      random.col <- sample(0:col_df, 1)
      
      batch_IMG_LIST[[i]] <- resize_img[random.row+1:img_crop,random.col+1:img_crop,,drop=FALSE]
      
    } else {
      
      if (VALID_CROP %in% 1:10) {
        
        if (VALID_CROP %in% c(1, 3, 6, 8)) {
          
          fixed.row <- 0
          
        } else if (VALID_CROP %in% c(2, 4, 7, 9)) {
          
          fixed.row <- row_df
          
        } else if (VALID_CROP %in% c(5, 10)) {
          
          fixed.row <- floor(row_df / 2)
          
        }
        
        if (VALID_CROP %in% c(1, 2, 6, 7)) {
          
          fixed.col <- 0
          
        } else if (VALID_CROP %in% c(3, 4, 8, 9)) {
          
          fixed.col <- col_df
          
        } else if (VALID_CROP %in% c(5, 10)) {
          
          fixed.col <- floor(col_df / 2)
          
        }
        
        if (VALID_CROP %in% c(1:5)) {
          
          fixed.flip <- FALSE
          
        } else {
          
          fixed.flip <- TRUE
          
        }
        
        batch_IMG_LIST[[i]] <- resize_img[fixed.row+1:img_crop,fixed.col+1:img_crop,,drop=FALSE]
        
      } else {
        
        stop('VALID_CROP should be ranged 1 to 10!')
        
      }
      
    }
    
  }
  
  batch_img_array <- abind(batch_IMG_LIST, along = 4)
  
  if (!is.null(VALID_CROP)) {
    
    if (fixed.flip) {
      
      batch_img_array <- batch_img_array[,dim(batch_img_array)[2]:1,,,drop = FALSE]
      
    }
    
    batch_img_array <- mx.nd.array(batch_img_array)
    
  } else {
    
    if (col_flip) {
      
      if (sample(c(TRUE, FALSE), 1)) {
        
        batch_img_array <- batch_img_array[,dim(batch_img_array)[2]:1,,,drop = FALSE]
        
      }
      
    }
    
    batch_img_array <- mx.nd.array(batch_img_array)
    
    if (add_sd != 0) {
      
      add_item <- mx.nd.array(array(runif(length(img_paths) * 1, min = -add_sd, max = add_sd), dim = c(1, 1, 1, length(img_paths))))
      batch_img_array <- mx.nd.broadcast.add(batch_img_array, add_item)
      
    }
    
    if (mul_sd != 0) {
      
      mul_item <- mx.nd.array(array(exp(runif(length(img_paths) * 1, min = -mul_sd, max = mul_sd)), dim = c(1, 1, 1, length(img_paths))))
      batch_img_array <- mx.nd.broadcast.mul(batch_img_array, mul_item)
      
    }
    
    if (pow_sd != 0) {
      
      pow_item <- mx.nd.array(array(exp(runif(length(img_paths) * 1, min = -pow_sd, max = pow_sd)), dim = c(1, 1, 1, length(img_paths))))
      batch_img_array <- mx.nd.broadcast.power(batch_img_array, pow_item)
      
    }
    
    batch_img_array <- mx.nd.broadcast.maximum(batch_img_array, mx.nd.array(array(0, dim = c(1, 1, 1, 1))))
    batch_img_array <- mx.nd.broadcast.minimum(batch_img_array, mx.nd.array(array(1, dim = c(1, 1, 1, 1))))
    
  }
  
  return(batch_img_array)
  
}
