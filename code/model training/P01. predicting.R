
source('./code/model training/M01. cxr_process_core.R')

# 0. Settings

batch_size <- 32
valid_type <- 1:10

used_gpu <- mx.gpu(3)

model_name <- 'OPSCAN'

in_path <- './data/model training/label.csv'
RData_dir <- './data/model training/RData/'

# 1. Read data

ori_dat <- read.csv(in_path)

# 2. Load model

pred_model <- mx.model.load(prefix = paste0('model/', model_name), iteration = 0)

pred_executor <- mx.simple.bind(symbol = pred_model$symbol,
                                data = c(224, 224, 3, batch_size), 
                                ctx = used_gpu)

# Validation process

mx.exec.update.arg.arrays(pred_executor, pred_model$arg.params, match.name = TRUE)
mx.exec.update.aux.arrays(pred_executor, pred_model$aux.params, match.name = TRUE)

current_n_batch <- ceiling(nrow(ori_dat) / batch_size)

pb <- txtProgressBar(max = current_n_batch, style = 3)

for (k in 1:current_n_batch) {
  
  sample_idx <- 1:batch_size + (k - 1) * batch_size
  sample_idx[sample_idx > nrow(ori_dat)] <- nrow(ori_dat)
  used_img_paths <- paste0(RData_dir, gsub('\\.png', '.RData', ori_dat[sample_idx,'NIH_Image_Index']))
  
  for (l in 1:length(valid_type)) {
    
    BATCH_CXR_ARRAY <- cxr_process_core(img_paths = used_img_paths, VALID_CROP = valid_type[l],
                                        col_flip = FALSE, add_sd = 0, mul_sd = 0, pow_sd = 0)
    
    mx.exec.update.arg.arrays(pred_executor, arg.arrays = list(data = BATCH_CXR_ARRAY), match.name = TRUE)
    mx.exec.forward(pred_executor, is.train = FALSE)
    
    current_pred <- as.array(pred_executor[['ref.outputs']][[1]])
    ori_dat[sample_idx,paste0('pred.', l)] <- current_pred %>% as.numeric()
    
  }
  
  setTxtProgressBar(pb, k)
  
}

close(pb)

ori_dat[,'final_pred'] <- apply(ori_dat[,paste0('pred.', 1:length(valid_type))], 1, mean)
