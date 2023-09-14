
library(pROC)

source('./code/model training/M01. cxr_process_core.R')
source('./code/model training/M02. architecture.R')

# 0. Settings

model_name <- 'OPSCAN'

in_path <- './data/model training/label.csv'
RData_dir <- './data/model training/RData/'

batch_size <- 32
valid_type <- 1:10

used_gpu <- mx.gpu(3)

# 1. Read data

ori_dat <- read.csv(in_path)

set.seed(0)

train_pos <- sample(1:nrow(ori_dat), nrow(ori_dat) * 0.6)

train_subdat <- ori_dat[train_pos,]
valid_subdat <- ori_dat[-train_pos,]

# 2. executor

my_executor <- mx.simple.bind(symbol = loss_symbol,
                              data = c(224, 224, 3, batch_size), label = c(1, batch_size),
                              ctx = used_gpu, grad.req = "write")

pred_executor <- mx.simple.bind(symbol = logistic_pred,
                                data = c(224, 224, 3, batch_size), 
                                ctx = used_gpu)

# 3. Training

## Initialize

mx.set.seed(0)
new_arg <- mxnet:::mx.model.init.params(symbol = loss_symbol,
                                        input.shape = list(data = c(224, 224, 3, batch_size), label = c(1, batch_size)),
                                        output.shape = NULL,
                                        initializer = mxnet:::mx.init.uniform(0.01),
                                        ctx = mx.cpu())

for (i in 1:length(new_arg$arg.params)) {
  pos <- which(names(res_model$arg.params) == names(new_arg$arg.params)[i])
  if (length(pos) == 1) {
    if (all.equal(dim(res_model$arg.params[[pos]]), dim(new_arg$arg.params[[i]])) == TRUE) {
      new_arg$arg.params[[i]] <- res_model$arg.params[[pos]]
    }
  }
}

for (i in 1:length(new_arg$aux.params)) {
  pos <- which(names(res_model$aux.params) == names(new_arg$aux.params)[i])
  if (length(pos) == 1) {
    if (all.equal(dim(res_model$aux.params[[pos]]), dim(new_arg$aux.params[[i]])) == TRUE) {
      new_arg$aux.params[[i]] <- res_model$aux.params[[pos]]
    }
  }
}

mx.exec.update.arg.arrays(my_executor, new_arg$arg.params, match.name = TRUE)
mx.exec.update.aux.arrays(my_executor, new_arg$aux.params, match.name = TRUE)

## Multi-LR strategy

current_epoch <- 0L
train_loss <- NULL
valid_auc <- NULL

t0 <- Sys.time()

for (lr in c(1e-3, 1e-4, 1e-5)) {
  
  my_optimizer <- mx.opt.create(name = "adam", learning.rate = lr, beta1 = 0.9, beta2 = 0.999, epsilon = 1e-08, wd = 1e-3)
  my_updater <- mx.opt.get.updater(optimizer = my_optimizer, weights = my_executor$ref.arg.arrays)
  
  # Training
  
  for (sub_epoch in 1:50) {
    
    current_epoch <- current_epoch + 1L
    batch_loss <- NULL
    
    for (batch in 1:10) {
      
      positive_id <- sample(which(train_subdat[,'OP'] == 1), size = batch_size / 2, replace = TRUE)
      negative_id <- sample(which(train_subdat[,'OP'] == 0), size = batch_size / 2, replace = TRUE)
      used_img_paths <- paste0(RData_dir, gsub('\\.png', '.RData', train_subdat[c(positive_id, negative_id),'NIH_Image_Index']))
      
      current_data <- cxr_process_core(img_paths = used_img_paths)
      current_label <- mx.nd.array(array(rep(1:0, each = batch_size / 2), dim = c(1, batch_size)))
      
      mx.exec.update.arg.arrays(my_executor, arg.arrays = list(data = current_data, label = current_label), match.name = TRUE)
      mx.exec.forward(my_executor, is.train = TRUE)
      mx.exec.backward(my_executor)
      update_args <- my_updater(weight = my_executor$ref.arg.arrays, grad = my_executor$ref.grad.arrays)
      mx.exec.update.arg.arrays(my_executor, update_args, skip.null = TRUE)
      batch_loss <- c(batch_loss, as.array(my_executor$ref.outputs$ce_loss_output))
      
      if (batch %% 5 == 0) {
        
        epoch_time <- as.double(difftime(Sys.time(), t0, units = "secs"))
        epoch_speed <- epoch_time / batch

        message("lr = ", lr, " Epoch [", current_epoch, "] Batch [", batch, "] Speed: ", formatC(epoch_speed, 2, format = "f"),
                " sec/batch loss = ", formatC(mean(batch_loss), format = "f", 4))
        
      }
      
    }
    
    # Overall
    
    epoch_time <- as.double(difftime(Sys.time(), t0, units = "mins"))
    epoch_speed <- epoch_time / current_epoch
    train_loss <- c(train_loss, mean(batch_loss))
    
    message("lr = ", lr, " Epoch [", current_epoch, "] Speed: ", formatC(epoch_speed, 2, format = "f"),
            " min/epoch loss = ", formatC(mean(batch_loss), format = "f", 4))
    
    # Save model
    
    pred_model <- mxnet:::mx.model.extract.model(symbol = logistic_pred, train.execs = list(my_executor))
    mx.model.save(model = pred_model, prefix = paste0('model/', model_name), iteration = current_epoch)
    
    # Validation process
    
    mx.exec.update.arg.arrays(pred_executor, pred_model$arg.params, match.name = TRUE)
    mx.exec.update.aux.arrays(pred_executor, pred_model$aux.params, match.name = TRUE)
    
    current_n_batch <- ceiling(nrow(valid_subdat) / batch_size)
    
    for (k in 1:current_n_batch) {
      
      sample_idx <- 1:batch_size + (k - 1) * batch_size
      sample_idx[sample_idx > nrow(valid_subdat)] <- nrow(valid_subdat)
      used_img_paths <- paste0(RData_dir, gsub('\\.png', '.RData', valid_subdat[sample_idx,'NIH_Image_Index']))
      
      for (l in 1:length(valid_type)) {
        
        BATCH_CXR_ARRAY <- cxr_process_core(img_paths = used_img_paths, VALID_CROP = valid_type[l],
                                            col_flip = FALSE, add_sd = 0, mul_sd = 0, pow_sd = 0)
        
        mx.exec.update.arg.arrays(pred_executor, arg.arrays = list(data = BATCH_CXR_ARRAY), match.name = TRUE)
        mx.exec.forward(pred_executor, is.train = FALSE)
        
        current_pred <- as.array(pred_executor[['ref.outputs']][[1]])
        valid_subdat[sample_idx,paste0('pred.', l)] <- current_pred %>% as.numeric()
        
      }
      
    }
    
    valid_subdat[,'final_pred'] <- apply(valid_subdat[,paste0('pred.', 1:length(valid_type))], 1, mean)
    roc_valid <- roc(OP ~ final_pred, data = valid_subdat)
    valid_auc <- c(valid_auc, roc_valid[['auc']])
    
    message("valid-auc = ", formatC(roc_valid[['auc']], format = "f", 4))
    
    save(train_loss, valid_auc, file = 'model/logger.RData')
    
    if (which.max(valid_auc) < current_epoch) {break}
    
  }
  
}

## Save the best model

best_iter <- which.max(valid_auc)

pred_model <- mx.model.load(prefix = 'model/OPSCAN', iteration = best_iter)
mx.model.save(pred_model, prefix = 'model/OPSCAN', iteration = 0)

