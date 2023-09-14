
library(magrittr)
library(mxnet)

res_model <- mx.model.load(prefix = "model/resnet-18", iteration = 0)
all_layers <- res_model$symbol$get.internals()
flatten0_output <- which(all_layers$outputs == 'flatten0_output') %>% all_layers$get.output()

fc1 <- mx.symbol.FullyConnected(data = flatten0_output, num_hidden = 1, name = 'fc1')
logistic_pred <- mx.symbol.sigmoid(data = fc1, name = 'logistic_pred')

label_sym <- mx.symbol.Variable('label')

eps <- 1e-8
ce_loss_pos <- mx.symbol.broadcast_mul(mx.symbol.log(logistic_pred + eps), label_sym)
ce_loss_neg <- mx.symbol.broadcast_mul(mx.symbol.log(1 - logistic_pred + eps), 1 - label_sym)
ce_loss_mean <- 0 - mx.symbol.mean(ce_loss_pos + ce_loss_neg)
loss_symbol <- mx.symbol.MakeLoss(ce_loss_mean, name = 'ce_loss')
