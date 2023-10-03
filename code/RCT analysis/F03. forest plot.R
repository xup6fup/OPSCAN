
library(magrittr)
library(ggplot2)
library(scales) 
library(cowplot)
library(gridExtra)

# 0. Data path

plot_path <-  './result/Fig 03.png'
data_path <-  './data/RCT analysis/RCT data.RData'

# 1. Load data

load(data_path)

follow_data <- follow_data[follow_data[,'Group'] %in% c('Screening', 'Control'),]
follow_data[,'Group'] <- as.character(follow_data[,'Group']) %>% factor(., levels = c('Screening', 'Control'))

# 2. Data list

data_list <- list()

# Group 1 DXA

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[DXA]']
data_list[[1]] <- current_data

# Group 1 DXA - Trial

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[DXA]']
current_data[current_data[,'with_response'] %in% 0,'status'] <- 0
data_list[[2]] <- current_data

# Group 1 DXA - Physicians

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[DXA]']
current_data[current_data[,'with_response'] %in% 1,'status'] <- 0
data_list[[3]] <- current_data

# Group 1 DXA - L-spine

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[DXA-L]']
data_list[[4]] <- current_data

# Group 1 DXA - Hip

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[DXA-H]']
data_list[[5]] <- current_data

# Group 2 osteoporosis

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[lowBMD]']
data_list[[6]] <- current_data

# Group 2 osteoporosis - L-spine

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[lowBMD-L]']
data_list[[7]] <- current_data

# Group 2 osteoporosis - Hip

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[lowBMD-H]']
data_list[[8]] <- current_data

# Group 3 All treatment

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[BMD-drug]']
data_list[[9]] <- current_data

# Group 4 Denosumab treatment

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[BMD-drug[Denosumab]]']
data_list[[10]] <- current_data

# Group 4 Romosozumab treatment

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[BMD-drug[Romosozumab]]']
data_list[[11]] <- current_data

# Group 3 Bisphosponate treatment

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[BMD-drug[Bisphosponate]]']
data_list[[12]] <- current_data

# Group 4 Parathyroid hormone treatment

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[BMD-drug[Parathyroid hormone]]']
data_list[[13]] <- current_data

# Group 4 SERM treatment

current_data <- follow_data
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[BMD-drug[SERM]]']
data_list[[14]] <- current_data

# 3. Forest plot

## Parameters

data_position <- c(16, 15, 14, 13, 12, 10, 9, 8, 6, 5, 4, 3, 2, 1)

data_names.1 <- expression(paste('All DXA examination'),
                           paste('New-onset osteoporosis (primary endpoints)'),
                           paste('All anti-osteoporotic medications'))

data_position.1 <- c(16, 10, 6)

data_names.2 <- expression(paste('DXA examination by OPSCAN'),
                           paste('DXA examination by physicians'),
                           paste('DXA examination for L-spine'),
                           paste('DXA examination for hip'),
                           paste('New-onset osteoporosis for L-spine'),
                           paste('New-onset osteoporosis for hip'),
                           paste('Denosumab'),
                           paste('Romosozumab'),
                           paste('Bisphosponate'),
                           paste('Parathyroid hormone'),
                           paste('Selective estrogen-receptor modulator'))

data_position.2 <- c(15, 14, 13, 12, 9, 8, 5, 4, 3, 2, 1)

group_names <- c('Examination', 'Result', 'Medication')

group_position <- c(17, 11, 7)

## For Odds ratio

OR_data_list <- list()

for (j in 1:length(data_list)) {
  
  sub_data <- data_list[[j]]
  
  glm_model <- try(summary(glm(as.formula(paste0('status ~ x')), data = sub_data, family = 'binomial')), silent = TRUE)
  
  n0 <- nrow(sub_data[sub_data[,'x'] %in% 'Control',])
  n1 <- nrow(sub_data[sub_data[,'x'] %in% 'Screening',])
  
  e0 <- sum(sub_data[sub_data[,'x'] %in% 'Control','status'])
  e1 <- sum(sub_data[sub_data[,'x'] %in% 'Screening','status'])
  
  p0 <- e0 / n0
  p1 <- e1 / n1
  
  if ('try-error' %in% class(glm_model) | (e0 %in% 0 & e1 %in% 0)) {
    
    y.OR <- NA
    se.OR <- NA
    OR_txt <- ''
    
  } else {
    
    y.OR <- glm_model[['coefficients']][2,1]
    se.OR <- glm_model[['coefficients']][2,2]
    
    OR_txt <- paste0(formatC(exp(y.OR), 2, format = 'f'), ' (',
                     formatC(exp(y.OR - qnorm(0.975) * se.OR), 2, format = 'f'), ', ',
                     formatC(exp(y.OR + qnorm(0.975) * se.OR), 2, format = 'f'), ')')
    
  }
  
  OR_data_list[[j]] <- data.frame(x = data_position[j],
                                  n = paste0(e0 + e1, '/', n0 + n1),
                                  n.1 = paste0(e1, '/', n1, ' (', formatC(p1 * 100, 1, format = 'f'), '%)'),
                                  n.0 = paste0(e0, '/', n0, ' (', formatC(p0 * 100, 1, format = 'f'), '%)'),
                                  y = exp(y.OR), yse = se.OR,
                                  ylo = exp(y.OR - qnorm(0.975) * se.OR),
                                  yhi = exp(y.OR + qnorm(0.975) * se.OR),
                                  txt.OR = OR_txt,
                                  col = '#000000',
                                  size = 1,
                                  stringsAsFactors = FALSE)
  
}

OR_data <- do.call('rbind', OR_data_list)

Inf_pos <- grep('(0.00, Inf)', OR_data[,'txt.OR'], fixed = TRUE)

if (length(Inf_pos) > 0) {
  
  y_min <- OR_data[-Inf_pos,'ylo'] %>% min()
  y_max <- OR_data[-Inf_pos,'yhi'] %>% max()
  
  OR_data[Inf_pos,'ylo'] <- y_min
  OR_data[Inf_pos,'yhi'] <- y_max
  
  positive_pos <- which(OR_data[Inf_pos,'y'] > 1)
  negative_pos <- which(OR_data[Inf_pos,'y'] <= 1)
  
  OR_data[Inf_pos[negative_pos],'y'] <- y_min
  OR_data[Inf_pos[positive_pos],'y'] <- y_max
  
  OR_data[Inf_pos[negative_pos],'txt.OR'] <- '0.00 (0.00, Inf)'
  OR_data[Inf_pos[positive_pos],'txt.OR'] <- 'Inf (0.00, Inf)'
  
}

y_lim <- c(min(OR_data[,'ylo'], na.rm = TRUE), max(OR_data[,'yhi'], na.rm = TRUE))
if (y_lim[1] > 0.9) {y_lim[1] <- 0.9}

OR_data[which(OR_data[,'ylo'] < y_lim[1]),'ylo'] <- y_lim[1]
OR_data[which(OR_data[,'yhi'] > y_lim[2]),'yhi'] <- y_lim[2]

OR_data[which(OR_data[,'y'] < y_lim[1]),'y'] <- y_lim[1]
OR_data[which(OR_data[,'y'] > y_lim[2]),'y'] <- y_lim[2]

diff_scale <- (y_lim[2] / y_lim[1])

y_revised_lim <- y_lim
y_revised_lim[1] <- y_revised_lim[1] * diff_scale^(-4.8)
y_revised_lim[2] <- y_revised_lim[2] * diff_scale^1.1

# Index

index_name <- 'OR (95% CI)'

# Plotting

forest_p <- ggplot(OR_data, aes(x = x, y = y, ymin = ylo, ymax = yhi))

forest_p <- forest_p + annotate(geom = 'segment', x = 0.4, y = 1, xend = max(group_position) + 1.5, yend = 1, colour = "black", size = 0.4, linetype = 'dotted')
forest_p <- forest_p + annotate(geom = 'line', x = c(0.5, 0.5), y = y_lim, colour = "black", size = 0.4)

scale_y <- c(0.1, 0.2, 0.5, 1, 2, 4, 8, 16, 32, 64)
scale_y <- scale_y[scale_y > y_lim[1] & scale_y < y_lim[2]]

for (j in 1:length(scale_y)) {
  
  forest_p <- forest_p + annotate(geom = 'segment', x = 0.3, y = scale_y[j], xend = 0.7, yend = scale_y[j], colour = "black", size = 0.4)
  forest_p <- forest_p + annotate(geom = "text", x = -0.3, y = scale_y[j], label = formatC(scale_y[j], 1, format = 'f'), size = 2.5, color = "black", angle = 45)
  
}

forest_p <- forest_p + annotate(geom = "text", x = -1.2, y = exp(mean(log(y_lim))), label = 'Screening vs. Control', size = 4.5, color = "black", fontface = "bold")

forest_p <- forest_p + geom_pointrange(shape = 15, size = 0.4, fill = 'white', position = position_dodge(width = 0.1), fatten = 5)
forest_p <- forest_p + geom_errorbar(size = 0.4, width = 0.4, cex = 1)
forest_p <- forest_p + coord_flip() + xlab('') + ylab('')
forest_p <- forest_p + scale_y_continuous(limits = y_revised_lim, trans = 'log2', breaks = NULL)
forest_p <- forest_p + theme(panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank(),
                             axis.text.x = element_blank(),
                             axis.text.y = element_blank(),
                             axis.line.x = element_blank(),
                             axis.line.y = element_blank(),
                             axis.title.x = element_blank(),
                             axis.ticks.y = element_blank(),
                             legend.position = "none")

forest_p <- forest_p + annotate(geom = "text", x = group_position, y = y_lim[1] * diff_scale^(-4.8), label = group_names, size = 4, color = "black", fontface = "bold", hjust = 0)
forest_p <- forest_p + annotate(geom = "text", x = data_position.1, y = y_lim[1] * diff_scale^(-4.6), label = data_names.1, size = 4, hjust = 0)
forest_p <- forest_p + annotate(geom = "text", x = data_position.2, y = y_lim[1] * diff_scale^(-4.4), label = data_names.2, size = 4, hjust = 0)
forest_p <- forest_p + annotate(geom = "text", x = OR_data[,'x'], y = y_lim[1] * diff_scale^(-1.8), label = OR_data[,'n.1'], size = 3.5)
forest_p <- forest_p + annotate(geom = "text", x = OR_data[,'x'], y = y_lim[1] * diff_scale^(-0.6), label = OR_data[,'n.0'], size = 3.5)
forest_p <- forest_p + annotate(geom = "text", x = OR_data[,'x'], y = y_lim[2] * diff_scale^(0.7), label = OR_data[,'txt.OR'], size = 3.5)

forest_p <- forest_p + annotate(geom = "text", x = max(group_position) + 1, y = y_lim[1] * diff_scale^(-1.8), label = "event/n (%)", size = 4, color = "black")
forest_p <- forest_p + annotate(geom = "text", x = max(group_position) + 2, y = y_lim[1] * diff_scale^(-1.8), label = "Screening", size = 4.5, color = "black", fontface = "bold")
forest_p <- forest_p + annotate(geom = "text", x = max(group_position) + 1, y = y_lim[1] * diff_scale^(-0.6), label = "event/n (%)", size = 4, color = "black")
forest_p <- forest_p + annotate(geom = "text", x = max(group_position) + 2, y = y_lim[1] * diff_scale^(-0.6), label = "Control", size = 4.5, color = "black", fontface = "bold")
forest_p <- forest_p + annotate(geom = "text", x = max(group_position) + 1.5, y = y_lim[2] * diff_scale^(0.7), label = index_name, size = 4.5, color = "black", fontface = "bold")

# 4. Merge plotting

if (!dir.exists(dirname(plot_path))) {dir.create(dirname(plot_path), recursive = TRUE)}

final_p <- ggdraw()
final_p <- final_p + draw_plot(forest_p, x = -0.06, y = -0.03, width = 1.08, height = 1.07)

png(plot_path, width = 800, height = 560)

print(final_p)

dev.off()
