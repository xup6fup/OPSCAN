
library(magrittr)
library(ggplot2)
library(scales) 
library(cowplot)
library(gridExtra)

# 0. Data path

plot_path <-  './result/Fig 04.png'
data_path <-  './data/RCT analysis/RCT data.RData'

# 1. Load data

load(data_path)

follow_data[,'ISCD indication'] <- factor((follow_data[,'GENDER'] %in% 'male' & follow_data[,'AGE'] >= 70) | (follow_data[,'GENDER'] %in% 'female' & follow_data[,'AGE'] >= 65))
levels(follow_data[,'ISCD indication']) <- c('Not meeting ISCD indication', 'Meeting ISCD indication')

# 2. Data

data_list <- list()

# Group 1 DXA

current_data <- follow_data[follow_data[,'ISCD indication'] %in% 'Not meeting ISCD indication',]
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[DXA]']
data_list[[1]] <- current_data

# Group 1 DXA

current_data <- follow_data[follow_data[,'ISCD indication'] %in% 'Meeting ISCD indication',]
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[DXA]']
data_list[[2]] <- current_data

# Group 2 Low T-score

current_data <- follow_data[follow_data[,'ISCD indication'] %in% 'Not meeting ISCD indication',]
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[lowBMD]']
data_list[[3]] <- current_data

# Group 2 Low T-score

current_data <- follow_data[follow_data[,'ISCD indication'] %in% 'Meeting ISCD indication',]
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[lowBMD]']
data_list[[4]] <- current_data

# Group 3 All treatment

current_data <- follow_data[follow_data[,'ISCD indication'] %in% 'Not meeting ISCD indication',]
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[BMD-drug]']
data_list[[5]] <- current_data

# Group 3 All treatment

current_data <- follow_data[follow_data[,'ISCD indication'] %in% 'Meeting ISCD indication',]
current_data[,'x'] <- factor(current_data[,'Group'], levels = c('Control', 'Screening'))
current_data[,'status'] <- current_data[,'event[BMD-drug]']
data_list[[6]] <- current_data

# 3. Forest plot

# Parameters

data_names <- expression(paste('Men aged ' < 70, ' y/o or women aged ' < 65, ' y/o'),
                         paste('Men aged ' >= 70, ' y/o or women aged ' >= 65, ' y/o'),
                         paste('Men aged ' < 70, ' y/o or women aged ' < 65, ' y/o'),
                         paste('Men aged ' >= 70, ' y/o or women aged ' >= 65, ' y/o'),
                         paste('Men aged ' < 70, ' y/o or women aged ' < 65, ' y/o'),
                         paste('Men aged ' >= 70, ' y/o or women aged ' >= 65, ' y/o'))

group_names <- c('All DXA examination', 'New-onset osteoporosis (primary endpoints)', 'New-onset usage of anti-osteoporotic medications')

int_position <- list(c(1:2), c(3:4), c(5:6))

# Positions

current_pos <- 0L
data_position <- NULL
group_position <- NULL

for (u in length(int_position):1) {
  
  data_position <- c(data_position, 1:length(int_position[[u]]) + current_pos)
  current_pos <- max(data_position) + 1L
  group_position <- c(group_position, current_pos)
  
}

data_position <- sort(data_position, decreasing = TRUE)
group_position <- sort(group_position, decreasing = TRUE)

p_int <- c()

# For loop interaction

for (j in 1:length(int_position)) {
  
  int_data_list <- list()
  
  check_data <- TRUE
  
  for (data_pos in int_position[[j]]) {
    
    int_data <- data_list[[data_pos]]
    
    if (nrow(int_data) == 0) {
      
      check_data <- FALSE
      break
      
    } else {
      
      int_data[,'stratified'] <- length(int_data_list) + 1
      int_data_list[[length(int_data_list) + 1]] <- int_data
      
    }
    
  }
  
  if (!check_data) {p_int <- c(p_int, '')} else {
    
    combine_data <- do.call('rbind', int_data_list)
    combine_data[,'stratified'] <- factor(combine_data[,'stratified'])
    
    int_model <- try(glm(as.formula(paste0('status ~ x + stratified + x:stratified')), data = combine_data, family = 'binomial'), silent = TRUE)
    
    if ('try-error' %in% class(int_model)) {stop('interaction model error')} else {
      
      COEF <- int_model[['coefficients']][-1]
      VCOV <- vcov(int_model)[-1,-1,drop = FALSE]
      
      pos_var <- (length(COEF) - length(levels(combine_data[,'stratified'])) + 2):length(COEF)
      DF <- length(pos_var)
      CHISQ <- t(COEF[pos_var]) %*% solve(VCOV[pos_var,pos_var]) %*% COEF[pos_var]
      PVALUE <- pchisq(CHISQ, df = DF, lower.tail = FALSE)
      
      p_int_txt <- paste0('p-value = ', formatC(PVALUE, 3, format = 'f'))
      if (p_int_txt %in% 'p-value = 0.000') {p_int_txt <- 'p-value < 0.001'}
      p_int <- c(p_int, p_int_txt)
      
    }
    
  }
  
}

# For Hazard ratio

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
y_revised_lim[1] <- y_revised_lim[1] * diff_scale^(-6.8)
y_revised_lim[2] <- y_revised_lim[2] * diff_scale^2.6

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

forest_p <- forest_p + annotate(geom = "text", x = group_position, y = y_lim[1] * diff_scale^(-6.5), label = group_names, size = 4, color = "black", fontface = "bold", hjust = 0)
forest_p <- forest_p + annotate(geom = "text", x = OR_data[,'x'], y = y_lim[1] * diff_scale^(-6.3), label = data_names, size = 4, hjust = 0)
forest_p <- forest_p + annotate(geom = "text", x = OR_data[,'x'], y = y_lim[1] * diff_scale^(-1.8), label = OR_data[,'n.1'], size = 3.5)
forest_p <- forest_p + annotate(geom = "text", x = OR_data[,'x'], y = y_lim[1] * diff_scale^(-0.6), label = OR_data[,'n.0'], size = 3.5)
forest_p <- forest_p + annotate(geom = "text", x = OR_data[,'x'], y = y_lim[2] * diff_scale^(0.7), label = OR_data[,'txt.OR'], size = 3.5)

forest_p <- forest_p + annotate(geom = "text", x = max(group_position) + 1, y = y_lim[1] * diff_scale^(-1.8), label = "event/n (%)", size = 4, color = "black")
forest_p <- forest_p + annotate(geom = "text", x = max(group_position) + 2, y = y_lim[1] * diff_scale^(-1.8), label = "Screening", size = 4.5, color = "black", fontface = "bold")
forest_p <- forest_p + annotate(geom = "text", x = max(group_position) + 1, y = y_lim[1] * diff_scale^(-0.6), label = "event/n (%)", size = 4, color = "black")
forest_p <- forest_p + annotate(geom = "text", x = max(group_position) + 2, y = y_lim[1] * diff_scale^(-0.6), label = "Control", size = 4.5, color = "black", fontface = "bold")
forest_p <- forest_p + annotate(geom = "text", x = max(group_position) + 1.5, y = y_lim[2] * diff_scale^(0.7), label = index_name, size = 4.5, color = "black", fontface = "bold")
forest_p <- forest_p + annotate(geom = "text", x = max(group_position) + 1.5, y = y_lim[2] * diff_scale^(2.1), label = "p for interaction", size = 4.5, color = "black", fontface = "bold")

forest_p <- forest_p + annotate(geom = "text", x = group_position, y = y_lim[2] * diff_scale^(2.1), label = p_int, size = 3, fontface = "italic")

# 4. Scatter plot

data_col <- c('#0000FF40', '#FF000040')

outcome_names <- c('Major osteoporotic risk' = 'major_osteoporotic', 'Hip fracture risk' = 'hip_fracture')
outcome_y_max <- c(100, 50)

line_data_list <- list()

x_seq <- seq(0, 100, by = 0.1)
RR_list <- c(0.25, 0.5, 1, 2, 4)

for (j in 1:length(RR_list)) {
  
  line_data_list[[length(line_data_list) + 1]] <- data.frame(x1 = x_seq, x2 = x_seq * RR_list[j])
  
}

scattor_p_list <- list()

sub_data <- follow_data[follow_data[,'Group'] %in% 'Screening' & follow_data[,'DXA_type'] %in% 2,]

for (i in 1:length(outcome_names)) {
  
  title_name <- paste0('FRAX score (', names(outcome_names)[i], ')')
  
  x_name <- paste0('before_', outcome_names[i])
  y_name <- paste0('after_', outcome_names[i])
  col_name <- 'ISCD indication'
  
  names(data_col) <- levels(sub_data[,col_name])
  
  x_title <- paste0('Risk before DXA (%)')
  y_title <- paste0('Risk after DXA (%)')
  
  x1 <- sub_data[,x_name]
  x2 <- sub_data[,y_name]
  
  scattor_dat <- data.frame(x = x1, y = x2, Stratification = sub_data[,col_name], RR = x2 / x1, Diff = x2 - x1)
  scattor_dat <- scattor_dat[apply(is.na(scattor_dat), 1, sum) == 0,]
  
  RR_means <- tapply(scattor_dat[,'RR'], scattor_dat[,'Stratification'], mean, na.rm = TRUE)
  RR_sds <- tapply(scattor_dat[,'RR'], scattor_dat[,'Stratification'], sd, na.rm = TRUE)
  
  p_val <- t.test(scattor_dat[,'RR'] ~ scattor_dat[,'Stratification'])[['p.value']]
  p_val_txt <- paste0('p = ', formatC(p_val, digit = 3, format = 'f'))
  if (p_val_txt == 'p = 0.000') {p_val_txt <- 'p < 0.001'}
  
  my_txt.0 <- 'Average of RR:'
  
  my_txt <- paste0(formatC(RR_means, format = 'f', 2),
                   ' (+/-) ',
                   formatC(RR_sds, format = 'f', 2), '\n(',
                   levels(scattor_dat[,'Stratification']), ')')
  
  
  my_txt <- paste0(formatC(RR_means, format = 'f', 2),
                   ' (+/-) ',
                   formatC(RR_sds, format = 'f', 2))
  
  scattor_p <- ggplot(data = scattor_dat, aes(x, y, color = Stratification))
  
  for (j in 1:length(line_data_list)) {
    
    scattor_p <- scattor_p + annotate(geom = "line",
                                      x = line_data_list[[j]][,'x1'],
                                      y = line_data_list[[j]][,'x2'],
                                      size = 0.5, linetype = 'dotted', colour = '#808080')
    
    pos_1 <- which.min(abs(line_data_list[[j]][,'x1'] - 0.7 * outcome_y_max[i]))
    pos_2 <- which.min(abs(line_data_list[[j]][,'x2'] - 0.7 * outcome_y_max[i]))
    
    pos <- min(pos_1, pos_2)
    
    scattor_p <-  scattor_p + annotate(geom = "text",
                                       x = line_data_list[[j]][pos,'x1'],
                                       y = line_data_list[[j]][pos,'x2'],
                                       label = paste0('RR:\n', formatC(RR_list[j], format = 'f', digits = 2)),
                                       size = 3, fontface = 2, colour = '#303030', hjust = 0.5, angle = 315)
    
  }
  
  scattor_p <- scattor_p + geom_point(aes(x, y, color = Stratification), size = 1)
  scattor_p <- scattor_p + scale_color_manual(values = data_col, name = NULL)
  
  scattor_p <- scattor_p + theme_bw()
  scattor_p <- scattor_p + xlim(c(0, outcome_y_max[i])) + ylim(c(0, outcome_y_max[i]))
  scattor_p <- scattor_p + xlab(x_title) + ylab(y_title) + ggtitle(title_name) 
  scattor_p <- scattor_p + coord_equal()
  
  scattor_p <-  scattor_p + annotate(geom = "text",
                                     x = outcome_y_max[i] * 0.00,
                                     y = outcome_y_max[i] * 1.00,
                                     label = my_txt.0,
                                     size = 4, fontface = 2, colour = '#000000', hjust = 0, angle = 0)    
  
  scattor_p <-  scattor_p + annotate(geom = "text",
                                     x = outcome_y_max[i] * 0.00,
                                     y = outcome_y_max[i] * 0.92,
                                     label = 'Not meeting ISCD indication',
                                     size = 4, fontface = 2, colour = substr(data_col, 1, 7)[1], hjust = 0, angle = 0)    
  
  scattor_p <-  scattor_p + annotate(geom = "text",
                                     x = outcome_y_max[i] * 0.00,
                                     y = outcome_y_max[i] * 0.84,
                                     label = 'Meeting ISCD indication',
                                     size = 4, fontface = 2, colour = substr(data_col, 1, 7)[2], hjust = 0, angle = 0)    
  
  scattor_p <-  scattor_p + annotate(geom = "text",
                                     x = outcome_y_max[i] * 0.80,
                                     y = outcome_y_max[i] * 0.92,
                                     label = my_txt[1],
                                     size = 4, fontface = 1, colour = substr(data_col, 1, 7)[1], hjust = 0.5, angle = 0)    
  
  scattor_p <-  scattor_p + annotate(geom = "text",
                                     x = outcome_y_max[i] * 0.80,
                                     y = outcome_y_max[i] * 0.84,
                                     label = my_txt[2],
                                     size = 4, fontface = 1, colour = substr(data_col, 1, 7)[2], hjust = 0.5, angle = 0)    
  
  scattor_p <-  scattor_p + annotate(geom = "text",
                                     x = outcome_y_max[i] * 0.80,
                                     y = outcome_y_max[i] * 1.00,
                                     label = p_val_txt,
                                     size = 4, fontface = 1, colour = '#000000', hjust = 0.5, angle = 0)    
  
  scattor_p <- scattor_p + theme(plot.title = element_text(color = "#000000", size = 14, face = 1, hjust = 0),
                                 legend.position = "bottom",
                                 axis.title.x = element_text(color = "#000000", size = 12),
                                 axis.title.y = element_text(color = "#000000", size = 12),
                                 axis.text.x = element_text(color = "#000000", angle = 0, hjust = 1, size = 8, face = "bold"),
                                 axis.text.y = element_text(color = "#000000", angle = 0, hjust = 1, size = 8, face = "bold"))
  
  scattor_p_list[[length(scattor_p_list) + 1]] <- scattor_p
  
}

# 5. Merge plotting

if (!dir.exists(dirname(plot_path))) {dir.create(dirname(plot_path), recursive = TRUE)}

final_p <- ggdraw()
final_p <- final_p + draw_plot(forest_p, x = -0.10, y = 0.53, width = 1.13, height = 0.47)
final_p <- final_p + draw_plot(scattor_p_list[[1]], x = 0, y = 0, width = 0.5, height = 0.53)
final_p <- final_p + draw_plot(scattor_p_list[[2]], x = 0.5, y = 0, width = 0.5, height = 0.53)
final_p <- final_p + draw_plot_label(c("A", "B"), c(0.005, 0.005), c(1, 0.55), size = 18, hjust = 0)

png(plot_path, width = 800, height = 800)

print(final_p)

dev.off()