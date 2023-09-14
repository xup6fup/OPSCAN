
library(magrittr)
library(ggplot2)
library(scales) 
library(cowplot)
library(gridExtra)

# 0. Data path

fill_color <- c('#074F57B2', '#077187B2', '#74A57FB2', 
                '#694F5D', '#EFC7C2', '#BFD3C1',
                '#A7BED3', '#C6E2E9')

plot_path <-  './result/Fig 02.png'
data_path <-  './data/RCT analysis/RCT data.RData'

# 1. Load data

load(data_path)

# 2. Bar plot

data_list <- list()

data_list[[1]] <- follow_data[follow_data[,'Group'] %in% 'Screening',]

names(data_list) <- c('All patients in the screening group')

# 3. Key plot

bar_p_list <- list()

for (m in 1:length(data_list)) {
  
  sub_data <- data_list[[m]]
  
  group_dat.1 <- data.frame(x = c(1, 1, 1),
                            x_name = paste0('AI-CXR high risk (n = ', nrow(sub_data), ')'),
                            group = c('No response', 'DXA by physicians', 'DXA by OPSCAN'),
                            num  = as.numeric(table(sub_data[,'DXA_type'])),
                            prop = as.numeric(prop.table(table(sub_data[,'DXA_type']))),
                            stringsAsFactors = FALSE)
  group_dat.1[,'y_pos'] <- cumsum(group_dat.1[,'prop']) - group_dat.1[,'prop'] / 2
  
  sub_sub_data <- sub_data[sub_data[,'DXA_type'] %in% 2L,]
  
  group_dat.2 <- data.frame(x = c(2, 2, 2),
                            x_name = paste0('DXA by OPSCAN (n = ', nrow(sub_sub_data), ')'),
                            group = c('Normal T-score', 'Osteopenia', 'Osteoporosis'),
                            num  = as.numeric(table(sub_sub_data[,'OP_type'])),
                            prop = as.numeric(prop.table(table(sub_sub_data[,'OP_type']))),
                            stringsAsFactors = FALSE)
  group_dat.2[,'y_pos'] <- cumsum(group_dat.2[,'prop']) - group_dat.2[,'prop'] / 2
  
  sub_sub_data <- sub_data[sub_data[,'DXA_type'] %in% 2L & sub_data[,'OP_type'] %in% 2L,]
  
  group_dat.3 <- data.frame(x = c(3, 3),
                            x_name = paste0('Osteoporosis (n = ', nrow(sub_sub_data), ')'),
                            group = c('No medication', 'Medication usage'),
                            num  = as.numeric(table(sub_sub_data[,'AOM_type'])),
                            prop = as.numeric(prop.table(table(sub_sub_data[,'AOM_type']))),
                            stringsAsFactors = FALSE)
  group_dat.3[,'y_pos'] <- cumsum(group_dat.3[,'prop']) - group_dat.3[,'prop'] / 2
  
  group_dat <- rbind(group_dat.1, group_dat.2, group_dat.3)
  
  group_dat[,'group'] <- factor(group_dat[,'group'], levels = c('DXA by OPSCAN', 'DXA by physicians', 'No response', 'Osteoporosis', 'Osteopenia', 'Normal T-score', 'Medication usage', 'No medication'))
  group_dat[,'col'] <- group_dat[,'group']
  levels(group_dat[,'col']) <- fill_color
  
  group_dat[,'prop_txt'] <- paste0('n = ', group_dat[,'num'], '\n(', formatC(group_dat[,'prop'] * 100, 1, format = 'f'), '%)')
  group_dat[group_dat[,'prop'] <= 0.05,'prop_txt'] <- ''
  
  # Plottings
  
  gg_p <- ggplot(data = group_dat, aes(x = x, y = prop, fill = group))
  gg_p <- gg_p + geom_bar(stat = "identity", width = 0.7)
  gg_p <- gg_p + theme_classic()
  gg_p <- gg_p + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05), breaks = seq(0, 1, by = 0.2), labels = paste0(seq(0, 100, by = 20), '%'))
  gg_p <- gg_p + scale_x_continuous(name = '', breaks = group_dat[,'x'], labels = group_dat[,'x_name'], limits = c(0.5, max(group_dat[,'x']) + 0.5))
  
  gg_p <- gg_p + ggtitle(names(data_list)[m])
  gg_p <- gg_p + ylab("Proportion")
  gg_p <- gg_p + scale_fill_manual(values = fill_color)
  
  gg_p <- gg_p + annotate(geom = "text", x = group_dat[,'x'], y = group_dat[,'y_pos'], label = group_dat[,'prop_txt'], size = 3.5, colour = '#000000')
  
  for (i in 1:2) {
    
    current_pos <- which(group_dat[,'group'] %in% c('DXA by OPSCAN', 'Osteoporosis') & group_dat[,'x'] %in% i)
    end_prop <- 1 - group_dat[current_pos, 'prop']
    polygon_color <- paste0(substr(group_dat[current_pos,'col'], 1, 7), 30)
    
    gg_p <- gg_p + annotate(geom = "polygon", x = c(i + 0.35, i + 0.65, i + 0.65, i + 0.35), y = c(end_prop, 0, 1, 1), colour = '#00000000', fill = polygon_color)
    
  }
  
  gg_p <- gg_p + guides(fill = guide_legend(nrow = 3, byrow = TRUE, title = ""))
  
  gg_p <- gg_p + theme(plot.title = element_text(color = "#000000", size = 15, face = "bold"),
                       panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.title.y = element_text(size = 12, face = "bold"),
                       axis.text.x = element_text(angle = 55, hjust = 1, size = 10, colour = 'black'),
                       axis.text.y = element_text(angle = 0, hjust = 1, size = 8, colour = 'black', face = "bold"))
  
  if (m %in% 1) {
    
    gg_p <- gg_p + theme(legend.position = "bottom")
    
  } else {
    
    gg_p <- gg_p + theme(legend.position = "none")
    
  }
  
  bar_p_list[[m]] <- gg_p
  
}

# 4. Merge plotting

if (!dir.exists(dirname(plot_path))) {dir.create(dirname(plot_path), recursive = TRUE)}

final_p <- ggdraw()
final_p <- final_p + draw_plot(bar_p_list[[1]], x = 0, y = 0, width = 1, height = 1)

png(plot_path, width = 400, height = 560)

print(final_p)

dev.off()
