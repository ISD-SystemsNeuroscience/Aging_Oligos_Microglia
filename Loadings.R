# Library Loadings ====================================================================
library(tidyverse)
library(Seurat)
library(ggplot2)
library(magrittr)
library(reticulate)
library(writexl)
library(RColorBrewer)
library(scales)

# Functions ====================================================================
qc_mad_limits <- function(vc) {
  md <- median(vc)
  abs_diff <- abs(vc - md)
  mad_lower <- md - 4 * median(abs_diff)
  mad_upper <- md + 4 * median(abs_diff)
  mad <- c(mad_lower, mad_upper)
  mad
}



plot_hist <- function(dataset, gene) {
  df <- FetchData(dataset, gene)
  df %>% 
    ggplot(., aes_string(gene)) +
    geom_histogram(bins = 100)
}


barplot_column <- function(df, column1, column2, variable, colors) {
  df %<>% filter(df[[column2]] == variable)
  
  prop_celltype <- table(df[[column1]]) / sum(table(df[[column1]]))
  
  #prop_celltype %<>% round(., 3)
  
  
  df_CellType <- as.data.frame(prop_celltype) %>% set_colnames(c("Stage", "Percentage"))
  
  df_CellType <- df_CellType %>% filter(Percentage != 0.0)
  df_CellType$Percentage %<>% str_remove(., "0+$") %>% as.numeric()
  #names(colors) <- unique(df[[column1]])
  library(scales) # to format labels in percent
  ggplot(df_CellType, aes(x = 2, y = Percentage, fill = Stage)) +
    geom_bar(stat = "identity", color = "white") +
    # coord_polar(theta = "y", start = 0)+
    geom_text(aes(group = Stage, label = percent(Percentage, trim = F)),
              position = position_stack(vjust = 0.5),
              color = "black"
    ) +
    scale_fill_manual("CellType", values = colors, drop = FALSE) +
    theme_void() +
    # xlim(0.5, 2.5) +
    ggtitle(label = paste(column2, " - ", variable, sep = "")) +
    theme(plot.title = element_text(size = 12, hjust = 0.5))
}

barplot_column_wrapper <- function(df, column1, column2, colors) {
  keeper <- list()
  i <- 1
  for (var in unique(df[[column2]])) {
    keeper[[i]] <- barplot_column(df, 
                                  column1, 
                                  column2, 
                                  variable = var,
                                  colors = colors
    )
    i <- i + 1
  }
  return(keeper)
}


barplot_row <- function(df, column1, column2, variable, colors, text.size = 20) {
  prop.df <- table(df[[column1]], df[[column2]]) / rowSums(table(df[[column1]], df[[column2]])) * 100
  
  gg.df <- as.data.frame(prop.df) %>% 
    set_colnames(c(column1, column2, "Percentage"))
  
  gg.df %<>% filter(gg.df[[column2]] == variable)
  
  ggplot(gg.df, aes_string(x = column1, y = 'Percentage', fill = column1)) +
    geom_bar(stat = 'identity', color = 'white') +
    scale_fill_manual(values = colors) +
    theme_classic() +
    theme(text = element_text(size = text.size)) +
    labs(title = paste0(column2, "_", variable))
  
}

barplot_FreqbyMetadata_zeros <- function(Sobj, metacolumn, cluster_var, colors) {
  # cluster_var is the metadata element that is grouped in the x-axis
  # metacolumn is the element that is stacked on the barplots
  
  mdf <- Sobj@meta.data
  tt <- table(mdf[[cluster_var]], mdf[[metacolumn]]) %>% 
    as.data.frame()
  
  mm <- Sobj[[metacolumn]] %>% table() %>% as.list()
  
  tt$total <- unlist(mm[tt$Var2])
  tt$Prop <- tt$Freq/tt$total
  
  dff <- data.frame()
  for (j in unique(tt$Var1)) {
    keeper <- subset(tt, Var1 == j)
    keeper$Prop_norm <- keeper$Prop/sum(keeper$Prop)
    
    dff <- rbind(dff, keeper)
  }
  colnames(dff)[1:2] <- c("cluster", "Metacolumn")
  
  
  ggplot(dff, aes(y=Prop_norm, x=cluster, fill=Metacolumn)) +
    geom_bar(stat='identity', color = 'black') +
    geom_text(aes(label = percent(Prop_norm, accuracy = 0.1)),
              position = position_stack(vjust = 0.5),
              color = 'black'
    ) +
    # geom_text(aes(label = Metacolumn),
    #           position = position_stack(vjust = 0.3),
    #           size = 3,
    #           color = 'black'
    # ) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = c(0,0)) +
    ylab('normalized proportion') +
    theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.text.x = element_text(angle=45, hjust=1,size = 15),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      axis.line.y = element_blank(),
      axis.line.x = element_blank()
    )
}

ave_heatmap_gg <- function(Sobj, column, genes, gene_clust.dict, col_palette, x.size = 4) {
  # Fetch the data
  df_fetch <- FetchData(
    Sobj,
    c(genes, column)
  )
  
  df_fetch %<>% gather(., key = "GeneName", value = "Expression", -column)
  
  df_fetch_mean <- df_fetch %>%
    group_by(GeneName, CellType) %>%
    mutate(AverageExpression = mean(Expression)) %>%
    select(-Expression) %>%
    distinct(GeneName, CellType, .keep_all = TRUE) %>%
    group_by(GeneName) %>%
    mutate_at(vars(matches("AverageExpression")), function(y) y / sum(y))
  
  df_fetch_mean_merged <- merge(df_fetch_mean, gene_clust.dict, by = "GeneName")
  
  df_fetch_mean_merged$GeneName <- factor(df_fetch_mean_merged$GeneName, levels = genes)
  
  # Draw heatmap
  p_df_fetch_mean_merged <- ggplot(
    data = df_fetch_mean_merged,
    mapping = aes(y = CellType, x = GeneName)
  ) +
    geom_tile(aes(fill = AverageExpression)) + # Expression
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = x.size, angle = 90),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank()
    ) +
    scale_fill_gradientn(colours = col_palette(256), limits = c(0, 1)) + # limits -3 to 3 when no normalized average # myPalette or 2
    facet_grid(CellType ~ cluster, scales = "free", space = "free")
  
  p_df_fetch_mean_merged
}


# Colours ====================================================================
col4groups_orj <- c("#66cd00", "#eec900", "#2d7dd2", "#ff0000") # Homeo., Act., WAM1, WAM2
col4groups <- col4groups_orj[c(4, 3, 2, 1)]
isdblue <- "#23395e"
isdorange <- "#f58220"
isdgrey <- "#8698a8"

palette.10x.clusters <- c("#a3e166","#66cd00","#51a400",
                          "#f4de66","#eec900",
                          "#2d7dd2","#e50000")
col4groups_modified <- c("#66cd00","#eec900","#2d7dd2","#e50000")
