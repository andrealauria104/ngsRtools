# Plot functions 

## ggpie with percentages ====

# df$main should contain observations of interest
# df$condition can optionally be used to facet wrap
# labels should be a character vector of the same length as group_by(df, main) or
# group_by(df, condition, main) if facet wrapping

pie_chart <- function(data, main, value = NULL, labels = NULL, condition = NULL, stat_title = NULL) {
  require(tidyverse)
  # convert the data into percentages. group by conditional variable if needed
  df <- dplyr::group_by_(data, .dots = c(condition, main)) %>%
    if(is.null(value)) {
      dplyr::summarize(counts = n())
    } else {
      dplyr::mutate(perc = value / sum(value))  
    }  %>%
    dplyr::mutate(perc = counts / sum(counts)) %>%
    dplyr::arrange(desc(perc)) %>%
    dplyr::mutate(label_pos = sum(perc) - cumsum(perc) + perc / 2,
                  #label_pos = cumsum(perc) - 0.65*perc, 
                  perc_text = paste0(round(perc * 100), "%"))
  # label_pos is a tricky variable to define...the one here will work fine in most, but not all, cases
  # if it doesn't work, try toying with the formula for the label_pos to get the desired result
  
  # reorder the category factor levels to order the legend
  df[[main]] <- factor(df[[main]], levels = unique(df[[main]]))
  
  # if labels haven't been specified, use what's already there
  if (is.null(labels)) labels <- as.character(df[[main]])
  
  p <- ggplot2::ggplot(data = df, aes_string(x = factor(1), y = "perc", fill = main)) +
    
    # make stacked bar chart with black border
    geom_bar(stat = "identity", color = "black", width = 2) +
    
    # add the percents to the interior of the chart
    #geom_text(aes(x = 1.10, y = label_pos, label = perc_text), size = 5) +
    geom_label(aes(x = 1.10, y = label_pos, label = perc_text), fontface = 'bold', color = 'black', size = 5, 
               show.legend = FALSE, inherit.aes = FALSE) +
    
    # convert to polar coordinates
    coord_polar(theta = "y") +
    
    # formatting
    scale_y_continuous(breaks = NULL) +
    scale_fill_discrete(name = "", labels = unique(labels)) +
    theme_grey() +
    theme(panel.grid  = element_blank(), 
          axis.ticks = element_blank(),  
          axis.title = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          strip.text.x = element_text(size = 14, face = "bold"),
          strip.text.y = element_text(size = 14, face = "bold"), 
          strip.text = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14, face = "bold"),
          legend.title.align = 0.5,
          legend.text.align = 0.5,
          legend.direction = 'horizontal', 
          legend.position = 'bottom',
          legend.key = element_rect(size = 5),
          legend.key.size = unit(1.5, 'lines'),
          legend.margin = margin(5,5,5,5),
          legend.box.margin = margin(5,5,5,5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          plot.subtitle = element_text(color = "black", size = 14, hjust = 0.5),
          plot.title = element_text(color = "black", size = 16, face = "bold", hjust = 0.5)) +
    guides(fill = guide_legend(override.aes = list(colour = NA))) + # remove black diagonal line from legend
    scale_fill_brewer(palette = "Dark2") + scale_colour_brewer(palette = "Dark2")
  
  # facet wrap if that's happening
  if (!is.null(condition)) {
    
    # create a dataframe on which chi-square tests will be carried out in case there is "condition" variable present
    df2 <- data %>% dplyr::select(condition, main)
    colnames(df2) <- c("col1", "col2")
    p <- p + facet_wrap(condition, labeller = "label_both") 
    p <- p + labs(subtitle = chi_subtitle(jmv::contTables(df2, rows = 'col1', cols = 'col2', phiCra = TRUE), 
                                          effect = stat_title))
  }
  
  return(p)
  
}

## Venn diagram ----
plot_venn_diagram <- function(vd_list
                              , pal    = NULL
                              , outfig = NULL
                              , plot.h = 3
                              , plot.w = 3
                              , ...)
{
  nc <- length(vd_list)
  if(is.null(pal)) {
    pal <- ggsci::pal_aaas()(nc) 
  } else if(is.character(pal)) {
    pal <- pal[1:nc]
  } else if(is.function(pal)) {
    pal <- pal(nc)
  }
  
  names(vd_list) <- gsub("_"," ", names(vd_list))
  
  vd <- VennDiagram::venn.diagram(
    x = vd_list
    
    # diagram
    , margin = 0.1
    , fill = pal
    , alpha = 0.5
    , col = rep('white',nc)
    , filename = NULL
    # , lwd = 0.5
    , lwd = 2
    , lty = 'blank'
    
    # numbers
    , fontfamily = "sans"
    , cex = .6
    
    # title
    , main.fontfamily = "sans"
    , main.cex = 0.6
    
    # names
    , cat.fontfamily = "sans"
    , cat.cex = .6
    , cat.default.pos = "outer"
    , cat.dist = 0.05
    # , rotation = 1
    , ... )
  
  if( !is.null(outfig) ) {
    message(" -- Save to: ", outfig)
    pdf(file=outfig, paper = "a4", h=unit(plot.h,'cm'), w=unit(plot.w,'cm'))
    grid.draw(vd)
    dev.off()
    
  }
  rmlog <- list.files(pattern = 'VennDiagram.*.log')
  unlink(rmlog)
  grid.draw(vd)
}

## Plot Grid ----
plot_grid <- function(pp, nrow, ncol, top_title = NULL, ...) 
{
  if(!is.null(title)) {
    args <- c(pp, list(nrow = nrow
                       , ncol = ncol
                       , top = textGrob(top_title, gp=gpar(fontsize=8,font=8)), ...))
  } else {
    args <- c(pp, list(nrow = nrow, ncol = ncol, ...))  
  }
  
  do.call(gridExtra::grid.arrange, args)
}

## Plot legend ----
get_legend <- function(a_gplot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a_gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
