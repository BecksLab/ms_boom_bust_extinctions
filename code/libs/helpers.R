# function to make cluster plots

make_cluster_plot <- function(hc,
                              title = NULL){
  
  dend <- ggdendro::dendro_data(hc)
  
  ggplot() +
    geom_segment(
      data = dend$segments,
      aes(x = x,
          y = y,
          xend = xend,
          yend = yend),
      linewidth = 0.8) +
    geom_label(
      data = dend$labels %>%
        left_join(col_df),
      aes(x = x,
          y = -0.02 * max(dend$segments$y),
          label = label,
          colour = contrast,
          fill = col),
      size = 4,
      hjust = 0,
      show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_identity() +
    coord_flip(clip = "off") +
    scale_y_reverse(expand = expansion(mult = c(0.05, 0.02))) +
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text( hjust = 0.5, face = "bold", size = 16),
          plot.margin = margin(t = 10, r = 80, b = 10, l = 10))
}