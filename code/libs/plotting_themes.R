library(ggplot2)
library(ggrepel)
library(ggtext)
library(patchwork)
library(RColorBrewer) 

# for accents
minni_silver <- "#687982"
minni_black <- "#1F2C33"

# ---- General Figure Theme ----

figure_theme = 
  theme_classic(16) +
  theme(panel.border = element_rect(colour = minni_black,
                                    fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = colorspace::lighten(minni_black, 0.7),
                                  linewidth = 0.3),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_blank(),
        text = element_text(color = minni_black),
        plot.margin = margin(10, 5, 5, 10),
        legend.margin = margin(1, 2, 1, 2)
  )

# this set the theme once so dont need to call every time when plotting
set_theme(figure_theme)

# ---- Colour Palettes ----

# for network types
net_type_col_pal <- c(
  "down_power" = "#A6192E", 
  "metaweb"    = "#154734", 
  #"niche".    = "#4A8770",
  "niche"      = "#C5E86C",
  "down_rand"  = "#662932",
  "down_link"  = "#332124",
  "down_niche" = "#9FC7B8",
  "atn"        = "#BFFFE7")

# General colour vibes
col_pal <- c(
  "iron_range_red"  = "#A6192E",
  "harvest_gold"    = "#EAAA00",
  "gold_metallic"   = "#A89968",
  "forest_green"    = "#154734",
  "minnesota_wheat" = "#DDCBA4",
  "light_green"     = "#C5E86C"
)

# extinction colours
extinction_pal <- c(
  "topo_realised"  = "#046A38", 
  "dyn_realised"   = "#FFB81C",
  "topo_creation"  = "#63d4a9"
)


# --- Colour blending utilities ----

# Helper function to blend RGB values for letter strings (e.g. "ab" -> mean of 'a' and 'b')

blend_cld_colors <- function(letter_string, color_map) {
  chars <- strsplit(tolower(trimws(letter_string)), "")[[1]]
  chars <- chars[chars %in% names(color_map)]
  
  if (length(chars) == 0) return("#CCCCCC") # Fallback grey
  
  rgb_matrix <- col2rgb(color_map[chars])
  avg_rgb <- rowMeans(rgb_matrix)
  
  rgb(avg_rgb[1], avg_rgb[2], avg_rgb[3], maxColorValue = 255)
}

get_cld_color <- Vectorize(blend_cld_colors, vectorize.args = "letter_string")

# ---- Shapes for network types ----

# we can specify chapes for metawebs, downsampled, and realised webs to help
# with a consistent visual language so to say

net_shapes <- c(
  "Downsampled web" = 15, 
  "Realised web"     = 16, 
  "Metaweb"    = 17)

# ---- Assign network types ----
# helper function that can be quickly deployed to set the netowrk types
assign_net_group <- function(net_type) {
  
  case_when(str_detect(net_type, "down") ~ "Downsampled web",
            net_type %in% c("niche", "atn") ~ "Realised web",
            net_type %in% c("metaweb") ~ "Metaweb",
            TRUE ~ net_type)
}

