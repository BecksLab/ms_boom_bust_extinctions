library(ggplot2)
library(ggrepel)
library(ggtext)
library(patchwork)
library(RColorBrewer) 

# for accents
minni_silver <- "#687982"
minni_black <-  "#292726"
minni_wheat <-  "#F2E4C7"
minni_tan   <-  "#5E4735"

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
        legend.margin = margin(1, 2, 1, 2),
        plot.title = element_text( hjust = 0.5, face = "bold", size = 16)
  )

# this set the theme once so dont need to call every time when plotting
set_theme(figure_theme)

# ---- Colour Palettes ----

# for network types
net_type_col_pal <- c(
  "down_power" = "#59BD96", 
  "metaweb"    = "#5DB6E7", 
  "niche"      = "#EAAA00",
  "down_rand"  = "#216F52",
  "down_link"  = "#44E4A7",
  "down_niche" = "#154734",
  "atn"        = "#A89968")

net_type_col_pal <- c(
  "down_power" = "#28634D",
  "metaweb"    = "#154734",
  "niche"      = "#A6192E",
  "down_rand"  = "#4C8064",
  "down_link"  = "#739675",
  "down_niche" = "#A0A47A",
  "atn"        = "#EAAA00"
)


net_type_col_pal_contrast <- c(
  "down_power" = minni_wheat, 
  "metaweb"    = minni_wheat, 
  "niche"      = minni_wheat,
  "down_rand"  = minni_wheat,
  "down_link"  = minni_wheat,
  "down_niche" = minni_black,
  "atn"        = minni_black)

col_df <-
  cbind(net_type_col_pal,net_type_col_pal_contrast) %>%
  as_tibble() %>%
  lowkey(col = net_type_col_pal,
         contrast = net_type_col_pal_contrast) %>%
  glow_up(label = names(net_type_col_pal_contrast))

# General colour vibes
col_pal <- c(
  "iron_range_red"  = "#A6192E",
  "harvest_gold"    = "#EAAA00",
  "gold_metallic"   = "#A89968",
  "forest_green"    = "#154734",
  "minnesota_wheat" = "#DDCBA4",
  "minnesota_sage"  = "#91A47A"
)

# extinction colours
extinction_pal <- c(
  "topo_realised"  = "#046A38", 
  "dyn_realised"   = "#FFB81C",
  "topo_creation"  = "#739675"
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
# helper function that can be quickly deployed to set the network types
assign_net_group <- function(net_type) {
  
  case_when(str_detect(net_type, "down") ~ "Downsampled web",
            net_type %in% c("niche", "atn") ~ "Realised web",
            net_type %in% c("metaweb") ~ "Metaweb",
            TRUE ~ net_type)
}

# helper function that can be quickly deployed to set the extinction scenario
# as well as network stage
assign_ext_scen <- function(ext_scen) {
  
  case_when(ext_scen == "dyn_realised" ~ "Realised: Dynamic extinctions",
            ext_scen == "topo_realised" ~ "Realised: Topological extinctions",
            ext_scen == "topo_creation" ~ "Creation: Topological extinctions")
  
}

ext_scen_labs = c("Creation: Topological extinctions",
                  "Realised: Topological extinctions",
                  "Realised: Dynamic extinctions")

# ---- Network ordering ----

net_levs <- c("metaweb",
              "down_rand",
              "down_power",
              "down_link",
              "down_niche",
              "niche",
              "atn")

# ---- Community ordering ----

comm_ord <- 
  topology %>% 
  yeet(stage == "creation", 
       net_type == "metaweb") %>% 
  vibe_check(community, connectance) %>% 
  distinct() %>% 
  slay(-connectance) %>%
  glow_up(print_name = paste0(stringr::str_to_title(community),
                              ", Co = ",
                              round(connectance, digits = 3)))

comm_levs <- c(comm_ord$print_name)
