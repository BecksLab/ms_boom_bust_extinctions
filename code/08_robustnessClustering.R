# clustering

library(tidyverse)
library(ggdendro)
library(patchwork)
library(genzplyr)

source("libs/plotting_themes.R")
source("libs/helpers.R")

# import emm means reaults
plot_tbl2 <- readRDS("outputs/derived/robustnessEmmTbl.rds")

cluster_tbl <-
  plot_tbl2 %>%
  group_by(ext_scen) %>%
  nest() %>%
  glow_up(## Build network type × (community × scenario) matrix
    mat = map(data, \(x){
      
      x %>% unite(condition,
                  community,
                  scenario,
                  remove = TRUE) %>%
        vibe_check(net_type, condition, emmean) %>%
        pivot_wider(names_from = condition,
                    values_from = emmean) %>%
        tibble::column_to_rownames("net_type") %>%
        as.matrix()
    }),
    
    ## Standardise each community/scenario
    X = map(mat, scale),
    ## Distance between network types
    dist = map(X,
               dist,
               method = "euclidean"),
    ## Hierarchical clustering
    hc = map(dist,
             hclust,
             method = "ward.D2"),
    ## Convert to dendrogram
    dend = map(hc,
               as.dendrogram),
    ## Colour branches (choose k)
    dend = map(dend,
               color_branches,
               k = 3),
    ## Nice labels
    dend = map(dend,
               set,
               "labels_cex",
               1.1),
    ext_scen_label = assign_ext_scen(ext_scen)) %>%
  slay(ext_scen)

plots <-
  
  purrr::map2(
    
    cluster_tbl$hc,
    cluster_tbl$ext_scen_label,
    
    \(hc, nm){
      
      dend <- ggdendro::dendro_data(hc)
      
      make_cluster_plot(hc,
                        title = nm)
      
    }
  )

p_clusters <-
  wrap_plots(plots,
             ncol = 1)

ggsave("../figures/robustnessClusters.png",
       p_clusters,
       width = 5000, 
       height = 6000, 
       units = "px", dpi = 500)


# ------------------------------------------------------------
# 2. Build clustering objects
# ------------------------------------------------------------

cluster_tbl_comm <-
  plot_tbl2 %>%
  glow_up(print_name = factor(print_name,
                              levels = comm_levs)) %>%
  group_by(print_name, ext_scen) %>%
  nest() %>%
  glow_up(
    # network type x scenario matrix
    mat = map(data, \(x){
      x %>% 
        unite(condition,
              community,
              scenario,
              remove = TRUE) %>%
        vibe_check(net_type, condition, emmean) %>%
        pivot_wider(names_from = condition,
                    values_from = emmean) %>%
        tibble::column_to_rownames("net_type") %>%
        as.matrix()
    }),
    # standardise
    X = map(mat, scale),
    # distance
    dist = map(X,
               dist,
               method = "euclidean"),
    # clustering
    hc = map(dist,
             hclust,
             method = "ward.D2"),
    # dendrogram
    dend = map(hc,
               as.dendrogram),
    # colour branches
    dend = map(dend,
               color_branches,
               k = 3),
    # labels
    dend = map(dend,
               set,
               "labels_cex",
               1.1),
    ext_scen_label = assign_ext_scen(ext_scen)) %>%
  ungroup()

cluster_plots <-
  cluster_tbl_comm %>%
  mutate(
    plot = map(
      hc,
      make_cluster_plot
    )
  )

# export all community plots

cluster_tbl_comm %>%
  periodt(print_name) %>%
  main_character(print_name) %>%
  walk(\(comm){
    
    p <-
      cluster_plots %>%
      yeet(print_name == comm) %>%
      slay(factor(ext_scen,
                  levels = names(ext_scen_labs))) %>%
      glow_up(plot = map2(plot,
                          print_name,
                          \(p, nm)
                          p + labs(title = nm))) %>%
      main_character(plot) %>%
      wrap_plots(ncol = 1) +
      plot_annotation(title = comm)
    
    short_name <-
      comm %>%
      str_extract("^[^,]+") %>%
      str_to_lower() %>%
      str_replace_all(" ", "_")
    
    ggsave(paste0("../figures/clustering/robustnessClusters_",
                  short_name,
                  ".png"),
           p,
           width = 5000,
           height = 6000,
           units = "px",
           dpi = 500)
    
  })

# write by extinction time point

cluster_tbl_comm %>%
  periodt(ext_scen_label) %>%
  main_character(ext_scen_label) %>%
  walk(\(time_pnt){
    
    p <-
      cluster_plots %>%
      yeet(ext_scen_label == time_pnt) %>%
      slay(print_name) %>%
      glow_up(plot = map2(plot,
                          print_name,
                          \(p, nm)
                          p + labs(title = nm))) %>%
      main_character(plot) %>%
      wrap_plots(ncol = 1) +
      plot_annotation(title = time_pnt)
    
    short_name <-
      cluster_plots %>%
      yeet(ext_scen_label == time_pnt) %>%
      periodt(ext_scen) %>%
      main_character(ext_scen)
    
    ggsave(
      paste0("../figures/clustering/robustnessClusters_",
             short_name,
             ".png"),
      p,
      width = 5000,
      height = 9000,
      units = "px",
      dpi = 500)
    
  })
