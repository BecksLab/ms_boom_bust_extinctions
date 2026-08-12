# --- 1. LIBRARIES ---
library(genzplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)

source("libs/plotting_themes.R")


# --- 3. DATA PREPARATION & MODELLING ---
df_long <- read.csv("outputs/extinction_summary.csv") %>%
  yeet(community != "russia") %>%
  vibe_check(!starts_with(c("C_", "S_"))) %>%
  pivot_longer(
    cols = -c(net_id, net_type, extinction_time, community),
    names_to = c("extinction", "time_pnt", "scenario"),
    names_pattern = "^(topo|dyn)_(creation|realised)_(.+)$"
  ) %>%
  glow_up(
    net_type = factor(net_type),
    scenario = factor(scenario),
    extinction = factor(extinction),
    net_id = factor(net_id),
    community = factor(community),
    time_pnt = factor(time_pnt)
  )

# Fit linear mixed-effects models
mods <- df_long %>%
  group_by(scenario, extinction, community, time_pnt) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(value ~ net_type, data = .x))
  )

anova_tbl <- mods %>%
  mutate(anova = map(model, broom::tidy)) %>%
  vibe_check(scenario, extinction, community, time_pnt, anova) %>%
  unnest(anova)

emm_tbl <- mods %>%
  mutate(emm = map(model, ~ emmeans(.x, "net_type") |> as.data.frame())) %>%
  vibe_check(scenario, extinction, community, time_pnt, emm) %>%
  unnest(emm)

letters_tbl <- mods %>%
  mutate(
    letters = map(model, ~ {
      cld(
        emmeans(.x, "net_type"),
        Letters = letters
      ) %>%
        as.data.frame()
    })
  ) %>%
  vibe_check(scenario, extinction, community, time_pnt, letters) %>%
  unnest(letters)

# --- 4. CLD PALETTE GENERATION ---

plot_tbl_base <- emm_tbl %>%
  left_join(
    letters_tbl %>% 
      vibe_check(scenario, extinction, community, time_pnt, net_type, .group),
    by = c("scenario", "extinction", "net_type", "community", "time_pnt")
  ) %>%
  mutate(cont_group = trimws(.group)) %>%
  group_by(scenario, extinction, community, time_pnt) %>%
  arrange(desc(emmean), 
          .by_group = TRUE) %>%
  mutate(net_type = factor(net_type, levels = net_type)) %>%
  ungroup()

# Extract unique letters present in your CLD results
all_letters <- unique(unlist(strsplit(tolower(plot_tbl_base$cont_group), "")))
all_letters <- sort(all_letters[all_letters %in% letters])

# Assign your named theme colors directly to individual letters (a, b, c, etc.)
base_palette <- setNames(
  unname(col_pal)[seq_along(all_letters)],
  all_letters
)


# Apply color blending algorithm to each cont_group string
plot_tbl <- plot_tbl_base %>%
  glow_up(blended_color = get_cld_color(cont_group, color_map = base_palette),
          ext_scen = paste0(extinction, "_", time_pnt)) %>%
  glow_up(net_type = factor(net_type, levels = net_levs),
          net_group = assign_net_group(net_type)) %>%
  left_join(comm_ord) %>%
  glow_up(print_name = factor(print_name, levels = comm_levs)) %>%
  glow_up(ext_scen = factor(paste0(extinction, "_", time_pnt),
                            levels = c("topo_creation", "topo_realised", "dyn_realised")),
          ext_scen_label = assign_ext_scen(ext_scen)) %>%
  slay(print_name, ext_scen)

# save this object because we need it for clustering
saveRDS(plot_tbl,
        file = "outputs/derived/robustnessEmmTbl.rds")


# ============================================================
# CONTRAST PLOTS — SAME BREAKDOWN SCHEMA AS CLUSTERING
# ============================================================

# 1. Build one contrast plot for each community × extinction
#    scenario

contrast_plots <-
  plot_tbl %>%
  mutate(
    plot = map2(
      print_name,
      ext_scen,
      \(comm, ext) {
        
        x <-
          plot_tbl %>%
          filter(
            print_name == comm,
            ext_scen == ext
          )
        
        ggplot(
          x,
          aes(
            x = net_type,
            y = emmean,
            colour = blended_color,
            shape = net_group
          )
        ) +
          geom_point(
            size = 3,
            show.legend = FALSE
          ) +
          geom_errorbar(
            aes(
              ymin = lower.CL,
              ymax = upper.CL
            ),
            width = 0.15
          ) +
          geom_text(
            aes(
              y = upper.CL + 0.015,
              label = cont_group
            ),
            size = 4.5,
            fontface = "bold",
            vjust = 0
          ) +
          scale_colour_identity() +
          scale_shape_manual(
            values = net_shapes
          ) +
          labs(
            x = "Network type",
            y = "Robustness",
            title = ext
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(
              angle = 45,
              hjust = 1
            )
          )
        
      }
    )
  ) %>%
  # Keep only one copy of each actual plot
  distinct(
    print_name,
    ext_scen,
    .keep_all = TRUE
  )


# ============================================================
# 2. Export one figure per COMMUNITY
# ============================================================

plot_tbl %>%
  periodt(print_name) %>%
  main_character(print_name) %>%
  walk(\(comm) {
    
    p <-
      plot_tbl %>%
      yeet(print_name == comm) %>%
      slay(
        factor(
          ext_scen,
          levels = names(ext_scen_labs)
        )
      ) %>%
      group_by(scenario) %>%
      nest() %>%
      glow_up(
        plot = map2(
          data,
          scenario,
          \(x, nm) {
            
            ggplot(
              x,
              aes(
                x = net_type,
                y = emmean,
                colour = blended_color,
                shape = net_group
              )
            ) +
              geom_point(
                size = 3,
                show.legend = FALSE
              ) +
              geom_errorbar(
                aes(
                  ymin = lower.CL,
                  ymax = upper.CL
                ),
                width = 0.15
              ) +
              geom_text(
                aes(
                  y = upper.CL + 0.015,
                  label = cont_group
                ),
                size = 4.5,
                fontface = "bold",
                vjust = 0
              ) +
              scale_colour_identity() +
              scale_shape_manual(
                values = net_shapes
              ) +
              facet_grid(rows = vars(ext_scen)) +
              labs(
                x = NULL,
                y = "Robustness",
                title = stringr::str_to_title(nm)
              ) +
              coord_cartesian(ylim = c(NA, 0.6)) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
          }
        )
      ) %>%
      main_character(plot) %>%
      wrap_plots(ncol = 2)+
      plot_annotation(title = comm)
    
    
    short_name <-
      comm %>%
      str_extract("^[^,]+") %>%
      str_to_lower() %>%
      str_replace_all(" ", "_")
    
    
    ggsave(
      paste0(
        "../figures/contrasts/robustnessContrasts_",
        short_name,
        ".png"
      ),
      p,
      width = 5000,
      height = 12000,
      units = "px",
      dpi = 500
    )
    
  })

# ============================================================
# 3. Export one figure per EXTINCTION / TIME POINT
# ============================================================


plot_tbl %>%
  periodt(ext_scen_label) %>%
  main_character(ext_scen_label) %>%
  walk(\(t_p) {
    
    p <-
      plot_tbl %>%
      yeet(ext_scen_label == t_p) %>%
      slay(print_name) %>%
      group_by(print_name) %>%
      nest() %>%
      glow_up(
        plot = map(
          data,
          \(x) {
            
            ggplot(
              x,
              aes(
                x = net_type,
                y = emmean,
                colour = blended_color,
                shape = net_group
              )
            ) +
              geom_point(
                size = 3,
                show.legend = FALSE
              ) +
              geom_errorbar(
                aes(
                  ymin = lower.CL,
                  ymax = upper.CL
                ),
                width = 0.15
              ) +
              geom_text(
                aes(
                  y = upper.CL + 0.015,
                  label = cont_group
                ),
                size = 4.5,
                fontface = "bold",
                vjust = 0
              ) +
              scale_colour_identity() +
              scale_shape_manual(values = net_shapes) +
              facet_wrap(vars(scenario),
                              ncol = 10) +
              labs(x = NULL,
                   y = "Robustness",
                   title = unique(print_name)) +
              coord_cartesian(ylim = c(NA, 0.6)) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
          }
        )
      ) %>%
      main_character(plot) %>%
      wrap_plots(ncol = 1) +
      plot_annotation(title = t_p)
    
    
    short_name <-
      plot_tbl %>%
      yeet(ext_scen_label == t_p) %>%
      pull(ext_scen) %>%
      unique() %>%
      paste(collapse = "_") %>%
      str_to_lower() %>%
      str_replace_all(
        "[^a-z0-9_]+",
        "_"
      )
    
    
    ggsave(
      paste0(
        "../figures/contrasts/robustnessContrasts_",
        short_name,
        ".png"
      ),
      p,
      width = 10000,
      height = 9000,
      units = "px",
      dpi = 500
    )
    
  })

# --- 6. PAIRWISE INDISTINGUISHABLE MATRIX ---

# Pairwise comparisons
pairs_tbl <- mods %>%
  mutate(
    pairs = map(model, ~ {
      emmeans(.x, "net_type") %>%
        contrast(method = "pairwise", adjust = "tukey") %>%
        as.data.frame()
    })
  ) %>%
  unnest(pairs) %>%
  separate(contrast,
           into = c("net_type_1", "net_type_2"),
           sep = " - ") %>%
  glow_up(net_type_1 = trimws(net_type_1),
          net_type_2 = trimws(net_type_2),
          indistinguishable = p.value >= 0.05)

indistinguishable_counts <- pairs_tbl %>%
  squad_up(extinction, time_pnt, community, net_type_1, net_type_2) %>%
  no_cap(
    count_indistinguishable = sum(indistinguishable),
    total_scenarios = n(),
    prop_indistinguishable = count_indistinguishable / total_scenarios,
    .groups = "drop"
  ) %>%
  glow_up(ext_scen = paste0(extinction, "_", time_pnt)) %>%
  left_join(comm_ord) %>%
  glow_up(print_name = factor(print_name, levels = comm_levs),
          ext_scen = factor(paste0(extinction, "_", time_pnt),
                            levels = c("topo_creation", "topo_realised", "dyn_realised")),
          ext_scen_label = assign_ext_scen(ext_scen)) %>%
  slay(ext_scen)


plot_indistinguishable_counts <- function(ext_scenario) {
  
  indistinguishable_counts %>%
    filter(ext_scen_label == ext_scenario) %>%
    ggplot(aes(x = net_type_1,
               y = net_type_2,
               fill = prop_indistinguishable)) +
    geom_tile(colour = "white", 
              linewidth = 0.5) +
    geom_text(aes(label = count_indistinguishable),
              colour = "black",
              size = 3.5) +
    facet_wrap(vars(print_name)) +
    scale_fill_gradient(
      high = "#046A38",
      low = "#DDCBA4",
      name = "Proportion\nScenarios\nIndistinguishable") +
    labs(title = stringr::str_to_title(ext_scenario),
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold",hjust = 0.5,size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
}

p_indistinguishable_counts <- indistinguishable_counts %>%
  distinct(ext_scen_label) %>%
  pull(ext_scen_label) %>%
  purrr::map(plot_indistinguishable_counts)

wrap_plots(p_indistinguishable_counts,
           ncol = 1,
           guides = "collect") &
  theme(legend.position = "right")

ggsave("../figures/robustnessIndistinguishable.png",
       width = 9000, 
       height = 7000, 
       units = "px", dpi = 500)

