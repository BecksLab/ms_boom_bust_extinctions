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
  group_by(extinction, time_pnt) %>%
  nest() %>%
  glow_up(model = map(data, ~ lmer(value ~ net_type + (1|community), 
                                   data = .x)))

anova_tbl <- mods %>%
  glow_up(anova = map(model,
                      ~ anova(.x, type = 3))) %>%
  vibe_check(extinction, time_pnt, anova) %>%
  unnest(anova)

emm_tbl <- mods %>%
  glow_up(
    emm = map(model,~ emmeans(.x, ~ net_type) %>%
                as.data.frame())) %>%
  vibe_check(extinction, time_pnt, emm) %>%
  unnest(emm)

letters_tbl <- mods %>%
  glow_up(letters = map(model, ~ 
                          {cld(
                            emmeans(.x, "net_type"),
                            Letters = letters) %>%
                              as.data.frame()
                          })) %>%
  vibe_check(extinction, time_pnt, letters) %>%
  unnest(letters)

# --- 4. CLD PALETTE GENERATION ---

plot_tbl_base <- emm_tbl %>%
  left_join(
    letters_tbl %>% 
      vibe_check(extinction, time_pnt, net_type, .group),
    by = c("net_type", "time_pnt", "extinction")
  ) %>%
  mutate(cont_group = trimws(.group)) %>%
  group_by(time_pnt, extinction) %>%
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
  glow_up(ext_scen = factor(paste0(extinction, "_", time_pnt),
                            levels = c("topo_creation", "topo_realised", "dyn_realised")),
          ext_scen_label = assign_ext_scen(ext_scen)) %>%
  slay(ext_scen)

# --- 5. PLOTTING CLD RESULTS (WITH BLENDED GRADIENTS) ---

plot_contrasts <- function(comm) {
  
  plot_tbl %>%
    filter(ext_scen_label == comm) %>%
    ggplot(aes(x = net_type,
               y = emmean,
               colour = blended_color,
               shape = net_group)) +
    geom_point(size = 3,
               show.legend = FALSE) +
    geom_errorbar(
      aes(ymin = asymp.LCL, 
          ymax = asymp.UCL),
      width = 0.15) + 
    geom_text(aes(y = asymp.UCL + 0.01, 
                  label = cont_group), 
              size = 4.5, fontface = "bold", vjust = 0 ) +
    scale_colour_identity() +
    scale_shape_manual(values = net_shapes) +
    ylim(0, 0.5) +
    labs(x = NULL,
         y = "Robustness",
         title = stringr::str_to_title(comm)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_contrast <- plot_tbl %>%
  distinct(ext_scen_label) %>%
  pull(ext_scen_label) %>%
  map(plot_contrasts)

wrap_plots(
  ggplot() +
    labs(title = "Creation: Topological extinctions") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)),
  p_contrast[[1]] +
    theme(plot.title = element_blank()) + plots[[1]] +
    theme(plot.title = element_blank()),
  ggplot() +
    labs(title = "Realised: Topological extinctions") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)),
  p_contrast[[2]] +
    theme(plot.title = element_blank()) + plots[[2]] +
    theme(plot.title = element_blank()),
  ggplot() +
    labs(title = "Realised: Dynamic extinctions") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)),
  p_contrast[[3]] +
    theme(plot.title = element_blank()) + plots[[3]] +
    theme(plot.title = element_blank()),
  ncol = 1,
  heights = c(0.05, 1, 0.05, 1, 0.05, 1)
)

ggsave("../figures/robustnessHighLevel.png",
       width = 5000, 
       height = 6000, 
       units = "px", dpi = 500)

