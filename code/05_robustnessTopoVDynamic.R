library(tidyverse)
library(broom)
library(emmeans)
library(genzplyr)

source("libs/plotting_themes.R")

# ------------------------------------------------------------------------------
# Data preparation
# ------------------------------------------------------------------------------

extinction_df <- read.csv("outputs/extinction_summary.csv") %>% 
  vibe_check(!starts_with(c("C_", "S_"))) %>% 
  yeet(community != "russia") %>% 
  pivot_longer(cols = -c(net_id, net_type, extinction_time, community), 
               names_to = c("extinction", "time_pnt", "scenario"), 
               names_pattern = "^(topo|dyn)_(creation|realised)_(.+)$") %>% 
  glow_up(extincion_point = factor(paste0(extinction, "_", time_pnt), 
                                   levels = c("topo_creation", 
                                              "topo_realised", 
                                              "dyn_realised")),
          net_group = assign_net_group(net_type))

# ------------------------------------------------------------------------------
# Calculate topology-dynamic difference
# ------------------------------------------------------------------------------

diff_df <-
  extinction_df %>% 
  yeet(time_pnt == "realised") %>% 
  vibe_check(community, net_id, net_type, scenario, extinction, value) %>% 
  pivot_wider(names_from = extinction, 
              values_from = value) %>% 
  glow_up(delta_R50 = dyn - topo)


# ------------------------------------------------------------------------------
# Models and estimated contrasts
# ------------------------------------------------------------------------------


mod_all <-
  lm(delta_R50 ~ scenario * net_type + community,
     data = diff_df)

car::Anova(mod_all, type = 3)
car::Anova(mod_all, type = 3)["scenario:net_type", ]

summary(mod_all)

glance(mod_all)

emm <-
  emmeans(mod_all,
          ~ net_type | scenario)
pairs(emm, adjust = "sidak")

cld(
  emm,
  Letters = letters,
  adjust = "sidak"
)

joint_tests(mod_all)
joint_tests(emm)

mods <-
  diff_df %>%
  group_nest(scenario) %>%
  mutate(
    model  = map(data, ~ lm(delta_R50 ~ community + net_type, data = .x)),
    anova  = map(model, ~ car::Anova(.x, type = 3)),
    glance = map(model, broom::glance),
    emm    = map(model, ~ emmeans(.x, ~ net_type)),
    pairs  = map(emm, ~ pairs(.x, adjust = "sidak")),
    cld    = map(emm, ~ cld(.x,
                            Letters = letters,
                            adjust = "sidak"))
  )


cld_df <-
  mods %>%
  vibe_check(scenario, cld) %>%
  unnest(cld)

anova_df <-
  mods %>%
  vibe_check(scenario, anova) %>%
  glow_up(anova = map(anova, broom::tidy)) %>%
  unnest(anova) %>%
  yeet(term == "net_type") %>%
  transmute(
    scenario,
    df1 = df,
    F = statistic,
    p = p.value
  )

glance_df <-
  mods %>%
  vibe_check(scenario, glance) %>%
  unnest(glance)

emm_df <-
  mods %>%
  vibe_check(scenario, emm) %>%
  glow_up(emm = map(emm, as.data.frame)) %>%
  unnest(emm)

pairs_df <-
  mods %>%
  vibe_check(scenario, pairs) %>%
  glow_up(pairs = map(pairs, as.data.frame)) %>%
  unnest(pairs)

cld_df
anova_df
glance_df
emm_df
pairs_df

anova_table <-
  mods %>%
  transmute(scenario,
            df_res = map_int(model, df.residual),
            anova = map(anova, broom::tidy)) %>%
  unnest(anova) %>%
  yeet(term == "net_type") %>%
  glow_up(F_text = sprintf("$F_{%d,%d}$ = %.1f",df,
                           df_res,
                           statistic),
          p.value = if_else(p.value < 0.05,
                            "<0.05",
                            as.character(p.value)))

write.csv(anova_table,
          "../tables/robustnessTopoVDynamic.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# CLD colours
# ------------------------------------------------------------------------------

all_letters <-
  unique(unlist(strsplit(tolower(cld_df$.group), ""))) %>%
  intersect(letters)

base_palette <-
  setNames(
    unname(col_pal)[seq_along(all_letters)],
    all_letters
  )


cld_df <-
  cld_df %>%
  glow_up(.group = trimws(.group),
          blended_color = get_cld_color(.group,
                                        color_map = base_palette)) %>%
  glow_up(delta_R50 = emmean) %>%
  glow_up(net_group = assign_net_group(net_type)) %>%
  squad_up(scenario) %>%
  arrange(net_group, 
          .by_group = TRUE) %>%
  glow_up(net_type = factor(net_type, levels = net_levs)) %>%
  disband()


# ------------------------------------------------------------------------------
# Figure panels
# ------------------------------------------------------------------------------

p_contrast <-
  ggplot(cld_df,
         aes(x = net_type,
             y = delta_R50,
             fill = blended_color,
             colour = blended_color,
             shape = net_group)) +
  geom_col() +
  geom_errorbar(
    aes(ymin = lower.CL, 
        ymax = upper.CL),
    width = 0.15) +
  geom_hline(
    yintercept = 0,
    colour = minni_black) +
  geom_text(aes(label = .group),
            fontface = "bold",
            vjust = 1,
            nudge_y = -0.025) +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_shape_manual(values = net_shapes,
                     guide = "none") +
  facet_wrap(~scenario,
             ncol = 2) +
  labs(x = "Network type",
       y = expression(Delta*R[50]~"(dynamic - topological)")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p_lollipop <-
  extinction_df %>%
  yeet(time_pnt == "realised") %>%
  group_by(net_type, scenario, extincion_point, net_group) %>%
  glow_up(net_group = assign_net_group(net_type)) %>%
  squad_up(scenario, extincion_point) %>%
  arrange(net_group, 
          .by_group = TRUE) %>%
  glow_up(net_type = factor(net_type, levels = net_levs)) %>%
  disband() %>%
  ggplot(aes(x = net_type, 
             y = value, 
             colour = extincion_point)) +
  geom_boxplot() +
  scale_colour_manual(values = extinction_pal) +
  scale_shape_manual(values = net_shapes,
                     guide = "none") +
  facet_wrap(~scenario,
             ncol = 2) +
  labs(x = "Network type",
       y = expression(R[50]),
       colour = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ------------------------------------------------------------------------------
# Combine and export
# ------------------------------------------------------------------------------

p_final <-
  p_lollipop + p_contrast +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


ggsave(
  "../figures/robustnessTopoVDynamic_allScenario.png",
  p_final,
  width = 10000,
  height = 10000,
  units = "px",
  dpi = 500
)


p_lollipop_collapsed <-
  extinction_df %>%
  yeet(time_pnt == "realised") %>%
  glow_up(net_group = assign_net_group(net_type)) %>%
  squad_up(extinction, community, scenario, net_id, time_pnt) %>%
  glow_up(net_type = factor(net_type, levels = net_levs),
          extinction = paste0(extinction, "_", time_pnt)) %>%
  slay(net_type) %>%
  disband() %>%
  ggplot(aes(x = net_type, 
             y = value, 
             colour = extinction)) +
  geom_boxplot(outliers = FALSE) +
  scale_colour_manual(values = extinction_pal) +
  labs(x = NULL,
       y = expression(Robustness~(R[50])),
       colour = "Extinction scenario") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "../figures/robustnessTopoVDynamic.png",
  p_lollipop_collapsed,
  width = 5000,
  height = 3000,
  units = "px",
  dpi = 500
)


# ------------------------------------------------------------------------------
# Look at topo creation-realised diff
# ------------------------------------------------------------------------------

diff_topo <-
  extinction_df %>% 
  yeet(extinction == "topo") %>% 
  vibe_check(community, net_id, net_type, scenario, time_pnt, value) %>% 
  pivot_wider(names_from = time_pnt, 
              values_from = value) %>%
  squad_up(scenario, community, net_type) %>%
  nest() %>%
  glow_up(model = map(data,~ tryCatch(t.test(Pair(realised, creation) ~ 1, 
                                             data = .x),
                                      error = function(e) NA)))

ttest_tbl <- diff_topo %>%
  glow_up(ttest = map(model, broom::tidy)) %>%
  vibe_check(scenario, community, net_type, ttest) %>%
  unnest(ttest) %>%
  yeet(p.value < 0.05) %>%
  left_join(comm_ord) %>%
  glow_up(net_type = factor(net_type, levels = net_levs),
          print_name = factor(print_name, levels = comm_levs))

extinction_df %>%
  yeet(extinction == "topo") %>%
  glow_up(net_group = assign_net_group(net_type)) %>%
  squad_up(extinction, community, scenario, net_id, time_pnt) %>%
  glow_up(net_type = factor(net_type, levels = net_levs),
          extinction = paste0(extinction, "_", time_pnt)) %>%
  left_join(comm_ord) %>%
  glow_up(print_name = factor(print_name, levels = comm_levs)) %>%
  slay(print_name, net_type) %>%
  disband() %>%
  ggplot(aes(x = net_type, 
             y = value, 
             colour = extinction)) +
  geom_boxplot(outliers = FALSE) +
  geom_text(data = ttest_tbl,
            aes(x = net_type,
                y = 0.52,
                label = "*"),
            colour = minni_tan) +
  scale_colour_manual(values = extinction_pal) +
  facet_grid(cols = vars(print_name),
             rows = vars(scenario)) +
  labs(x = NULL,
       y = expression(Robustness~(R[50])),
       colour = "Extinction scenario") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plotting data -----------------------------------------------------------

ttest_plot_tbl <-
  extinction_df %>%
  yeet(extinction == "topo") %>%
  glow_up(
    net_group = assign_net_group(net_type)
  ) %>%
  squad_up(
    extinction,
    community,
    scenario,
    net_id,
    time_pnt
  ) %>%
  glow_up(
    net_type = factor(net_type, levels = net_levs),
    extinction = paste0(extinction, "_", time_pnt)
  ) %>%
  left_join(comm_ord) %>%
  glow_up(
    print_name = factor(print_name, levels = comm_levs)
  ) %>%
  slay(print_name, net_type) %>%
  disband()


# export one plot per community x scenario -------------------------------

ttest_plot_tbl %>%
  group_split(scenario) %>%
  walk(\(x) {
    
    scen <- x$scenario[[1]]
    
    p <-
      x %>%
      ggplot(aes(x = net_type,
                 y = value,
                 colour = extinction)) +
      geom_boxplot(outliers = FALSE) +
      geom_text(data = ttest_tbl %>%
                  yeet(scenario == scen),
                aes(x = net_type,
                    y = 0.52,
                    label = "*"),
                colour = minni_tan) +
      scale_colour_manual(values = extinction_pal) +
      facet_wrap(vars(print_name)) +
      labs(x = NULL,
           y = expression(Robustness~(R[50])),
           colour = "Extinction scenario",
           title = scen) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    scenario_name <-
      scen %>%
      str_to_lower() %>%
      str_replace_all(" ", "_")
    
    ggsave(paste0("../figures/ttestTopo/ttestTopo_",
                  scenario_name,
                  ".png"),
           p,
           width = 5000,
           height = 3000,
           units = "px",
           dpi = 500
    )
  })

