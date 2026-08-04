library(tidyverse)
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
                                              "dyn_realised")))

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

mods <-
  diff_df %>%
  group_nest(scenario) %>%
  glow_up(
    model = map(data, ~ lm(delta_R50 ~ community + net_type, data = .x)),
    emm = map(model, ~ emmeans(.x, ~ net_type)),
    cld = map(emm, ~ cld(.x, Letters = letters, adjust = "sidak"))
  )


cld_df <-
  mods %>%
  vibe_check(scenario, cld) %>%
  unnest(cld)


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
  glow_up(
    .group = trimws(.group),
    blended_color = get_cld_color(.group,
                                  color_map = base_palette)
  ) %>%
  left_join(
    diff_df %>%
      group_by(net_type, scenario) %>%
      summarise(delta_R50 = mean(delta_R50),
                .groups = "drop")
  )


# ------------------------------------------------------------------------------
# Figure panels
# ------------------------------------------------------------------------------

p_contrast <-
  ggplot(cld_df,
         aes(x = net_type,
             y = delta_R50,
             fill = blended_color,
             colour = blended_color)) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    colour = minni_black) +
  geom_text(aes(label = .group),
            fontface = "bold",
            vjust = 1,
            nudge_y = -0.03) +
  scale_fill_identity() +
  scale_colour_identity() +
  facet_wrap(~scenario) +
  labs(x = "Network type",
       y = expression(Delta*R[50]~"(dynamic - topological)")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p_lollipop <-
  extinction_df %>%
  yeet(time_pnt == "realised") %>%
  group_by(net_type, scenario, extincion_point) %>%
  summarise(value = mean(value),
            .groups = "drop") %>%
  ggplot(aes(x = net_type, 
             y = value, 
             colour = extincion_point)) +
  geom_line(aes(group = net_type),
            colour = minni_silver) +
  geom_point(size = 3) +
  scale_colour_manual(values = extinction_pal) +
  facet_wrap(~scenario) +
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
  "../figures/robustnessTopoVDynamic.png",
  p_final,
  width = 9000,
  height = 6000,
  units = "px",
  dpi = 500
)
