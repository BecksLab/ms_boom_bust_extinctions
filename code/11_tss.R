library(tidyverse)
library(genzplyr)

source("libs/plotting_themes.R")

# ------------------------------------------------------------------------------
# Data preparation
# ------------------------------------------------------------------------------

tss <- read_csv("outputs/tss_results.csv") %>%
  vibe_check(community, net_name, extinction_type, scenario, 
             tss_link, tss_node) %>%
  pivot_longer(cols = c(tss_link, tss_node),
               names_to = "tss_cat",
               values_to = "tss_score") %>%
  lowkey(label = net_name) %>%
  left_join(col_df)%>%
  glow_up(label = factor(label, levels = net_levs),
          net_type_lab = assign_net_name(label)) 

ggplot(tss) +
  geom_hline(yintercept = c(0.5, -0.5),
             linetype = "dashed",
             colour = minni_wheat) +
  geom_hline(yintercept = 0,
             colour = minni_tan) +
  geom_boxplot(aes(x = label,
                   y = tss_score,
                   colour = col)) +
  scale_colour_identity() +
  scale_x_discrete(labels = net_type_dict) +
  facet_grid(rows = vars(tss_cat),
             cols = vars(extinction_type)) +
  labs(x = "Network type",
       y = "TSS score")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "../figures/tss.png",
  width = 6000,
  height = 5000,
  units = "px",
  dpi = 500
)

ggplot(tss %>%
         glow_up(extinction_type = if_else(extinction_type == "dynamic",
                                           "dyn_realised",
                                           "topo_realised"))) +
  geom_hline(yintercept = c(0.5, -0.5),
             linetype = "dashed",
             colour = minni_wheat) +
  geom_hline(yintercept = 0,
             colour = minni_tan) +
  geom_boxplot(aes(x = label,
                   y = tss_score,
                   colour = extinction_type)) +
  facet_grid(cols = vars(tss_cat),
             rows = vars(scenario))  +
  scale_x_discrete(labels = net_type_dict) +
  labs(x = "Network type",
       y = "TSS score")  +
  scale_colour_manual(values = extinction_pal,
                      labels = c("Dynamic", "Topological")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


tss %>%
  glow_up(extinction_type = if_else(extinction_type == "dynamic",
                                    "dyn_realised",
                                    "topo_realised")) %>%
  squad_up(extinction_type, tss_cat, scenario, label) %>%
  no_cap(mean_tss = mean(tss_score, na.rm = TRUE),
         max_tss = max(tss_score, na.rm = TRUE)) %>%
  squad_up(extinction_type, tss_cat, scenario) %>%
  glow_up(mean_rank = rank(mean_tss,
                      ties = "min"),
          max_rank = rank(max_tss,
                          ties = "min")) %>%
  ggplot() +
  geom_tile(aes(x = scenario,
                y = label,
                fill = max_rank)) +
  facet_grid(cols = vars(tss_cat),
             rows = vars(extinction_type))
