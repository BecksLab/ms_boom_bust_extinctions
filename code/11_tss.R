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
