# --- 1. LIBRARIES ---
library(genzplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)


# --- 1. DATA PREPARATION & MODELLING ---

# robustness data
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

# degree of downsampling
topology <-
  read_csv("outputs/paleo_topology.csv") %>%
  yeet(community != "russia") %>%
  drop_na() %>%
  glow_up(stage = factor(stage,
                         levels = c("creation", "burnin")))

meta_connectance <-
  topology %>%
  yeet(net_type == "metaweb",
       stage == "creation") %>%
  transmute(net_id, community,
            co_metaweb = connectance)

delta_co <-
  topology %>%
  yeet(stage == "creation") %>%
  left_join(meta_connectance,
            by = c("net_id", "community")) %>%
  glow_up(delta_co = case_when(net_type %in% c("atn", "niche") ~ NA_real_,
                               net_type == "metaweb" ~ NA_real_,
                               TRUE ~ co_metaweb - connectance),
          net_group = assign_net_group(net_type)) %>%
  vibe_check(net_id, stage, community, net_type, net_group, delta_co)

# ------------------------------------------------------------------------------
# Distance in robustness from reference network states
# ------------------------------------------------------------------------------

topo_df <-
  extinction_df %>%
  yeet(extinction != "dyn") %>%
  vibe_check(net_id, community, scenario, net_type, value, time_pnt)

# reference robustness values
ref_df <-
  topo_df %>%
  yeet(net_type %in% c("metaweb", "atn", "niche")) %>%
  vibe_check(net_id, community, scenario, net_type, value, time_pnt) %>%
  pivot_wider(names_from = net_type,
              values_from = value)


# downsampled networks only
dist_df <-
  topo_df %>%
  filter(str_detect(net_type, "^down")) %>%
  left_join(delta_co) %>%
  left_join(ref_df) %>%
  glow_up(d_meta  = abs(value - metaweb),
          d_atn   = abs(value - atn),
          d_niche = abs(value - niche))

dist_long <-
  dist_df %>%
  pivot_longer(
    starts_with("d_"),
    names_to = "reference",
    values_to = "distance"
  ) %>%
  mutate(
    reference = recode(
      reference,
      d_meta = "metaweb",
      d_atn = "atn",
      d_niche = "niche"
    )
  )


mod_downsamp <-
  dist_long %>%
  filter(time_pnt == "creation") %>%
  group_nest(scenario) %>%
  glow_up(model = map(data, ~lmer(distance ~ delta_co * net_type * reference +
                                    (1|community) + (1|net_id),
                                  data = .x)),
          anova = map(model, ~car::Anova(.x, type = 3) %>%
                        as.data.frame() %>%
                        rownames_to_column("term")),
          singular = map_lgl(model, isSingular),
          emm = map(model, ~emmeans(.x, ~ net_type * reference)),
          trends = map(model, ~emtrends(.x, ~ net_type * reference, 
                                        var = "delta_co")),
          emip = map(model, ~emmip(.x, reference ~ delta_co | net_type)))

mod_downsamp_anova <-
  mod_downsamp %>%
  vibe_check(scenario, anova) %>%
  unnest(anova)

mod_downsamp_emm <-
  mod_downsamp %>%
  vibe_check(scenario, emm) %>%
  mutate(emm = map(emm, as.data.frame)) %>%
  unnest(emm)

mod_burnin <-
  dist_long %>%
  group_nest(scenario) %>%
  glow_up(model = map(data, ~lmer(distance ~ delta_co * net_type * reference * time_pnt +
                                    (1|community) + (1|net_id),
                                  data = .x)),
          singular = map_lgl(model, isSingular),
          anova = map(model, ~car::Anova(.x, type = 3) %>%
                        as.data.frame() %>%
                        rownames_to_column("term")),
          emm = map(model, ~emmeans(.x, ~ net_type * reference * time_pnt)),
          trends = map(model, ~emtrends(.x, ~ net_type * reference * time_pnt,
                                        var = "delta_co")),
          emip = map(model, ~emmip(.x, reference ~ delta_co | net_type)))

mod_burnin_anova <-
  mod_burnin %>%
  vibe_check(scenario, anova) %>%
  unnest(anova)

mod_burnin_emm <-
  mod_burnin %>%
  vibe_check(scenario, emm) %>%
  glow_up(emm = map(emm, as.data.frame)) %>%
  unnest(emm)

mod_burnin_trends <-
  mod_burnin %>%
  vibe_check(scenario, trends) %>%
  mutate(trends = map(trends, as.data.frame)) %>%
  unnest(trends)

burnin_change <-
  mod_burnin_trends %>%
  vibe_check(scenario, net_type, reference, time_pnt, delta_co.trend) %>%
  pivot_wider(names_from = time_pnt,
              values_from = delta_co.trend) %>%
  glow_up(slope_change = realised - creation)

pred_creation <-
  mod_downsamp %>%
  vibe_check(scenario, model) %>%
  glow_up(pred = map(model, ~emmeans(.x, ~ delta_co * net_type * reference,
                                     at = list(delta_co = c(0, 0.25, 0.5))) %>% 
                       as.data.frame())) %>%
  unnest(pred)


p_creation <- 
  ggplot(pred_creation,
         aes(delta_co,
             emmean,
             colour = reference)) +
  geom_line(show.legend = FALSE) +
  facet_grid(scenario~net_type) +
  labs(y = "Predicted robustness distance",
       x = "Degree of downsampling (Δ connectance)")

p_realisation <-
  ggplot(mod_burnin_trends,
         aes(time_pnt,
             delta_co.trend,
             colour = reference)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             colour = minni_silver) +
  geom_point(size=3) +
  geom_line(aes(group = interaction(reference,net_type)))+
  scale_colour_manual(values = net_type_col_pal) +
  facet_grid(scenario~net_type) +
  labs(
    x = "Assembly stage",
    y = "Effect of burn-in of distance in robustness",
    colour = "Reference state"
  )

p_creation + p_realisation +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggplot(dist_long,
       aes(delta_co,
           distance,
           colour = reference,
           linetype = time_pnt)) +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = net_type_col_pal) +
  facet_grid(scenario ~ net_type)

mod_downsamp_trends <-
  mod_downsamp %>%
  transmute(scenario,
            trends = map(model, ~
                           emtrends(.x,
                                    ~ net_type * reference,
                                    var = "delta_co") %>%
                           summary(infer = TRUE) %>%
                           as.data.frame())) %>%
  unnest(trends) %>%
  glow_up(sig = p.value < 0.05,
          alpha = if_else(sig, 1, 0.25))

pred_creation <-
  pred_creation %>%
  left_join(mod_downsamp_trends %>%
              vibe_check(scenario, net_type, reference, alpha),
            by = c("scenario","net_type","reference"))

p_creation <-
  ggplot(pred_creation,
       aes(delta_co,
           emmean,
           colour = reference,
           alpha = alpha,
           group = reference)) +
  geom_hline(yintercept = 0,
             colour = minni_tan) +
  geom_line(linewidth = 1) +
  scale_alpha_identity() +
  facet_grid(scenario ~ net_type)+
  scale_colour_manual(values = net_type_col_pal) +
  labs(y = "Predicted robustness distance",
       x = "Degree of downsampling (Δ connectance)")

mod_burnin_trends <-
  mod_burnin %>%
  transmute(scenario,
            trends = map(model, ~
                           emtrends(.x,
                                    ~ net_type * reference * time_pnt,
                                    var = "delta_co") %>%
                           summary(infer = TRUE) %>%
                           as.data.frame())) %>%
  unnest(trends) %>%
  glow_up(sig = p.value < 0.05,
          alpha = if_else(sig, 1, 0.25))

pred_burnin <-
  mod_burnin %>%
  transmute(scenario,
            pred = map(model, ~ emmeans(.x,
                                        ~ delta_co | net_type * reference * time_pnt,
                                        at = list(delta_co = seq(0, 0.5, length.out = 50))) %>%
                         as.data.frame())) %>%
  unnest(pred) %>%
  left_join(mod_burnin_trends %>%
              vibe_check(scenario, net_type, reference, alpha))


p_realisation <- 
  ggplot(pred_burnin,
       aes(delta_co,
           emmean,
           colour = reference,
           group = interaction(reference, time_pnt, alpha))) +
  geom_hline(yintercept = 0,
             colour = minni_tan) +
  geom_line(aes(alpha = alpha,
                linetype = time_pnt),
            linewidth = 1) +
  scale_alpha_identity() +
  scale_colour_manual(values = net_type_col_pal) +
  facet_grid(scenario ~ net_type)

p_creation + p_realisation +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "../figures/robustnessDownsampling.png",
  width = 9500,
  height = 10000,
  units = "px",
  dpi = 500
)
