library(tidyverse)
library(genzplyr)
library(patchwork)
library(ggridges)

source("libs/plotting_themes.R")

trait_df <- read.csv("outputs/S1_bodysize/species_metadata.csv") %>%
  glow_up(stage = if_else(stage == "post_burn_in",
                          "burnin",
                          stage)) %>%
  left_join(read.csv("outputs/S1_bodysize/spp_degrees.csv") %>%
              lowkey(species_id = spp_id)) %>%
  distinct() %>%
  glow_up(bm_specification = if_else(str_detect(net_type, "_empiricalbm"),
                                     "Empirical",
                                     "Internal"),
          net_type = str_remove(net_type, "_empiricalbm"))

survivors <- trait_df %>%
  vibe_check(net_id, stage, net_type, original_id, community, bm_specification) %>%
  yeet(stage == "burnin") %>%
  glow_up(survived = "yes") %>%
  vibe_check(-stage) %>%
  yeet(community != "russia")

trait_df %>%
  yeet(stage == "creation") %>%
  yeet(net_type != "niche") %>%
  left_join(survivors) %>%
  glow_up(survived = factor(if_else(is.na(survived), "no", "yes")),
          stage = factor(stage,
                         levels = c("creation", "burnin"))) %>%
  pivot_longer(-c(net_id, net_type, stage, species_id,
                  original_id, trophic_class, survived,
                  t_val, community, bm_specification, can, bodysize)) %>%
  ggplot() +
  geom_density_ridges(aes(x = value,
                          y = net_type,
                          fill = survived),
                      alpha = 0.7) +
  facet_grid(
    cols = vars(name),
    rows = vars(bm_specification),
    scales = "free_x") +
  scale_fill_manual(values = c("yes" = "#046A38",
                               "no" = "#EAAA00")) +
  labs(y = NULL,
       x = "Value",
       fill = "Survived")

ggsave("../figures/bodysize/distributions.png",
       width = 10000,
       height = 3500,
       units = "px",
       dpi = 500)


trait_df %>%
  left_join(survivors) %>%
  glow_up(survived = factor(if_else(is.na(survived), "no", "yes")),
          stage = factor(stage, levels = c("creation", "burnin"))) %>%
  yeet(stage == "creation") %>%
  squad_up(net_type, survived, community, bm_specification) %>%
  tally() %>%
  ggplot(aes(x = net_type, 
             y = n, 
             fill = survived)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("yes" = "#046A38",
                               "no" = "#EAAA00")) +
  facet_wrap(vars(bm_specification)) +
  labs(y = NULL,
       x = NULL,
       title = "Survivorship") +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))

ggsave("../figures/bodysize/survivorship.png",
       width = 5000,
       height = 3500,
       units = "px",
       dpi = 500)

trait_df %>%
  left_join(survivors) %>%
  glow_up(survived = factor(if_else(is.na(survived), "no", "yes")),
          stage = factor(stage, levels = c("creation", "burnin")))  %>%
  ggplot() +
  geom_boxplot(aes(x = net_type,
                   y = log10(bodysize),
                   colour = trophic_class),
               outliers = FALSE) +
  facet_grid(cols = vars(bm_specification)) +
  scale_colour_manual(values = c("basal" = "#046A38",
                                 "top" = "#A6192E",
                                 "intermediate" = "#EAAA00",
                                 "isolated" = "#F2E4C7")) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
  labs(y = "Log10 Body size",
       x = NULL,
       title = "Body size distribution across trophic classes")

trait_df %>%
  left_join(survivors) %>%
  yeet(stage == "creation") %>%
  glow_up(survived = factor(if_else(is.na(survived), "no", "yes"))) %>%
  ggplot() +
  geom_point(aes(y = bodysize,
                 x = trophic_level, 
                 colour = survived),
             alpha = 0.6) +
  scale_colour_manual(values = c("yes" = "#046A38",
                                 "no" = "#EAAA00")) +
  facet_grid(cols = vars(bm_specification),
             rows = vars(net_type))

trait_df %>%
  left_join(survivors) %>%
  glow_up(survived = factor(if_else(is.na(survived), "no", "yes")),
          stage = factor(stage, levels = c("creation", "burnin")))  %>%
  ggplot() +
  geom_boxplot(aes(x = net_type,
                   y = metabolism,
                   colour = survived),
               outliers = FALSE) +
  facet_grid(cols = vars(bm_specification)) +
  scale_colour_manual(values = c("yes" = "#046A38",
                                 "no" = "#EAAA00")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Metabolism",
       x = NULL,
       title = "Metabolism and survivorship")


trait_df %>%
  left_join(survivors) %>%
  yeet(stage == "creation") %>%
  glow_up(survived = factor(if_else(is.na(survived), "no", "yes"))) %>%
  ggplot() +
  geom_boxplot(aes(y = bodysize,
                   x = net_type, 
                   colour = survived),
               outliers = FALSE) +
  scale_colour_manual(values = c("yes" = "#046A38",
                                 "no" = "#EAAA00")) +
  facet_grid(cols = vars(bm_specification))


trait_df %>%
  left_join(survivors) %>%
  glow_up(survived = factor(if_else(is.na(survived), "no", "yes")),
          bm_expected = 10^(trophic_level - 1),
          log10_bm = log10(bodysize),
          log10_bm_expected = trophic_level - 1,
          bm_residual = log10_bm - log10_bm_expected) %>%
  ggplot(aes(x = trophic_level, 
             y = log10(bodysize))) +
  geom_point(aes(shape = trophic_class,
                 colour = survived),
             alpha = 0.6) +
  geom_abline(intercept = -1,
              slope = 1,
              linetype = "dashed",
              colour = minni_tan) +
  facet_grid(cols = vars(bm_specification),
             rows = vars(net_type)) +
  scale_colour_manual(values = c("yes" = "#046A38",
                                 "no" = "#EAAA00")) +
  labs(x = "Trophic level",
       y = "log10(body size)",
       title = "Empirical body size vs Z = 10 expectation")

ggsave("../figures/bodysize/expectedBM.png",
       width = 6000,
       height = 5500,
       units = "px",
       dpi = 500)

trait_df %>%
  left_join(survivors) %>%
  yeet(stage == "creation") %>%
  glow_up(survived = factor(if_else(is.na(survived), "no", "yes"))) %>%
  ggplot() +
  geom_density_ridges(aes(x = log10(bodysize),
                          y = net_type,
                          fill = survived),
                      alpha = 0.7) +
  scale_fill_manual(values = c("yes" = "#046A38",
                               "no" = "#EAAA00")) +
  facet_grid(cols = vars(bm_specification))

trait_df %>%
  left_join(survivors) %>%
  yeet(is.na(survived)) %>%
  squad_up(net_id, net_type, original_id) %>%
  no_cap(
    presence = case_when(any(bm_specification == "Internal") &
                           any(bm_specification == "Empirical") ~ "Both",
                         any(bm_specification == "Internal") ~ "Internal",
                         any(bm_specification == "Empirical") ~ "Empirical"),
    .groups = "drop") %>%
  count(net_id, net_type, presence)

trait_df %>%
  left_join(survivors) %>%
  group_by(net_id, net_type, original_id, bm_specification) %>%
  summarise(
    survived = any(survived == "yes"),
    .groups = "drop"
  ) %>%
  group_by(net_type, original_id, bm_specification) %>%
  summarise(
    prop_survived = mean(survived),
    .groups = "drop"
  ) %>%
  group_by(net_type, original_id) %>%
  summarise(
    presence = case_when(
      all(prop_survived == 0) ~ "Always extinct",
      all(prop_survived == 1) ~ "Always survives",
      n_distinct(prop_survived) > 1 ~ "BM specification dependent",
      TRUE ~ "Intermediate"
    ),
    .groups = "drop"
  ) %>%
  squad_up(net_type, presence) %>%
  tally()


status <- trait_df %>%
  yeet(net_type != "niche") %>%
  vibe_check(net_id, net_type, original_id) %>%
  distinct() %>%
  crossing(
    bm_specification = c("Internal", "Empirical")
  ) %>%
  left_join(
    survivors %>%
      vibe_check(net_id, net_type, original_id, bm_specification) %>%
      distinct() %>%
      mutate(survived = TRUE)
  ) %>%
  mutate(
    survived = replace_na(survived, FALSE),
    extinct = !survived
  )

comparison <- status %>%
  vibe_check(net_id, net_type, original_id, bm_specification, extinct) %>%
  pivot_wider(
    names_from = bm_specification,
    values_from = extinct
  ) %>%
  mutate(
    extinction_pattern = case_when(
      Internal == FALSE & Empirical == FALSE ~ "Survives both",
      Internal == TRUE  & Empirical == TRUE  ~ "Extinct both",
      Internal == TRUE  & Empirical == FALSE ~ "Extinct Internal only",
      Internal == FALSE & Empirical == TRUE  ~ "Extinct Empirical only"
    )
  )

subset_test <- comparison %>%
  group_by(net_id, net_type) %>%
  summarise(
    internal_extinct = sum(Internal),
    internal_also_empirical = sum(Internal & Empirical),
    proportion_internal_in_empirical =
      internal_also_empirical / internal_extinct,
    .groups = "drop"
  )

ggplot(
  subset_test,
  aes(x = net_type, 
      y = proportion_internal_in_empirical)) +
  geom_boxplot(colour = "#EAAA00") +
  geom_jitter(width = 0.1, 
              alpha = 0.4,
              colour = "#EAAA00") +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent) +
  labs(x = NULL,
       y = "Proportion of Internal extinctions\nalso extinct under Empirical") +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))

ggsave("../figures/bodysize/propOverlap.png",
       width = 5000,
       height = 3000,
       units = "px",
       dpi = 500)

trait_df %>%
  yeet(stage == "creation") %>%
  yeet(net_type != "niche") %>%
  left_join(survivors) %>%
  glow_up(survived = factor(if_else(is.na(survived), "no", "yes")),
          stage = factor(stage,
                         levels = c("creation", "burnin"))) %>%
  vibe_check(net_id, net_type, stage, species_id, original_id, survived,
             bm_specification, S4_consumer, S5_consumer, S4_resource,
             S5_resource) %>%
  pivot_longer(-c(net_id, net_type, stage, species_id, original_id, survived,
                  bm_specification)) %>%
  ggplot() +
  geom_boxplot(aes(x = net_type,
                   y = value,
                   colour = survived))  +
  facet_grid(cols = vars(name),
             rows = vars(bm_specification))

mod_dat <-
  trait_df %>%
  left_join(survivors) %>%
  glow_up(survived = if_else(is.na(survived), 0, 1),
          bm_expected = 10^(trophic_level - 1),
          log10_bm = log10(bodysize),
          log10_bm_expected = trophic_level - 1,
          bm_residual = log10_bm - log10_bm_expected) %>%
  glow_up(metabolism_z = as.numeric(scale(metabolism)),
          bm_deviation_z = as.numeric(scale(bm_residual)),
          S4_consumer_z = as.numeric(scale(S4_consumer)),
          S4_resource_z = as.numeric(scale(S4_resource)))

library(lme4)

m1 <- glmer(survived ~
              metabolism_z +
              S4_consumer_z +
              S4_resource_z +
              net_type +
              bm_specification +
              (1 | net_id),
            data = mod_dat,
            family = binomial)

m_emp <- glmer(survived ~
                 metabolism_z +
                 bm_deviation_z +
                 S4_consumer_z +
                 S4_resource_z +
                 net_type +
                 (1 | net_id),
               data = mod_dat %>%
                 filter(bm_specification == "Empirical"),
               family = binomial)

m_emp_int <- glmer(survived ~
                     metabolism_z +
                     bm_deviation_z +
                     S4_consumer_z * bm_deviation_z +
                     S4_resource_z * bm_deviation_z +
                     net_type +
                     (1 | net_id),
                   data = mod_dat %>%
                     filter(bm_specification == "Empirical"),
                   family = binomial)

anova(m_emp, m_emp_int, test = "Chisq")

aic_tab <- AIC(m1, m_emp, m_emp_int)

aic_tab %>%
  glow_up(delta_AIC = AIC - min(AIC),
          weight = exp(-0.5 * delta_AIC) /
            sum(exp(-0.5 * delta_AIC))) %>%
  slay(AIC)
