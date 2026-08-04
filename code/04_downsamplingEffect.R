library(tidyverse)
library(FactoMineR)
library(vegan)
library(lme4)
library(emmeans)
library(multcomp)

source("libs/plotting_themes.R")

# ------------------------------------------------------------------------------
# Data preparation
# ------------------------------------------------------------------------------

topology <-
  read_csv("outputs/paleo_topology.csv") %>%
  yeet(community != "russia") %>%
  drop_na() %>%
  glow_up(stage = factor(stage,
                         levels = c("creation", "burnin")))

topo_vars <-
  topology %>%
  vibe_check(-net_id, -community, -net_type, -stage) %>%
  names()

# ------------------------------------------------------------------------------
# Connectance difference from the ancestral (metaweb) network
# ------------------------------------------------------------------------------

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
                               TRUE ~ co_metaweb - connectance))

# ------------------------------------------------------------------------------
# PCA of topology
# (ATN and niche projected as supplementary individuals)
# ------------------------------------------------------------------------------

meta <-
  topology %>%
  vibe_check(net_id, community, net_type, stage)

X <-
  topology %>%
  vibe_check(all_of(topo_vars))

supp_rows <-
  which(meta$net_type %in% c("niche", "atn"))

pca <-
  PCA(
    X,
    scale.unit = TRUE,
    ind.sup = supp_rows,
    graph = FALSE
  )

scores <-
  bind_rows(
    bind_cols(meta[-supp_rows, ],
              as.data.frame(pca$ind$coord)),
    bind_cols(meta[supp_rows, ],
              as.data.frame(pca$ind.sup$coord)))

# ------------------------------------------------------------------------------
# Community centroids in PCA space
# ------------------------------------------------------------------------------

centroids <-
  scores %>%
  squad_up(community, stage, net_type) %>%
  no_cap(across(starts_with("Dim"), mean),
         .groups = "drop")

# ------------------------------------------------------------------------------
# Construct trajectories:
# metaweb (creation) -> creation -> burnin
# ------------------------------------------------------------------------------

meta_start <-
  centroids %>%
  yeet(net_type == "metaweb",
       stage == "creation") %>%
  vibe_check(community, Dim.1, Dim.2) %>%
  crossing( net_type = unique(
    centroids$net_type[!centroids$net_type %in%
                         c("metaweb", "niche", "atn")])) %>%
  glow_up(time_stage = "start")

traj <-
  centroids %>%
  glow_up(time_stage = as.character(stage)) %>%
  bind_rows(meta_start) %>%
  glow_up(time_stage = factor(time_stage,
                              levels = c("start", "creation", "burnin"))) %>%
  arrange(community, net_type, time_stage)

# ------------------------------------------------------------------------------
# PCA plot
# ------------------------------------------------------------------------------

ggplot(scores,
       aes(Dim.1,
           Dim.2,
           colour = net_type)) +
  geom_vline(xintercept = 0,
             colour = minni_silver) +
  geom_hline(yintercept = 0,
             colour = minni_silver) +
  stat_ellipse(
    data = yeet(scores,
                !net_type %in% c("niche", "atn"))) +
  geom_point(aes(shape = stage),
             alpha = 0.4) +
  scale_colour_manual(values = net_type_col_pal)

# ------------------------------------------------------------------------------
# Trajectories through PCA space
# ------------------------------------------------------------------------------

p_traj <-
  ggplot(traj,
         aes(Dim.1,
             Dim.2,
             colour = net_type,
             group = interaction(community, net_type))) +
  geom_vline(xintercept = 0,
             colour = minni_silver) +
  geom_hline(yintercept = 0,
             colour = minni_silver) +
  geom_path(arrow = arrow(length = unit(0.25, "cm")),
            linewidth = 1.2,
            alpha = 0.8) +
  facet_wrap(~community) +
  scale_colour_manual(values = net_type_col_pal)

p_traj_all <- 
  ggplot(traj %>%
           yeet(stage %in% c("creation", "burnin")),
         aes(Dim.1,
             Dim.2,
             colour = net_type,
             group = interaction(community, net_type))) +
  geom_vline(xintercept = 0,
             colour = minni_silver) +
  geom_hline(yintercept = 0,
             colour = minni_silver) +
  geom_path(arrow = arrow(length = unit(0.25, "cm")),
            linewidth = 1.2,
            alpha = 0.8) +
  labs(title = "Change in structure during burnin") +
  scale_colour_manual(values = net_type_col_pal)

# ------------------------------------------------------------------------------
# Trajectory test
# ------------------------------------------------------------------------------

burnin_traj <- scores %>%
  vibe_check(net_id, community, net_type, stage, Dim.1, Dim.2, Dim.3) %>%
  pivot_wider(names_from = stage,
              values_from = c(Dim.1, Dim.2, Dim.3)) %>%
  glow_up(displacement = sqrt((Dim.1_burnin - Dim.1_creation)^2 +
                                (Dim.2_burnin - Dim.2_creation)^2 +
                                (Dim.3_burnin - Dim.3_creation)^2))

m_disp <- lmer(log(displacement + 0.01) ~ net_type + (1|community),
               data = burnin_traj)

anova(m_disp)

disp_emm <- emmeans(m_disp, ~ net_type)

cld_disp <- cld(disp_emm,
                Letters = letters,
                adjust = "sidak")

# compact letter display
cld_tble <-
  cld_disp %>%
  as_tibble() %>%
  left_join(as_tibble(disp_emm ))

# Extract unique letters present in your CLD results
all_letters <- unique(unlist(strsplit(tolower(cld_tble$.group), "")))
all_letters <- sort(all_letters[all_letters %in% letters])

# Assign your named theme colors directly to individual letters (a, b, c, etc.)
base_palette <- setNames(
  unname(col_pal)[seq_along(all_letters)],
  all_letters
)

# Apply color blending algorithm to each cont_group string
cld_tble <- cld_tble %>%
  glow_up(blended_color = get_cld_color(.group, color_map = base_palette),
          .group = trimws(.group)) %>%
  arrange(desc(.group), 
          .by_group = TRUE) %>%
  mutate(net_type = factor(net_type, levels = net_type))

p_pca_contrast <- 
  ggplot(cld_tble,
         aes(y = net_type,
             x = emmean,
             colour = blended_color)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(xmin = lower.CL, 
        xmax = upper.CL),
    width = 0.15) + 
  geom_text(aes(label = .group), 
            size = 4.5, 
            fontface = "bold", 
            vjust = 0,
            nudge_y = 0.2) +
  scale_colour_identity() +
  labs(y = "Network type",
       x = "Displacement",
       title = "Differences in displacement during burn in")

# ------------------------------------------------------------------------------
# Multivariate tests (creation networks only)
# ------------------------------------------------------------------------------

topology_creation <-
  topology %>%
  yeet(stage == "creation",
       net_type != "metaweb")

topology_scaled <-
  topology_creation %>%
  vibe_check(all_of(topo_vars)) %>%
  scale()

# Partial RDA

rda_dat <-
  delta_co %>%
  yeet(!is.na(delta_co))

topology_scaled <-
  rda_dat %>%
  vibe_check(all_of(topo_vars)) %>%
  scale()

rda_model <-
  rda(topology_scaled ~
        net_type * delta_co +
        Condition(community),
      data = rda_dat)

# Overall test
anova(rda_model)

# Marginal effects
anova(rda_model, by = "term")

anova(rda_model, by = "margin")

# Canonical axes
anova(rda_model, by = "axis")


rda_sites <-
  scores(rda_model,
         display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("obs") %>%
  bind_cols(rda_dat)

# biplot arrows
rda_bp <-
  scores(rda_model,
         display = "bp") %>%
  as.data.frame() %>%
  rownames_to_column("variable")

# factor centroids
rda_cn <-
  scores(rda_model,
         display = "cn") %>%
  as.data.frame() %>%
  rownames_to_column("factor")

p_rda <-
  ggplot(rda_sites,
         aes(RDA1,
             RDA2,
             colour = net_type)) +
  geom_hline(yintercept = 0,
             colour = minni_silver) +
  geom_vline(xintercept = 0,
             colour = minni_silver) +
  stat_ellipse(level = .68,
               linewidth = 1,
               show.legend = FALSE) +
  geom_point(size = 2,
             alpha = .6,
             show.legend = FALSE) +
  geom_segment(data = rda_bp,
               aes(x = 0,
                   y = 0,
                   xend = RDA1,
                   yend = RDA2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(.25, "cm")),
               show.legend = FALSE) +
  geom_label_repel(data = rda_bp,
                   aes(RDA1,
                       RDA2,
                       label = variable),
                   inherit.aes = FALSE) +
  labs(title = "Effect of Downsampling") +
  scale_colour_manual(values = net_type_col_pal)


rda_scores <- scores(rda_model, display = "sites", choices = 1:3)

rda_dat$rda1 <- rda_scores[,1]
rda_dat$rda2 <- rda_scores[,2]
rda_dat$rda3 <- rda_scores[,3]

m <- lmer(rda1 ~ net_type * delta_co + (1|community),
          data = rda_dat)

anova(m)
# estimate slopes
trends <- emtrends(m, ~ net_type, var = "delta_co")

# compact letter display
cld_tble <-
  cld(trends, Letters = letters, adjust = "sidak") %>%
  as_tibble() %>%
  left_join(as_tibble(trends))

# Extract unique letters present in your CLD results
all_letters <- unique(unlist(strsplit(tolower(cld_tble$.group), "")))
all_letters <- sort(all_letters[all_letters %in% letters])

# Assign your named theme colors directly to individual letters (a, b, c, etc.)
base_palette <- setNames(
  unname(col_pal)[seq_along(all_letters)],
  all_letters
)

# Apply color blending algorithm to each cont_group string
cld_tble <- cld_tble %>%
  glow_up(blended_color = get_cld_color(.group, color_map = base_palette),
          .group = trimws(.group))

p_contrast <- 
  ggplot(cld_tble,
         aes(y = net_type,
             x = delta_co.trend,
             colour = blended_color)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(xmin = lower.CL, 
        xmax = upper.CL),
    width = 0.15) + 
  geom_text(aes(x = delta_co.trend, 
                y = net_type,
                label = .group), 
            size = 4.5, 
            fontface = "bold", 
            vjust = 0,
            nudge_y = 0.15) +
  scale_colour_identity() +
  labs(y = "Network type",
       x = "delta_co.trend",
       title = "Difference in downsample approach and connectance change")


design <- "
  12
  34
"

p_rda + p_contrast +
  p_traj_all + p_pca_contrast +
  plot_annotation(tag_levels = 'A') +
  plot_layout(design = design,
              guides = "collect") &
  theme(legend.position='bottom')

ggsave("../figures/downsamplingEffect.png",
       width = 6000, 
       height = 6000, 
       units = "px", dpi = 500)
