# Make a figure of model predictions with data

source("header.R")

dat <- readRDS("dat_v2.rds")
fit <- readRDS("fit_v3.rds")

# Make data for plotting ----
df <- fit$draws("beta") %>%
  as_draws_df() %>%
  mutate(
    log_asv_abundance_1 = `beta[1,1]`,
    log_sample_count_1 = `beta[2,1]`,
    log_asv_abundance_2 = `beta[1,1]` + `beta[1,2]`,
    log_sample_count_2 = `beta[2,1]` + `beta[2,2]`,
    log_asv_abundance_3 = `beta[1,1]` + `beta[1,3]`,
    log_sample_count_3 = `beta[2,1]` + `beta[2,3]`,
    log_asv_abundance_4 = `beta[1,1]` + `beta[1,2]` + `beta[1,3]` + `beta[1,4]`,
    log_sample_count_4 = `beta[2,1]` + `beta[2,2]` + `beta[2,3]` + `beta[2,4]`,
  ) %>%
  select(starts_with("log")) %>%
  pivot_longer(everything(), values_to = "estimate") %>%
  group_by(name) %>%
  point_interval() %>%
  mutate(
    group = str_extract(name, "[1-4]{1}$"),
    name = str_remove(name, "_[1-4]{1}$")
  ) %>%
  pivot_wider(values_from = c(estimate, .lower, .upper))

# only including core and host-specific ASVs
df_mut <- df[2:3,]

# make a column in the data df with only 2 color options
# aka removing instances of "black" 
dat_no_others = dat[
  (dat$soil_host_core_color == "darkgoldenrod1") | 
    (dat$host_spec_color == "slateblue1"),]

dat_no_others$plot_color = dat_no_others$soil_host_core_color
dat_no_others$plot_color = str_replace(dat_no_others$plot_color, 
                                       "black", "slateblue1")


## Make figure ----
ggplot(df_mut, aes(
  x = estimate_log_sample_count, y = estimate_log_asv_abundance,
  xmin = .lower_log_sample_count, ymin = .lower_log_asv_abundance,
  xmax = .upper_log_sample_count, ymax = .upper_log_asv_abundance,
  color = as.factor(group)
)) +
  facet_grid(group ~ ., labeller = labeller(
    group = c(
      #`1` = "Others",
      `2` = "Core Host + Soil",
      `3` = "Host-specific"
      #`4` = "both"
    ))) +
  geom_point(
    data = dat_no_others,
    mapping = aes(x = log_sample_count, y = log_asv_abundance,
                  color = as.factor(group)), 
    inherit.aes = FALSE, alpha = 0.2
  ) +
  geom_linerange(orientation = "x", size = 2) +
  geom_linerange(orientation = "y", size = 2) +
  geom_point(size = 5, shape = 21, fill = "white") +
  scale_color_manual(values = c("#e65262", "#39acd5")) +
  xlab("log(number_of_samples)") +
  ylab("log(asv_abundance)") +
  theme_cowplot() +
  theme(legend.position = "none")

## Save figure ----
ggsave("Fig_6_abundance_occupancy_scatterplot.png", width = 6, height = 5)
