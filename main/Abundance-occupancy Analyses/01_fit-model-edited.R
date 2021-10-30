# This script fits model.stan to ASV abundance-occupancy data

## Source header ----
source("header.R")

## Read data and edit ----

dat <- read.csv("data_for_model.csv") %>%
  mutate(
    log_asv_abundance = log(asv_abundance),
    log_sample_count = log(sample_count),
    group = case_when(
      soil_host_core_color == "black" & host_spec_color == "black" ~ 1,
      soil_host_core_color == "darkgoldenrod1" & host_spec_color == "black" ~ 2,
      soil_host_core_color == "black" & host_spec_color == "slateblue1" ~ 3,
      soil_host_core_color == "darkgoldenrod1" & host_spec_color == "slateblue1" ~ 4,
    )
  )

## Save processed data ----
write_rds(dat, "dat_v2.rds")

## Convert data to list for Stan ----

dat_list <- list(
  x = model.matrix(lm(log_asv_abundance ~ soil_host_core_color *
                        host_spec_color,  data = dat)),
  y = as.matrix(select(dat, log_asv_abundance, log_sample_count))
  )

dat_list$y1 <- dat$asv_abundance
dat_list$y2 <- dat$sample_count

dat_list$n_obs <- nrow(dat_list$x)

dat_list$J <- ncol(dat_list$x)
dat_list$K <- ncol(dat_list$y)

## Compile Stan model (this only needs to be done once) ----
mod <- cmdstan_model("model.stan")

## Fit model ----
fit <- mod$sample(
  data = dat_list,
  seed = 89638177,
  chains = 4,
  parallel_chains = 4,
  refresh = 500
)

## Save fitted model ----
fit$save_object(file = "fit_v3.rds")
