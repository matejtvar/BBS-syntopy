library(brms)
library(tidyverse)
library(GGally)

# targets::tar_load_everything()
# list2env(c(data_import, syntopy_predictors), envir = environment())


# 1. Predictors check --------------------------------------------------

# Visual exploration of predictors
# par(mfrow = c(2,2))
# hist(range_metrics$sympatry, breaks = 12, main="Degree of Sympatry", xlab="Range overlap (%)", col="darkseagreen3")
# hist(range_metrics$symmetry, breaks = 12, main="Degree of Range Symmetry", xlab="Range symmetry (from 0 to 0.5)", col="darkseagreen4")
# hist(range_metrics$age_myr, breaks = 12, main="Evolutionary age", xlab="Evolutionary age (myr)", col="steelblue")
# hist(range_metrics$centroid_dist_km, breaks = 12, main="Distance Between Species Range Centroids", xlab="Geographical Distance (Km)", col="steelblue2")
# par(mfrow = c(1,1))
#
# predictors <- range_metrics[, c("sympatry", "symmetry", "age_myr", "centroid_dist_km")]
# cor_matrix <- cor(predictors)
# round(cor_matrix, 2)
#
# cen_dist <- log(range_metrics$centroid_dist_km)
# age <- log(range_metrics$age_myr)
#
# GGally::ggpairs(cbind(predictors[,1:2], cen_dist, age))

# Visual comparison of environmental predictors

# Visual comparions of niche-related predictors

# 2. Route-scale model ----------------------------------------------------

route_predictors <- range_metrics |>
  dplyr::select(sympatry, symmetry, centroid_dist_km)

route_data <- route_data |>
  dplyr::left_join(route_predictors, by = "sympatry") |>
  dplyr::select(sp1, sp2, sympatry, symmetry, centroid_dist_km, age_myr, N, k, m1, m2)

head(route_data)


# Prepare and scale data
route_data_model <- route_data |>
  mutate(
    k = as.integer(k),
    m1 = as.integer(m1),
    m2 = as.integer(m2),
    N = as.integer(N),
    z_age = as.numeric(scale(age_myr)),
    z_sympatry = as.numeric(scale(sympatry)),
    z_symmetry = as.numeric(scale(symmetry)),
    z_cen_dist = as.numeric(scale(symmetry))
  )


# Sanity checks for physical bounds
stopifnot(all(route_data_model$k >= (route_data_model$m1 + route_data_model$m2 - route_data_model$N)))
stopifnot(all(route_data_model$k <= route_data_model$m1))
stopifnot(all(route_data_model$k <= route_data_model$m2))
stopifnot(all(route_data_model$N > 0))


# Define Fisher hyper likelihood in Stan
stan_fisher <- stanvar(
  scode = "
  real fisher_hyper_lpmf(int k, real mu, int m1, int m2, int N) {
    if (mu <= 0) return not_a_number();

    // Bounds for hypergeometric distribution (renamed to avoid reserved keywords)
    int lower_bound = max({0, m1 + m2 - N});
    int upper_bound = min({m1, m2});

    // Log-numerator for the observed value k
    real log_num = lchoose(m1, k) + lchoose(N - m1, m2 - k) + k * log(mu);

    // Normalization constant calculation
    int range_size = upper_bound - lower_bound + 1;
    vector[range_size] terms;

    for (i in 0:(range_size - 1)) {
      int curr_i = lower_bound + i;
      terms[i + 1] = lchoose(m1, curr_i) + lchoose(N - m1, m2 - curr_i) + curr_i * log(mu);
    }

    return log_num - log_sum_exp(terms);
  }
  ",
  block = "functions",
  name = "fisher_hyper_lpmf"
)


# Define custom family for brms
fisher_hyper <- custom_family(
  "fisher_hyper",
  dpars = "mu",    # mu represents affinity (exp(log_alpha))
  links = "log",   # models log(alpha)
  type = "int",    # response 'k' is an integer
  vars = c("vint1[n]", "vint2[n]", "vint3[n]")  # maps to m1, m2, N
)


# Fit the model in brms
brm.fit <- brm(
  formula = k | vint(m1, m2, N) ~ z_age + z_sympatry + z_symmetry + z_cen_dist,
  data = route_data_model,
  family = fisher_hyper,
  stanvars = stan_fisher,
  prior = c(
    prior(normal(0, 3), class = "b"),        # Prior for predictors
    prior(normal(0, 3), class = "Intercept") # Prior for the Intercept
  ),
  init = 0,     # Initializes all parameters at 0 for safe start
  chains = 4,
  iter = 3000,
  warmup = 2000,
  cores = 4
)

summary(brm.fit)
fixef(brm.fit)
plot(brm.fit)


# 3. Stop-scale model --------------------------------------------------------

# Join stop_data with range_metrics by species pair (sp1 and sp2)
stop_model_data <- stop_data |>
  dplyr::left_join(range_metrics, by = c("sp1", "sp2")) |>
  dplyr::filter(m1 > 0 & m2 > 0)

# Scale predictors and ensure integer constraints
stop_model_data <- stop_model_data |>
  mutate(
    pair_id = paste(sp1, sp2, sep = "__"),
    k  = as.integer(k),
    m1 = as.integer(m1),
    m2 = as.integer(m2),
    N  = as.integer(N),
    z_age      = as.numeric(scale(age_myr)),
    z_sympatry = as.numeric(scale(sympatry)),
    z_symmetry = as.numeric(scale(symmetry))
  )

# Sanity checks for physical bounds
stopifnot(all(stop_model_data$k >= (stop_model_data$m1 + stop_model_data$m2 - stop_model_data$N)))
stopifnot(all(stop_model_data$k <= stop_model_data$m1))
stopifnot(all(stop_model_data$k <= stop_model_data$m2))
stopifnot(all(stop_model_data$N > 0))

# Define custom family
fisher_hyper <- custom_family(
  "fisher_hyper",
  dpars = "mu",
  links = "log",
  type = "int",
  vars = c("vint1[n]", "vint2[n]", "vint3[n]")
)

brm_stop_fit <- brm(
  formula = k | vint(m1, m2, N) ~ z_age + z_sympatry + z_symmetry + (1 | pair_id),
  data = stop_model_data,
  family = fisher_hyper,
  stanvars = stan_fisher,
  prior = c(
    prior(normal(0, 3), class = "b"),
    prior(normal(0, 3), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  chains = 4,
  iter = 3000,
  warmup = 2000,
  cores = 4
)

# Summary of results
summary(brm_stop_fit)
fixef(brm_stop_fit)
plot(brm_stop_fit)

# Productivity
stop_data_2 <- stop_data |>
  dplyr::left_join(range_metrics, by = c("sp1", "sp2")) |>
  dplyr::filter(m1 > 0 & m2 > 0) |>
  dplyr::filter(is.na(NDVI_filled) == FALSE)

hist(stop_data_2$NDVI_filled)

# Scale predictors and ensure integer constraints
stop_data_2 <- stop_data_2 |>
  mutate(
    pair_id = paste(sp1, sp2, sep = "__"),
    k  = as.integer(k),
    m1 = as.integer(m1),
    m2 = as.integer(m2),
    N  = as.integer(N),
    z_age      = as.numeric(scale(age_myr)),
    z_sympatry = as.numeric(scale(sympatry)),
    z_symmetry = as.numeric(scale(symmetry)),
    z_ndvi = as.numeric(scale(NDVI_filled))
  )

# Sanity checks for physical bounds
stopifnot(all(stop_data_2$k >= (stop_data_2$m1 + stop_data_2$m2 - stop_data_2$N)))
stopifnot(all(stop_data_2$k <= stop_data_2$m1))
stopifnot(all(stop_data_2$k <= stop_data_2$m2))
stopifnot(all(stop_data_2$N > 0))

# Define custom family
fisher_hyper <- custom_family(
  "fisher_hyper",
  dpars = "mu",
  links = "log",
  type = "int",
  vars = c("vint1[n]", "vint2[n]", "vint3[n]")
)

stop_model_2 <- brm(
  formula = k | vint(m1, m2, N) ~ z_age + z_sympatry + z_symmetry + z_prod + (1 | pair_id),
  data = stop_data_2,
  family = fisher_hyper,
  stanvars = stan_fisher,
  prior = c(
    prior(normal(0, 3), class = "b"),
    prior(normal(0, 3), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  chains = 4,
  iter = 3000,
  warmup = 2000,
  cores = 4
)

summary(stop_model_2)
fixef(stop_model_2)
plot(stop_model_2)



list(
  route_model = brm.fit,
  stop_model = brm_stop_fit
)
