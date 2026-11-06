library(brms)
library(tidyverse)

# 1. Predictors check --------------------------------------------------

# Visual exploration of predictors
par(mfrow = c(2,2))
hist(range_metrics$sympatry, breaks = 12, main="Degree of Sympatry", xlab="Range overlap (%)", col="darkseagreen3")
hist(range_metrics$symmetry, breaks = 12, main="Degree of Range Symmetry", xlab="Range symmetry (from 0 to 0.5)", col="darkseagreen4")
hist(range_metrics$age_myr, breaks = 12, main="Evolutionary age", xlab="Evolutionary age (myr)", col="steelblue")
hist(range_metrics$centroid_dist_km, breaks = 12, main="Distance Between Species Range Centroids", xlab="Geographical Distance (Km)", col="steelblue2")
par(mfrow = c(1,1))

predictors <- range_metrics[, c("sympatry", "symmetry", "age_myr", "centroid_dist_km")]
cor_matrix <- cor(predictors)
round(cor_matrix, 2)

pairs(predictors,
      main = "Predictor Correlation Matrix",
      pch = 19,
      col = adjustcolor("steelblue", alpha.f = 0.6))

# Visual comparison of environmental predictors

# Visual comparions of niche-related predictors

# 2. Route-scale model ----------------------------------------------------

# 1. Prepare and scale data
route_data_model <- route_data |>
  mutate(
    k = as.integer(k),
    m1 = as.integer(m1),
    m2 = as.integer(m2),
    N = as.integer(N),
    z_age = as.numeric(scale(age_myr)),
    z_sympatry = as.numeric(scale(sympatry))
  )


# Sanity checks for physical bounds
stopifnot(all(route_data_model$k >= (route_data_model$m1 + route_data_model$m2 - route_data_model$N)))
stopifnot(all(route_data_model$k <= route_data_model$m1))
stopifnot(all(route_data_model$k <= route_data_model$m2))
stopifnot(all(route_data_model$N > 0))


# 2. Define Fisher hyper likelihood in Stan
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


# 3. Define custom family for brms
fisher_hyper <- custom_family(
  "fisher_hyper",
  dpars = "mu",    # mu represents affinity (exp(log_alpha))
  links = "log",   # models log(alpha)
  type = "int",    # response 'k' is an integer
  vars = c("vint1[n]", "vint2[n]", "vint3[n]")  # maps to m1, m2, N
)


# 4. Fit the model in brms
brm.fit <- brm(
  formula = k | vint(m1, m2, N) ~ z_age + z_sympatry,
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


# 5. Model Summary and Diagnostics
summary(brm.fit)
fixef(brm.fit)
plot(brm.fit)


library(brms)
library(tidyverse)

# 1. Extract posterior draws for the coefficients
draws <- as_draws_df(brm.fit)

# 2. Create a prediction grid for scaled sympatry (holding z_age at 0)
sympatry_seq <- seq(min(route_data_model$z_sympatry), max(route_data_model$z_sympatry), length.out = 100)

# 3. Calculate the linear predictor for every MCMC draw across the grid
pred_grid <- expand_grid(
  z_sympatry = sympatry_seq,
  draw_id = seq_len(nrow(draws))
) %>%
  mutate(
    # Linear predictor equation: Intercept + b_z_sympatry * sympatry + b_z_age * 0
    fit_draw = draws$b_Intercept[draw_id] +
      draws$b_z_sympatry[draw_id] * z_sympatry +
      draws$b_z_age[draw_id] * 0
  )

# 4. Summarize mean and 95% credible intervals across all posterior draws
pred_summary <- pred_grid %>%
  group_by(z_sympatry) %>%
  summarise(
    fit = mean(fit_draw),
    lwr = quantile(fit_draw, 0.025),
    upr = quantile(fit_draw, 0.975),
    .groups = "drop"
  )

# 5. Plot using ggplot2
ggplot() +
  # Raw empirical route points
  geom_point(data = route_data_model, aes(x = z_sympatry, y = empirical_log_OR),
             alpha = 0.3, color = "gray40") +
  # Model fit regression line
  geom_line(data = pred_summary, aes(x = z_sympatry, y = fit), color = "darkblue", linewidth = 1) +
  # 95% Credible Interval ribbon
  geom_ribbon(data = pred_summary, aes(x = z_sympatry, ymin = lwr, ymax = upr),
              fill = "darkblue", alpha = 0.2) +
  labs(
    title = "Predicted Log-Odds Affinity vs. Sympatry",
    subtitle = "Calculated manually from brms posterior draws",
    x = "Degree of Sympatry (Scaled)",
    y = "Predicted Log-Odds Ratio (Affinity)"
  ) +
  theme_minimal()


# Define parameter vector for clean reuse in diagnostics/plots
target_pars <- c("b0", "b_cen_dist", "b_age", "b_symmetry", "b_sympatry")

# 6. Diagnostics & Summaries
# par(mfrow = c(1, 1))
# traceplot(out_routes, pars = target_pars)
#
# print(precis(out_routes, pars = target_pars))

# 7. Plotting Posterior Estimates vs True Parameters
# Order matches target_pars: b0, b_cen_dist, b_age, b_symmetry, b_sympatry

# plot(precis(out_routes, pars = target_pars),
#      xlim = c(-1, 1),
#      main = "Posterior Estimates vs. True Parameters")
#
# abline(v = 0, lty = 2, col = "gray50")

# 3. Stop-scale model --------------------------------------------------------

# 1. Join stop_data with range_metrics by species pair (sp1 and sp2)
stop_model_data <- stop_data |>
  left_join(range_metrics, by = c("sp1", "sp2"))

# 2. Scale predictors and ensure integer constraints
stop_model_data <- stop_model_data |>
  mutate(
    pair_id = paste(sp1, sp2, sep = "__"),
    k  = as.integer(k),
    m1 = as.integer(m1),
    m2 = as.integer(m2),
    N  = as.integer(N),
    z_age      = as.numeric(scale(age_myr)),
    z_sympatry = as.numeric(scale(sympatry))
  )

# 3. Sanity checks for physical bounds
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
  formula = k | vint(m1, m2, N) ~ z_age + z_sympatry + (1 | pair_id),
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
plot(brf_stop_fit)

# 1. Extract posterior draws
draws <- as_draws_df(brm_stop_fit)

# 2. Create a prediction grid across scaled sympatry (holding other predictors at 0)
sympatry_seq <- seq(min(stop_model_data$z_sympatry, na.rm = TRUE),
                    max(stop_model_data$z_sympatry, na.rm = TRUE), length.out = 100)

pred_grid <- expand_grid(
  z_sympatry = sympatry_seq,
  draw_id = seq_len(nrow(draws))
) %>%
  mutate(
    fit_draw = draws$b_Intercept[draw_id] +
      draws$b_z_sympatry[draw_id] * z_sympatry +
      draws$b_z_age[draw_id] * 0 # adjusted if other covariates are included
  )

# 3. Summarize mean and 95% credible intervals
pred_summary <- pred_grid %>%
  group_by(z_sympatry) %>%
  summarise(
    fit = mean(fit_draw),
    lwr = quantile(fit_draw, 0.025),
    upr = quantile(fit_draw, 0.975),
    .groups = "drop"
  )

# 4. Plot using ggplot2
ggplot() +
  geom_point(data = stop_model_data, aes(x = z_sympatry, y = empirical_log_OR),
             alpha = 0.2, color = "gray40") +
  geom_line(data = pred_summary, aes(x = z_sympatry, y = fit), color = "darkblue", linewidth = 1) +
  geom_ribbon(data = pred_summary, aes(x = z_sympatry, ymin = lwr, ymax = upr),
              fill = "darkblue", alpha = 0.2) +
  labs(
    title = "Stop-Scale Log-Odds Affinity vs. Sympatry",
    subtitle = "Bayesian NCHg model fit with 95% credible intervals",
    x = "Degree of Sympatry (Scaled)",
    y = "Predicted Log-Odds Ratio (Affinity)"
  ) +
  theme_minimal()
