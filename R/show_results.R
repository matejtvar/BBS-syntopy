library(bayesplot)
color_scheme_set("lightblue")
mcmc_hist(brm_stop_fit, pars = c("b_Intercept", "b_z_age", "b_z_sympatry"))


# 1. Extract posterior draws from your 4-predictor route model
draws <- as_draws_df(brm.fit)

# 2. Set up the original-scale sequence for the predictor you want to plot (e.g., sympatry)
# (You can swap 'sympatry' for 'age_myr', 'symmetry', or 'centroid_dist_km' to plot the others)
x_var_name <- "sympatry"
x_orig_seq <- seq(min(route_data_model[[x_var_name]], na.rm = TRUE),
                  max(route_data_model[[x_var_name]], na.rm = TRUE),
                  length.out = 100)

# Get mean and sd used to scale that specific variable originally
x_mean <- mean(route_data_model[[x_var_name]], na.rm = TRUE)
x_sd   <- sd(route_data_model[[x_var_name]], na.rm = TRUE)

# 3. Create prediction grid across draws
# All other scaled covariates are set to 0 (their mean value)
pred_grid_route <- expand_grid(
  x_orig = x_orig_seq,
  draw_id = seq_len(nrow(draws))
) |>
  mutate(
    z_val = (x_orig - x_mean) / x_sd,

    # Linear predictor calculation including all four predictors
    fit_draw = draws$b_Intercept[draw_id] +
      draws$b_z_sympatry[draw_id] * z_val +  # (matches your target variable)
      draws$b_z_age[draw_id] * 0 +          # held at mean
      draws$b_z_symmetry[draw_id] * 0 +       # held at mean
      draws$b_z_cen_dist[draw_id] * 0         # held at mean
  )

# 4. Summarize mean and 95% credible intervals across posterior draws
pred_summary_route <- pred_grid_route %>%
  group_by(x_orig) %>%
  summarise(
    fit = mean(fit_draw),
    lwr = quantile(fit_draw, 0.025),
    upr = quantile(fit_draw, 0.975),
    .groups = "drop"
  )


ggplot() +
  # Zero reference line (independence)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  # Raw route points (or pair-level summaries if you prefer)
  geom_point(data = route_data_model,
             aes(.data[[x_var_name]], empirical_log_OR),
             color = "gray40", size = 1.5, alpha = 0.4) +
  # Bayesian regression line from the full model
  geom_line(data = pred_summary_route,
            aes(x = x_orig, y = fit),
            color = "#2b7bb9", linewidth = 1) +
  # 95% Credible Interval ribbon
  geom_ribbon(data = pred_summary_route,
              aes(x = x_orig, ymin = lwr, ymax = upr),
              fill = "#2b7bb9", alpha = 0.2) +
  labs(
    x = "Sympatry (%)", # Change label depending on the active variable
    y = "Syntopy (Log-Odds Ratio / Affinity)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid = element_blank()
  )
