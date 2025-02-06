library(tidyverse)

# Load data
data(mtcars)

# Fit a linear model
lm_fit <- lm(mpg ~ wt, data = mtcars)

# Summarize the model
summary(lm_fit)

# Visualize
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Frequentist Linear Regression",
       x = "Weight",
       y = "Miles Per Gallon")

bayesian_fit <- brm(mpg ~ wt, data = mtcars, prior = c(set_prior("normal(0, 10)", class = "b")))


# Create new data for predictions
new_data <- tibble(wt = seq(min(mtcars$wt), max(mtcars$wt), length.out = 100))

# Generate predictions for new data
preds <- posterior_predict(bayesian_fit, newdata = new_data)

# Compute means and credible intervals
new_data <- new_data %>%
  mutate(
    pred_mean = colMeans(preds),
    pred_lower = apply(preds, 2, quantile, probs = 0.025),
    pred_upper = apply(preds, 2, quantile, probs = 0.975)
  )

# Plot
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point() +
  geom_line(data = new_data, aes(x = wt, y = pred_mean), color = "blue") +
  geom_ribbon(data = new_data, aes(x = wt, ymin = pred_lower, ymax = pred_upper), alpha = 0.2, inherit.aes = FALSE) +
  labs(title = "Bayesian Linear Regression with Credible Intervals",
       x = "Weight (wt)",
       y = "Miles Per Gallon (mpg)")


# Extract posterior samples of coefficients
posterior_draws <- as_draws_df(bayesian_fit) %>%
  select(b_Intercept, b_wt)

# Create a grid of wt values
wt_grid <- tibble(wt = seq(min(mtcars$wt), max(mtcars$wt), length.out = 100))

# Generate predictions for each posterior draw
posterior_lines <- posterior_draws %>%
  expand_grid(wt_grid) %>%
  mutate(predicted_mpg = b_Intercept + b_wt * wt)

# Plot
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point() +
  geom_line(data = posterior_lines, aes(x = wt, y = predicted_mpg, group = interaction(b_Intercept, b_wt)), 
            alpha = 0.1, color = "blue") +
  labs(title = "Bayesian Linear Regression with Posterior Lines",
       x = "Weight (wt)",
       y = "Miles Per Gallon (mpg)")
