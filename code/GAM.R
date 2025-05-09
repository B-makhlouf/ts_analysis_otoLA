# GAM Smoothing for Otolith LA-ICPMS Data
# ==================================
# This script implements a basic Generalized Additive Model (GAM) for smoothing 
# strontium isotope time series from otolith LA-ICPMS measurements

# Load required libraries
library(mgcv)     # For GAM modeling
library(ggplot2)  # For visualization

# Read in the simulated data
data <- read.csv("simulated_isotope_data.csv")

# Look at the data structure
head(data)
summary(data$Iso)

# Fit GAM model to isotope data
# The 'k' parameter controls smoothness (higher k = more flexible curve)
# Method "REML" (Restricted Maximum Likelihood) helps prevent overfitting
gam_model <- gam(Iso ~ s(Time, k=30), data = data, method = "REML")

# Display basic model summary
summary(gam_model)

# Extract model predictions
predictions <- predict(gam_model, se.fit = TRUE)

# Add predictions to the original data
data$fitted <- predictions$fit
data$se <- predictions$se.fit
data$lower_ci <- data$fitted - 1.96 * data$se
data$upper_ci <- data$fitted + 1.96 * data$se

# Simple visualization
p <- ggplot(data) +
  # Raw data points in gray
  geom_point(aes(x = Time, y = Iso), alpha = 0.3, color = "gray40") +
  # GAM fit in blue
  geom_line(aes(x = Time, y = fitted), color = "blue", size = 1) +
  # 95% confidence interval as shaded region
  geom_ribbon(aes(x = Time, ymin = lower_ci, ymax = upper_ci), 
              alpha = 0.2, fill = "blue") +
  # If true signal is available in the data, add it for comparison
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(x = Time, y = BaseSignal), color = "red", linetype = "dashed")
  } +
  # Labels and theme
  labs(title = "GAM Smoothing of Otolith LA-ICPMS Data",
       subtitle = "Blue line: GAM estimate with 95% CI, Red dashed line: True signal (if available)",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

# Display the plot
print(p)

# Optional: Save the plot
# ggsave("gam_smoothing_plot.png", p, width = 10, height = 6)

# Alternative visualization with different smoothing parameters
# Try a few different values for the smoothing parameter k
k_values <- c(10, 30, 60, 100)

# Create an empty list to store the models
models_list <- list()

# Fit models with different k values
for (k in k_values) {
  models_list[[as.character(k)]] <- gam(Iso ~ s(Time, k=k), data = data, method = "REML")
}

# Extract predictions from each model
for (k in k_values) {
  k_str <- as.character(k)
  pred <- predict(models_list[[k_str]], se.fit = TRUE)
  data[[paste0("fitted_k", k)]] <- pred$fit
}

# Plot comparing different smoothing levels
p_compare <- ggplot(data) +
  geom_point(aes(x = Time, y = Iso), alpha = 0.1, color = "gray40") +
  geom_line(aes(x = Time, y = fitted_k10, color = "k=10")) +
  geom_line(aes(x = Time, y = fitted_k30, color = "k=30")) +
  geom_line(aes(x = Time, y = fitted_k60, color = "k=60")) +
  geom_line(aes(x = Time, y = fitted_k100, color = "k=100")) +
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(x = Time, y = BaseSignal, color = "True Signal"), linetype = "dashed")
  } +
  scale_color_manual(name = "Smoothing Level", 
                     values = c("k=10" = "blue", "k=30" = "green", 
                                "k=60" = "purple", "k=100" = "orange",
                                "True Signal" = "red")) +
  labs(title = "Comparison of GAM Smoothing Levels",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

# Display the comparison plot
print(p_compare)

# Optional: Save the comparison plot
# ggsave("gam_smoothing_comparison.png", p_compare, width = 10, height = 6)