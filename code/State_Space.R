# State Space Modeling for Otolith LA-ICPMS Data
# ============================================
# This script implements a simple state space model to smooth 
# strontium isotope time series from otolith LA-ICPMS measurements

# Load required libraries
library(KFAS)      # For Kalman filtering and state space modeling
library(ggplot2)   # For visualization

# Read in the simulated data
data <- read.csv("simulated_isotope_data.csv")

# Look at the data structure
head(data)
summary(data$Iso)

# Prepare data for the state space model
# Create a time series object
iso_ts <- ts(data$Iso)

# Create a simple state space model structure using KFAS
# This model represents the data as a local level model 
# with random noise components for both the level and the observations
model_structure <- SSModel(
  iso_ts ~ SSMtrend(1, Q = list(NA)), 
  H = NA  # Observation variance (to be estimated)
)

# Estimate model parameters
model_fit <- fitSSM(model_structure, inits = c(0, 0))

# Extract estimated parameters
estimated_params <- model_fit$optim.out$par
names(estimated_params) <- c("log_observation_variance", "log_level_variance")
cat("Estimated parameters (log scale):", "\n")
print(estimated_params)
cat("Estimated parameters (original scale):", "\n")
print(exp(estimated_params))

# Update the model with estimated parameters
model_updated <- SSModel(
  iso_ts ~ SSMtrend(1, Q = list(exp(estimated_params[2]))), 
  H = exp(estimated_params[1])
)

# Run Kalman filter and smoother
kfs_result <- KFS(
  model_updated,
  simplify = TRUE,
  smoothing = c("state", "mean", "disturbance")
)

# Extract smoothed state estimates and confidence intervals
smoothed_mean <- kfs_result$alphahat[, 1]
smoothed_var <- kfs_result$V[1, 1, ]
confidence_interval <- 1.96 * sqrt(smoothed_var)

# Add results to data frame
data$smoothed <- smoothed_mean
data$lower_ci <- smoothed_mean - confidence_interval
data$upper_ci <- smoothed_mean + confidence_interval

# Simple visualization
p <- ggplot(data) +
  # Raw data points in gray
  geom_point(aes(x = Time, y = Iso), alpha = 0.3, color = "gray40") +
  # Smoothed state in blue
  geom_line(aes(x = Time, y = smoothed), color = "blue", size = 1) +
  # 95% confidence interval as shaded region
  geom_ribbon(aes(x = Time, ymin = lower_ci, ymax = upper_ci), 
              alpha = 0.2, fill = "blue") +
  # If true signal is available in the data, add it for comparison
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(x = Time, y = BaseSignal), color = "red", linetype = "dashed")
  } +
  # Labels and theme
  labs(title = "State Space Model Smoothing of Otolith LA-ICPMS Data",
       subtitle = "Blue line: Smoothed state with 95% CI, Red dashed line: True signal (if available)",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

# Display the plot
print(p)

# Optional: Save the plot
# ggsave("state_space_smoothing_plot.png", p, width = 10, height = 6)

# Try a different state space model configuration
# Local linear trend model (includes slope component)
model_structure2 <- SSModel(
  iso_ts ~ SSMtrend(2, Q = list(NA, NA)), 
  H = NA
)

# Estimate parameters for the new model
model_fit2 <- fitSSM(model_structure2, inits = c(0, 0, 0))

# Extract estimated parameters
estimated_params2 <- model_fit2$optim.out$par
names(estimated_params2) <- c("log_observation_variance", 
                              "log_level_variance", 
                              "log_slope_variance")
cat("Estimated parameters for trend model (log scale):", "\n")
print(estimated_params2)

# Update the model with estimated parameters
model_updated2 <- SSModel(
  iso_ts ~ SSMtrend(2, Q = list(exp(estimated_params2[2]), 
                                exp(estimated_params2[3]))), 
  H = exp(estimated_params2[1])
)

# Run Kalman filter and smoother for the new model
kfs_result2 <- KFS(
  model_updated2,
  simplify = TRUE,
  smoothing = c("state", "mean", "disturbance")
)

# Extract smoothed state estimates (level component)
smoothed_mean2 <- kfs_result2$alphahat[, 1]
smoothed_var2 <- kfs_result2$V[1, 1, ]
confidence_interval2 <- 1.96 * sqrt(smoothed_var2)

# Add results to data frame
data$smoothed2 <- smoothed_mean2
data$lower_ci2 <- smoothed_mean2 - confidence_interval2
data$upper_ci2 <- smoothed_mean2 + confidence_interval2

# Compare the two state space models
p_compare <- ggplot(data) +
  geom_point(aes(x = Time, y = Iso), alpha = 0.1, color = "gray40") +
  geom_line(aes(x = Time, y = smoothed, color = "Local Level")) +
  geom_line(aes(x = Time, y = smoothed2, color = "Local Linear Trend")) +
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(x = Time, y = BaseSignal, color = "True Signal"), linetype = "dashed")
  } +
  scale_color_manual(name = "Model Type", 
                     values = c("Local Level" = "blue", 
                                "Local Linear Trend" = "green",
                                "True Signal" = "red")) +
  labs(title = "Comparison of State Space Models",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

# Display the comparison plot
print(p_compare)

# Optional: Save the comparison plot
# ggsave("state_space_model_comparison.png", p_compare, width = 10, height = 6)