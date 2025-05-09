# Hidden Markov Model (HMM) for Otolith LA-ICPMS Data
# =================================================
# This script implements a Hidden Markov Model to detect discrete 
# states in strontium isotope time series from otolith LA-ICPMS data

# Load required libraries
library(depmixS4)  # For HMM modeling
library(ggplot2)   # For visualization
library(dplyr)     # For data manipulation

# Read in the simulated data
data <- read.csv("simulated_isotope_data.csv")

# Look at the data structure
head(data)
summary(data$Iso)

# Visualize the raw data
p_raw <- ggplot(data, aes(x = Time, y = Iso)) +
  geom_point(alpha = 0.3, color = "gray40") +
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(y = BaseSignal), color = "red", linetype = "dashed")
  } +
  labs(title = "Raw Isotope Data",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

print(p_raw)

# 1. HMM with a specified number of states
# ----------------------------------------
# We'll try 5 states to match our simulated data
# (in a real analysis, you might need to try different values)
n_states <- 5

# Set up and fit the HMM
# response = gaussian() means we're modeling the isotope values as normal distributions
hmm_model <- depmix(Iso ~ 1, data = data, nstates = n_states, 
                    family = gaussian(), ntimes = nrow(data))

# Fit the model
# seed is set for reproducibility
# maxit is the maximum number of iterations
set.seed(123)
hmm_fit <- fit(hmm_model, verbose = TRUE, maxit = 1000)

# Print model summary
summary(hmm_fit)

# Get state assignments
state_probs <- posterior(hmm_fit)
data$hmm_state <- state_probs$state

# Extract model parameters
# Get means and standard deviations for each state
state_params <- data.frame(
  State = 1:n_states,
  Mean = sapply(1:n_states, function(i) getpars(hmm_fit)[grep(paste0("^Re.*\\[", i), getpars(hmm_fit, which = "response"))][1]),
  SD = sapply(1:n_states, function(i) getpars(hmm_fit)[grep(paste0("^Re.*\\[", i), getpars(hmm_fit, which = "response"))][2])
)

# Sort states by mean value for better interpretation
state_params <- state_params[order(state_params$Mean), ]
state_params$SortedState <- 1:n_states

# Create a mapping from original to sorted states
state_mapping <- state_params$SortedState
names(state_mapping) <- state_params$State

# Apply the mapping to get sorted states
data$sorted_state <- state_mapping[data$hmm_state]

# Print state parameters
print("Estimated state parameters:")
print(state_params)

# Visualize the HMM results
p_hmm <- ggplot(data) +
  geom_point(aes(x = Time, y = Iso, color = factor(sorted_state)), alpha = 0.6) +
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(x = Time, y = BaseSignal), color = "black", linetype = "dashed")
  } +
  scale_color_brewer(name = "HMM State", palette = "Set1") +
  labs(title = paste("HMM Classification with", n_states, "States"),
       subtitle = "Colors indicate different detected states",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

print(p_hmm)

# 2. Try a range of state numbers to find the optimal model
# --------------------------------------------------------
# Define range of states to try
state_range <- 2:8

# Initialize storage for BIC values (for model comparison)
bic_values <- numeric(length(state_range))

cat("Fitting models with different numbers of states...\n")

# Fit models with different numbers of states
for(i in seq_along(state_range)) {
  n <- state_range[i]
  cat("Fitting model with", n, "states...\n")
  
  # Set up model
  temp_model <- depmix(Iso ~ 1, data = data, nstates = n, 
                       family = gaussian(), ntimes = nrow(data))
  
  # Fit the model
  set.seed(123)
  temp_fit <- try(fit(temp_model, verbose = FALSE, maxit = 1000), silent = TRUE)
  
  # Store BIC value if model fit successfully
  if(!inherits(temp_fit, "try-error")) {
    bic_values[i] <- BIC(temp_fit)
    cat("BIC:", bic_values[i], "\n")
  } else {
    bic_values[i] <- NA
    cat("Model fitting failed\n")
  }
}

# Create a data frame for BIC values
bic_df <- data.frame(
  States = state_range,
  BIC = bic_values
)

# Find the number of states with the lowest BIC
best_n_states <- state_range[which.min(bic_values)]
cat("Best model according to BIC:", best_n_states, "states\n")

# Plot BIC values
p_bic <- ggplot(bic_df, aes(x = States, y = BIC)) +
  geom_line() +
  geom_point() +
  geom_point(data = bic_df[which.min(bic_values), ], color = "red", size = 3) +
  labs(title = "Model Selection: BIC by Number of States",
       subtitle = paste("Optimal number of states:", best_n_states),
       x = "Number of States",
       y = "BIC (lower is better)") +
  theme_minimal()

print(p_bic)

# 3. Fit the optimal model
# ------------------------
if(best_n_states != n_states) {
  cat("Fitting the optimal model with", best_n_states, "states...\n")
  
  # Set up and fit the optimal model
  optimal_model <- depmix(Iso ~ 1, data = data, nstates = best_n_states, 
                          family = gaussian(), ntimes = nrow(data))
  
  # Fit the model
  set.seed(123)
  optimal_fit <- fit(optimal_model, verbose = TRUE, maxit = 1000)
  
  # Get state assignments
  optimal_state_probs <- posterior(optimal_fit)
  data$optimal_state <- optimal_state_probs$state
  
  # Extract model parameters
  optimal_state_params <- data.frame(
    State = 1:best_n_states,
    Mean = sapply(1:best_n_states, function(i) 
      getpars(optimal_fit)[grep(paste0("^Re.*\\[", i), 
                                getpars(optimal_fit, which = "response"))][1]),
    SD = sapply(1:best_n_states, function(i) 
      getpars(optimal_fit)[grep(paste0("^Re.*\\[", i), 
                                getpars(optimal_fit, which = "response"))][2])
  )
  
  # Sort states by mean value
  optimal_state_params <- optimal_state_params[order(optimal_state_params$Mean), ]
  optimal_state_params$SortedState <- 1:best_n_states
  
  # Create mapping from original to sorted states
  optimal_state_mapping <- optimal_state_params$SortedState
  names(optimal_state_mapping) <- optimal_state_params$State
  
  # Apply the mapping
  data$optimal_sorted_state <- optimal_state_mapping[data$optimal_state]
  
  # Print optimal state parameters
  print("Estimated state parameters for optimal model:")
  print(optimal_state_params)
  
  # Visualize the optimal HMM results
  p_optimal <- ggplot(data) +
    geom_point(aes(x = Time, y = Iso, color = factor(optimal_sorted_state)), alpha = 0.6) +
    {if("BaseSignal" %in% names(data)) 
      geom_line(aes(x = Time, y = BaseSignal), color = "black", linetype = "dashed")
    } +
    scale_color_brewer(name = "HMM State", palette = "Set1") +
    labs(title = paste("Optimal HMM Classification with", best_n_states, "States"),
         subtitle = "Colors indicate different detected states",
         x = "Time (Distance from core)",
         y = "87Sr/86Sr") +
    theme_minimal()
  
  print(p_optimal)
}

# 4. Compare with true states (if available in simulated data)
# -----------------------------------------------------------
if("TrueState" %in% names(data)) {
  # Floor the true state to get the base state (removing transition markers)
  data$true_base_state <- floor(data$TrueState)
  
  # Create a comparative visualization
  p_compare <- ggplot(data) +
    geom_point(aes(x = Time, y = Iso, color = factor(true_base_state)), 
               alpha = 0.3, size = 2) +
    geom_point(aes(x = Time, y = Iso - 0.001, shape = factor(optimal_sorted_state)), 
               alpha = 0.7, size = 2) +
    scale_color_brewer(name = "True State", palette = "Set1") +
    scale_shape_discrete(name = "HMM State") +
    labs(title = "Comparison of True States vs. HMM-Detected States",
         subtitle = "Colors show true states, shapes show HMM-detected states",
         x = "Time (Distance from core)",
         y = "87Sr/86Sr") +
    theme_minimal()
  
  print(p_compare)
  
  # Calculate the confusion matrix to evaluate model performance
  # Create a table of true vs. HMM states
  confusion <- table(True = data$true_base_state, 
                     HMM = data$optimal_sorted_state)
  
  print("Confusion matrix (rows = true states, columns = HMM states):")
  print(confusion)
  
  # Calculate accuracy
  accuracy <- sum(diag(confusion)) / sum(confusion)
  cat("Accuracy:", round(accuracy * 100, 1), "%\n")
}

# Optional: Save key plots
# ggsave("hmm_optimal_classification.png", p_optimal, width = 10, height = 6)
# ggsave("hmm_model_selection.png", p_bic, width = 8, height = 5)