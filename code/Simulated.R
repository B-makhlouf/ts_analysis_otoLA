# Simulation of isotope data with known states
# Based on characteristics of otolith LA-ICPMS data
# Now with 5 states: core, three middle states, and marine

library(ggplot2)

set.seed(123) # For reproducibility

# Parameters for simulation
n_points <- 1172           # Same number of points as your data
time_start <- 177.636      # Starting time from your data
time_end <- 791.24         # Approximate ending time
time_step <- (time_end - time_start) / (n_points - 1)

# Create time sequence
times <- seq(from = time_start, by = time_step, length.out = n_points)

# Define the states (mean values)
marine_value <- 0.7092     # Marine value
# Create a more varied pattern for the three middle states
state1_mean <- 0.7065      # First middle state
state2_mean <- 0.7055      # Second middle state (lower than first)
state3_mean <- 0.7075      # Third middle state (higher than both previous)
core_value <- 0.7078       # Core value - halfway between state1 and marine

# Define transition points (in terms of index)
core_end <- round(n_points * 0.1)      # Core to State 1 transition
transition1 <- round(n_points * 0.3)    # State 1 to State 2 transition
transition2 <- round(n_points * 0.5)    # State 2 to State 3 transition
transition3 <- round(n_points * 0.7)    # State 3 to Marine transition

# Define transition widths (how many points for the transition to occur)
transition_width <- 80     # Points for gradual transition between states

# Function to create sigmoid transition
sigmoid <- function(x) {
  1 / (1 + exp(-0.1 * x))
}

# Create the base signal with smooth transitions
base_signal <- numeric(n_points)
true_state <- numeric(n_points)

for (i in 1:n_points) {
  if (i < core_end - transition_width/2) {
    # Core state
    base_signal[i] <- core_value
    true_state[i] <- 0
  }
  else if (i >= core_end - transition_width/2 && i <= core_end + transition_width/2) {
    # Transition from core to first state (sigmoid)
    t <- (i - (core_end - transition_width/2)) / transition_width
    s <- sigmoid(20 * (t - 0.5))  # Steepness of 20 gives a nice transition
    base_signal[i] <- core_value + s * (state1_mean - core_value)
    true_state[i] <- 0.5  # Mark transition zone
  }
  else if (i > core_end + transition_width/2 && i < transition1 - transition_width/2) {
    # First state
    base_signal[i] <- state1_mean
    true_state[i] <- 1
  } 
  else if (i >= transition1 - transition_width/2 && i <= transition1 + transition_width/2) {
    # Transition from first to second state (sigmoid)
    t <- (i - (transition1 - transition_width/2)) / transition_width
    s <- sigmoid(20 * (t - 0.5))  # Steepness of 20 gives a nice transition
    base_signal[i] <- state1_mean + s * (state2_mean - state1_mean)
    true_state[i] <- 1.5  # Mark transition zone
  }
  else if (i > transition1 + transition_width/2 && i < transition2 - transition_width/2) {
    # Second state
    base_signal[i] <- state2_mean
    true_state[i] <- 2
  }
  else if (i >= transition2 - transition_width/2 && i <= transition2 + transition_width/2) {
    # Transition from second to third state (sigmoid)
    t <- (i - (transition2 - transition_width/2)) / transition_width
    s <- sigmoid(20 * (t - 0.5))
    base_signal[i] <- state2_mean + s * (state3_mean - state2_mean)
    true_state[i] <- 2.5  # Mark transition zone
  }
  else if (i > transition2 + transition_width/2 && i < transition3 - transition_width/2) {
    # Third state
    base_signal[i] <- state3_mean
    true_state[i] <- 3
  }
  else if (i >= transition3 - transition_width/2 && i <= transition3 + transition_width/2) {
    # Transition from third state to marine (sigmoid)
    t <- (i - (transition3 - transition_width/2)) / transition_width
    s <- sigmoid(20 * (t - 0.5))
    base_signal[i] <- state3_mean + s * (marine_value - state3_mean)
    true_state[i] <- 3.5  # Mark transition zone
  }
  else {
    # Marine state
    base_signal[i] <- marine_value
    true_state[i] <- 4
  }
}

# Add two types of noise
# 1. White noise (random variation around the mean)
white_noise_sd <- 0.0005 
white_noise <- rnorm(n_points, mean = 0, sd = white_noise_sd)

# 2. Autocorrelated noise (to mimic instrument drift or natural processes)
# This creates more realistic isotope data
ar_coefficient <- 0.7       # Autocorrelation coefficient
ar_noise_sd <- 0.0005       # Standard deviation of the AR process innovations
ar_noise <- numeric(n_points)
ar_noise[1] <- rnorm(1, 0, ar_noise_sd)

for (i in 2:n_points) {
  ar_noise[i] <- ar_coefficient * ar_noise[i-1] + rnorm(1, 0, ar_noise_sd)
}

# Combine base signal with both noise components
iso_values <- base_signal + white_noise + ar_noise

# Create the simulated dataset
sim_data <- data.frame(
  Cycle = 1:n_points,
  Time = times,
  Iso = iso_values,
  BaseSignal = base_signal,
  TrueState = true_state
)

# Create a factor version of the state variable for plotting
sim_data$StateFactor <- factor(floor(sim_data$TrueState), 
                               levels = c("0", "1", "2", "3", "4"),
                               labels = c("Core", "State 1", "State 2", "State 3", "Marine"))

# Plot the simulated data
p1 <- ggplot(sim_data, aes(x = Time, y = Iso)) +
  geom_point(color = "gray40") +
  geom_point(aes(y = BaseSignal), color = "blue", size = 1) +
  geom_vline(xintercept = c(times[core_end], times[transition1], times[transition2], times[transition3]), 
             linetype = "dashed", color = "red") +
  labs(title = "Simulated Otolith Isotope Data with 5 States", 
       subtitle = "Gray: Noisy signal, Blue: True underlying states, Red: Transition points",
       x = "Time (Distance from core)", y = "87Sr/86Sr") +
  theme_minimal()

print(p1)

# Plot showing states as colors
p2 <- ggplot(sim_data, aes(x = Time, y = Iso, color = StateFactor)) +
  geom_point() +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red"), 
                     name = "State") +
  labs(title = "Simulated Otolith Isotope Data with 5 States", 
       subtitle = paste("Core:", core_value, "→ State 1:", state1_mean, 
                        "→ State 2:", state2_mean, "→ State 3:", state3_mean, 
                        "→ Marine:", marine_value),
       x = "Time (Distance from core)", y = "87Sr/86Sr") +
  theme_minimal()

print(p2)

# Save the simulated data to a CSV file
write.csv(sim_data, "/Users/benjaminmakhlouf/Research_repos/ts_analysis_otoLA/data/simulated Data/sim01.csv", row.names = FALSE)

