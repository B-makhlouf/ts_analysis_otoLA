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
write.csv(sim_data, "simulated_isotope_data.csv", row.names = FALSE)

# Calculate summary statistics
cat("Summary statistics for simulated data:\n")
print(summary(sim_data$Iso))
cat("\nStandard deviation of Iso values:", sd(sim_data$Iso), "\n")

cat("\nTime range:", min(sim_data$Time), "to", max(sim_data$Time), "\n")
cat("Transition points at time:", 
    times[core_end], "(Core → State 1),",
    times[transition1], "(State 1 → State 2),", 
    times[transition2], "(State 2 → State 3),",
    times[transition3], "(State 3 → Marine)\n")


# Advanced option: Adding machine blurring effect
# This simulates how analytical instruments might blur sharp transitions
# due to signal integration, washout effects, etc.

# Create a version with machine blurring
convolve_window <- function(x, window_size) {
  n <- length(x)
  y <- numeric(n)
  
  half_window <- floor(window_size/2)
  
  for (i in 1:n) {
    start_idx <- max(1, i - half_window)
    end_idx <- min(n, i + half_window)
    y[i] <- mean(x[start_idx:end_idx])
  }
  
  return(y)
}

# Apply machine blurring to the simulated isotope values
window_size <- 5  # Size of the moving average window
blurred_iso <- convolve_window(iso_values, window_size)

# Add to the dataset
sim_data$BlurredIso <- blurred_iso

# Plot the data with machine blurring
p3 <- ggplot(sim_data) +
  geom_point(aes(x = Time, y = Iso), color = "gray40", alpha = 0.5) +
  geom_point(aes(x = Time, y = BlurredIso), color = "red") +
  geom_point(aes(x = Time, y = BaseSignal), color = "blue", size = 1) +
  labs(title = "Simulated Isotope Data with Machine Blurring", 
       subtitle = "Gray: Raw signal, Red: Blurred signal, Blue: True underlying states",
       x = "Time (Distance from core)", y = "87Sr/86Sr") +
  theme_minimal()

print(p3)

# Save this version as well
write.csv(sim_data, "simulated_isotope_data_with_blurring.csv", row.names = FALSE)

# Create versions with varying signal-to-noise ratios
create_simulation_with_SNR <- function(base_signal, snr) {
  # Calculate the average signal power
  signal_power <- mean(base_signal^2)
  
  # Calculate required noise power for the desired SNR
  noise_power <- signal_power / snr
  noise_sd <- sqrt(noise_power)
  
  # Generate noise
  white_noise <- rnorm(length(base_signal), mean = 0, sd = noise_sd * 0.5)
  
  # AR noise component
  ar_noise <- numeric(length(base_signal))
  ar_noise[1] <- rnorm(1, 0, noise_sd * 0.5)
  for (i in 2:length(base_signal)) {
    ar_noise[i] <- ar_coefficient * ar_noise[i-1] + rnorm(1, 0, noise_sd * 0.5)
  }
  
  # Combine
  return(base_signal + white_noise + ar_noise)
}

# Create versions with different SNRs
snr_values <- c(50, 20, 10, 5)
for (snr in snr_values) {
  sim_data[[paste0("Iso_SNR_", snr)]] <- create_simulation_with_SNR(base_signal, snr)
}

# Plot the different SNR versions
snr_data <- reshape2::melt(sim_data[, c("Time", paste0("Iso_SNR_", snr_values))], 
                           id.vars = "Time", 
                           variable.name = "SNR", 
                           value.name = "Iso")

p4 <- ggplot(snr_data, aes(x = Time, y = Iso, color = SNR)) +
  geom_point() +
  facet_wrap(~SNR, ncol = 1) +
  labs(title = "Simulated Isotope Data with Different Signal-to-Noise Ratios", 
       x = "Time (Distance from core)", y = "87Sr/86Sr") +
  theme_minimal() +
  theme(legend.position = "none")

print(p4)

# Save the data with different SNRs
write.csv(sim_data, "simulated_isotope_data_varying_SNR.csv", row.names = FALSE)

# Add a specific visualization of the otolith pattern with arrows showing movement
p5 <- ggplot() +
  geom_point(data = sim_data, aes(x = Time, y = BaseSignal), color = "blue", size = 1) +
  geom_hline(yintercept = marine_value, pointtype = "dashed", color = "red") +
  geom_hline(yintercept = core_value, pointtype = "dashed", color = "purple") +
  # Add arrows to show the non-sequential pattern
  geom_segment(aes(x = times[core_end+50], y = core_value, 
                   xend = times[core_end+70], yend = state1_mean),
               arrow = arrow(length = unit(0.3, "cm")), color = "darkgray") +
  geom_segment(aes(x = times[transition1+50], y = state1_mean, 
                   xend = times[transition1+70], yend = state2_mean),
               arrow = arrow(length = unit(0.3, "cm")), color = "darkgray") +
  geom_segment(aes(x = times[transition2+50], y = state2_mean, 
                   xend = times[transition2+70], yend = state3_mean),
               arrow = arrow(length = unit(0.3, "cm")), color = "darkgray") +
  geom_segment(aes(x = times[transition3+50], y = state3_mean, 
                   xend = times[transition3+70], yend = marine_value),
               arrow = arrow(length = unit(0.3, "cm")), color = "darkgray") +
  # Add annotations
  annotate("text", x = min(sim_data$Time) + 30, y = core_value + 0.0003, 
           label = paste("Core value (", sprintf("%.4f", core_value), ")"), color = "purple") +
  annotate("text", x = max(sim_data$Time) - 50, y = marine_value + 0.0003, 
           label = paste("Marine value (", sprintf("%.4f", marine_value), ")"), color = "red") +
  annotate("text", x = times[round(mean(c(core_end, transition1)))], y = state1_mean + 0.0003,
           label = paste("State 1 (", sprintf("%.4f", state1_mean), ")"), color = "blue") +
  annotate("text", x = times[round(mean(c(transition1, transition2)))], y = state2_mean - 0.0003,
           label = paste("State 2 (", sprintf("%.4f", state2_mean), ")"), color = "green4") +
  annotate("text", x = times[round(mean(c(transition2, transition3)))], y = state3_mean + 0.0003,
           label = paste("State 3 (", sprintf("%.4f", state3_mean), ")"), color = "orange3") +
  labs(title = "Simulated Otolith 87Sr/86Sr Pattern with Non-Sequential States", 
       subtitle = "The fish moves through habitats with varying 87Sr/86Sr values before reaching the marine environment",
       x = "Time (Distance from core)", 
       y = "87Sr/86Sr") +
  theme_minimal() +
  ylim(min(c(state1_mean, state2_mean, state3_mean, core_value, marine_value)) - 0.0005,
       max(c(state1_mean, state2_mean, state3_mean, core_value, marine_value)) + 0.0005)

print(p5)

# Save this visualization
ggsave("otolith_simulated_pattern.png", p5, width = 10, height = 6)

# Print a message about what was created
cat("\nSimulation complete. Files created:\n")
cat("1. simulated_isotope_data.csv - Basic simulation\n")
cat("2. simulated_isotope_data_with_blurring.csv - With machine blurring effect\n") 
cat("3. simulated_isotope_data_varying_SNR.csv - With different signal-to-noise ratios\n")

# Print information about the states
cat("\nState values used in simulation:\n")
cat("Core state:", core_value, "(Initial otolith core)\n")
cat("State 1:", state1_mean, "(First habitat)\n")
cat("State 2:", state2_mean, "(Second habitat)\n")
cat("State 3:", state3_mean, "(Third habitat)\n")
cat("Marine state:", marine_value, "(Final marine environment)\n")
