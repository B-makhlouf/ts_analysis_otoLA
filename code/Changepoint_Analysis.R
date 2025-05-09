# Changepoint Analysis for Otolith LA-ICPMS Data
# =============================================
# This script implements changepoint detection methods to identify 
# transitions between different states in otolith isotope data

# Load required libraries
library(changepoint)  # For changepoint detection
library(ggplot2)      # For visualization
library(dplyr)        # For data manipulation

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

# 1. Basic Changepoint Analysis
# -----------------------------
# Convert isotope data to a numeric vector
iso_values <- data$Iso

# Method 1: Changes in Mean (AMOC - At Most One Change)
# Detect a single changepoint in the mean
cpt_mean_single <- cpt.mean(iso_values, method = "AMOC")
cat("Single changepoint in mean at position:", cpt.mean(iso_values, method = "AMOC")@cpts, "\n")

# Method 2: Changes in Mean (BinSeg - Binary Segmentation)
# Detect multiple changepoints in the mean
# Q parameter controls the maximum number of changepoints to detect
cpt_mean_multiple <- cpt.mean(iso_values, method = "BinSeg", Q = 10)
print(cpt_mean_multiple)
cat("Multiple changepoints in mean at positions:", cpt_mean_multiple@cpts, "\n")

# Method 3: Changes in Mean and Variance
# Detect changes in both mean and variance simultaneously
cpt_meanvar <- cpt.meanvar(iso_values, method = "BinSeg", Q = 10)
cat("Changepoints in mean and variance at positions:", cpt_meanvar@cpts, "\n")

# Convert changepoint positions to time values
cp_times_mean <- data$Time[cpt_mean_multiple@cpts]
cp_times_meanvar <- data$Time[cpt_meanvar@cpts]

# Visualize the detected changepoints
p_changepoints <- ggplot(data, aes(x = Time, y = Iso)) +
  geom_point(alpha = 0.3, color = "gray40") +
  geom_vline(xintercept = cp_times_mean, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = cp_times_meanvar, color = "red", linetype = "dotted") +
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(y = BaseSignal), color = "black", linetype = "solid")
  } +
  labs(title = "Changepoint Detection Results",
       subtitle = "Blue dashed: Mean changes, Red dotted: Mean+Variance changes",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

print(p_changepoints)

# 2. PELT Algorithm for Optimal Segmentation
# ------------------------------------------
# PELT (Pruned Exact Linear Time) algorithm is more efficient for finding
# the optimal segmentation than BinSeg for many changepoints

# Using PELT to detect changes in mean
cpt_pelt_mean <- cpt.mean(iso_values, method = "PELT")
cat("PELT changepoints in mean at positions:", cpt_pelt_mean@cpts, "\n")

# Using PELT to detect changes in mean and variance
cpt_pelt_meanvar <- cpt.meanvar(iso_values, method = "PELT")
cat("PELT changepoints in mean and variance at positions:", cpt_pelt_meanvar@cpts, "\n")

# Convert PELT changepoint positions to time values
cp_times_pelt_mean <- data$Time[cpt_pelt_mean@cpts]
cp_times_pelt_meanvar <- data$Time[cpt_pelt_meanvar@cpts]

# Visualize the PELT changepoints
p_pelt <- ggplot(data, aes(x = Time, y = Iso)) +
  geom_point(alpha = 0.3, color = "gray40") +
  geom_vline(xintercept = cp_times_pelt_mean, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = cp_times_pelt_meanvar, color = "red", linetype = "dotted") +
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(y = BaseSignal), color = "black", linetype = "solid")
  } +
  labs(title = "PELT Changepoint Detection Results",
       subtitle = "Blue dashed: Mean changes, Red dotted: Mean+Variance changes",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

print(p_pelt)

# 3. Create Segments Based on Changepoints
# ----------------------------------------
# We'll use the PELT mean changepoints for segmentation
changepoints <- c(0, cpt_pelt_mean@cpts)
segments <- data.frame()

# Calculate segment statistics
for (i in 1:(length(changepoints))) {
  if (i < length(changepoints)) {
    segment_start <- changepoints[i] + 1
    segment_end <- changepoints[i + 1]
    
    # Extract segment data
    segment_data <- iso_values[segment_start:segment_end]
    
    # Calculate statistics
    segment_row <- data.frame(
      SegmentID = i,
      StartIdx = segment_start,
      EndIdx = segment_end,
      StartTime = data$Time[ifelse(segment_start == 0, 1, segment_start)],
      EndTime = data$Time[segment_end],
      Length = segment_end - segment_start + 1,
      Mean = mean(segment_data),
      SD = sd(segment_data),
      Min = min(segment_data),
      Max = max(segment_data)
    )
    
    segments <- rbind(segments, segment_row)
  }
}

# Print segment information
print("Segment information based on PELT mean changepoints:")
print(segments)

# Assign segment IDs to original data
data$segment <- NA
for (i in 1:nrow(segments)) {
  idx_start <- segments$StartIdx[i]
  idx_end <- segments$EndIdx[i]
  data$segment[idx_start:idx_end] <- i
}

# Visualize segmentation
p_segments <- ggplot(data, aes(x = Time, y = Iso, color = factor(segment))) +
  geom_point(alpha = 0.6) +
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(y = BaseSignal, color = NULL), color = "black", linetype = "dashed")
  } +
  scale_color_brewer(name = "Segment", palette = "Set1") +
  labs(title = "Data Segmentation Based on PELT Changepoints",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

print(p_segments)

# 4. Advanced: Changepoint Detection with Penalty Selection
# --------------------------------------------------------
# The penalty parameter controls the trade-off between model fit and complexity
# (how many changepoints are detected)

# Try a range of penalty values to see how they affect changepoint detection
penalty_values <- c("None", "SIC", "BIC", "AIC", "MBIC", "Manual")
penalty_manual <- c("Manual1" = 0.1, "Manual2" = 0.5, "Manual3" = 1, "Manual4" = 5)

# Initialize list to store results
cpt_results <- list()

# Run changepoint detection with different penalties
for (pen in penalty_values) {
  if (pen == "Manual") {
    for (i in seq_along(penalty_manual)) {
      pen_name <- names(penalty_manual)[i]
      cpt_results[[pen_name]] <- cpt.mean(iso_values, method = "PELT", 
                                          penalty = "Manual", 
                                          pen.value = penalty_manual[i])
    }
  } else {
    cpt_results[[pen]] <- cpt.mean(iso_values, method = "PELT", penalty = pen)
  }
}

# Compile results
penalty_summary <- data.frame(
  Penalty = character(),
  NumChangepoints = integer(),
  stringsAsFactors = FALSE
)

for (name in names(cpt_results)) {
  penalty_summary <- rbind(
    penalty_summary,
    data.frame(
      Penalty = name,
      NumChangepoints = length(cpt_results[[name]]@cpts),
      stringsAsFactors = FALSE
    )
  )
}

# Display penalty comparison
print("Comparison of different penalty methods:")
print(penalty_summary)

# Choose one penalty method for visualization (e.g., BIC)
selected_penalty <- "BIC"
selected_cpts <- cpt_results[[selected_penalty]]@cpts
selected_times <- data$Time[selected_cpts]

# Visualize with selected penalty
p_penalty <- ggplot(data, aes(x = Time, y = Iso)) +
  geom_point(alpha = 0.3, color = "gray40") +
  geom_vline(xintercept = selected_times, color = "blue", linetype = "dashed") +
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(y = BaseSignal), color = "black", linetype = "solid")
  } +
  labs(title = paste("Changepoint Detection with", selected_penalty, "Penalty"),
       subtitle = paste(length(selected_cpts), "changepoints detected"),
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

print(p_penalty)

# 5. Compare with true state changes (if available)
# ------------------------------------------------
if ("TrueState" %in% names(data)) {
  # Identify indices where true state changes
  true_change_indices <- which(diff(floor(data$TrueState)) != 0)
  true_change_times <- data$Time[true_change_indices + 1]  # +1 because diff reduces length by 1
  
  cat("True state changes at times:", true_change_times, "\n")
  
  # Compare true and detected changepoints
  p_compare <- ggplot(data, aes(x = Time, y = Iso)) +
    geom_point(alpha = 0.3, color = "gray40") +
    geom_vline(xintercept = selected_times, color = "blue", linetype = "dashed") +
    geom_vline(xintercept = true_change_times, color = "red", linetype = "dotted") +
    {if("BaseSignal" %in% names(data)) 
      geom_line(aes(y = BaseSignal), color = "black", linetype = "solid")
    } +
    labs(title = "Comparison of Detected vs. True Changepoints",
         subtitle = "Blue dashed: Detected changepoints, Red dotted: True state changes",
         x = "Time (Distance from core)",
         y = "87Sr/86Sr") +
    theme_minimal()
  
  print(p_compare)
  
  # Calculate accuracy of changepoint detection
  # For each true changepoint, find the closest detected changepoint
  accuracy_distance <- numeric(length(true_change_times))
  for (i in seq_along(true_change_times)) {
    closest_cp <- selected_times[which.min(abs(selected_times - true_change_times[i]))]
    accuracy_distance[i] <- abs(closest_cp - true_change_times[i])
  }
  
  cat("Mean distance between true and detected changepoints:", 
      mean(accuracy_distance), "time units\n")
  
  # Detect if we missed any true changepoints or had false positives
  missed <- sum(sapply(true_change_times, function(t) 
    min(abs(selected_times - t)) > (max(data$Time) - min(data$Time))/50))
  
  false_positives <- sum(sapply(selected_times, function(t) 
    min(abs(true_change_times - t)) > (max(data$Time) - min(data$Time))/50))
  
  cat("Missed changepoints:", missed, "\n")
  cat("False positive changepoints:", false_positives, "\n")
}

# 6. Optional: Behavioral Changepoint Analysis (not fully implemented)
# -------------------------------------------------------------------
# This is a conceptual extension of traditional changepoint analysis
# that could be implemented for more complex behavioral data.
# For a full implementation, you would need to:
# 1. Calculate additional metrics (e.g. derivatives, moving variances)
# 2. Apply multivariate changepoint detection
# 3. Use domain knowledge to interpret the biological meaning of changes

# Simple example: Calculate first differences and analyze those
data$iso_diff <- c(NA, diff(data$Iso))

# Detect changepoints in the rate of change (first differences)
iso_diff_values <- data$iso_diff[-1]  # Remove NA at start
cpt_diff <- cpt.mean(iso_diff_values, method = "PELT", penalty = "BIC")
cat("Changepoints in isotope rate of change at positions:", cpt_diff@cpts, "\n")

# Convert to time values
cp_times_diff <- data$Time[cpt_diff@cpts + 1]  # +1 for NA offset

# Visualize the results
p_diff <- ggplot() +
  geom_point(data = data[-1,], aes(x = Time, y = iso_diff), alpha = 0.3, color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = cp_times_diff, color = "purple", linetype = "dashed") +
  labs(title = "Changepoints in Rate of Isotope Change",
       x = "Time (Distance from core)",
       y = "First Difference of 87Sr/86Sr") +
  theme_minimal()

print(p_diff)

# Optional: Save key plots
# ggsave("changepoint_segmentation.png", p_segments, width = 10, height = 6)
# ggsave("changepoint_comparison.png", p_compare, width = 10, height = 6)