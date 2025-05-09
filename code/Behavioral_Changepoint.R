# Behavioral Changepoint Analysis for Otolith LA-ICPMS Data
# =========================================================
# This script implements various behavioral changepoint detection methods
# to identify transitions between different behavioral states in otolith data

# Load required libraries
library(ggplot2)           # For visualization
library(dplyr)             # For data manipulation
library(changepoint)       # Basic changepoint analysis
library(ecp)               # Energy-based changepoint detection
library(bcp)               # Bayesian changepoint detection
library(strucchange)       # Structural change models
library(mvcp)              # Multivariate changepoint detection

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

# Prepare data for behavioral analysis
# ====================================
# For behavioral analysis, we often want to look at multiple metrics
# derived from the raw data, such as:
# 1. Moving averages (to reduce noise)
# 2. First and second derivatives (rate of change)
# 3. Moving variance (stability of signal)
# 4. Complexity measures

# Calculate moving average (simple low-pass filter)
window_size <- 5
data$ma <- stats::filter(data$Iso, rep(1/window_size, window_size), sides = 2)

# Calculate first derivative (rate of change)
data$derivative <- c(NA, diff(data$Iso) / diff(data$Time))

# Calculate second derivative (acceleration of change)
data$second_deriv <- c(NA, NA, diff(data$derivative, differences = 1) / diff(data$Time[-1]))

# Calculate moving variance (as a proxy for signal stability)
calculate_moving_var <- function(x, window) {
  result <- numeric(length(x))
  for(i in 1:length(x)) {
    start_idx <- max(1, i - floor(window/2))
    end_idx <- min(length(x), i + floor(window/2))
    result[i] <- var(x[start_idx:end_idx], na.rm = TRUE)
  }
  return(result)
}

data$moving_var <- calculate_moving_var(data$Iso, window_size * 2)

# Create a summary plot of derived metrics
p_metrics <- ggplot(data, aes(x = Time)) +
  geom_line(aes(y = Iso), color = "black", alpha = 0.3) +
  geom_line(aes(y = ma), color = "blue") +
  geom_line(aes(y = derivative * 1000 + 0.7), color = "red") +  # Scaled for visibility
  geom_line(aes(y = moving_var * 10000 + 0.7), color = "green") +  # Scaled for visibility
  labs(title = "Derived Behavioral Metrics",
       subtitle = "Black: raw isotope data, Blue: moving average, Red: first derivative, Green: moving variance",
       x = "Time (Distance from core)",
       y = "Value (scaled)") +
  theme_minimal()

print(p_metrics)

# 1. Energy-Based Changepoint Detection (ecp package)
# ==================================================
# This method uses energy statistics to detect changes in distribution
# It can handle multivariate data and is non-parametric

# Prepare a cleaner dataset (remove NAs from calculations)
clean_data <- na.omit(data[,c("Time", "Iso", "ma", "derivative", "moving_var")])

# Scale the variables to comparable ranges
scaled_data <- clean_data
scaled_data[,-1] <- scale(clean_data[,-1])

# Method 1: Univariate energy-based changepoint detection
set.seed(123)
ecp_result <- e.divisive(scaled_data$Iso, sig.lvl = 0.05, min.size = 30, R = 199)

# Print results
cat("Energy-based changepoints at positions:", ecp_result$estimates, "\n")

# Convert positions to time values
ecp_times <- clean_data$Time[ecp_result$estimates]
cat("Energy-based changepoints at times:", ecp_times, "\n")

# Method 2: Multivariate energy-based changepoint detection
# Using isotope values, moving average, derivative, and variance
set.seed(123)
ecp_multi <- e.divisive(scaled_data[,c("Iso", "ma", "derivative", "moving_var")], 
                        sig.lvl = 0.05, min.size = 30, R = 199)

# Print results
cat("Multivariate energy-based changepoints at positions:", ecp_multi$estimates, "\n")

# Convert positions to time values
ecp_multi_times <- clean_data$Time[ecp_multi$estimates]
cat("Multivariate energy-based changepoints at times:", ecp_multi_times, "\n")

# Visualize energy-based changepoints
p_ecp <- ggplot(data, aes(x = Time, y = Iso)) +
  geom_point(alpha = 0.3, color = "gray40") +
  geom_vline(xintercept = ecp_times, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = ecp_multi_times, color = "red", linetype = "dotted") +
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(y = BaseSignal), color = "black", linetype = "solid")
  } +
  labs(title = "Energy-Based Changepoint Detection",
       subtitle = "Blue dashed: Univariate, Red dotted: Multivariate",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

print(p_ecp)

# 2. Bayesian Changepoint Analysis (bcp package)
# =============================================
# This approach uses Bayesian inference to detect changepoints
# and provides posterior probabilities

# Run Bayesian changepoint analysis on isotope data
bcp_result <- bcp(data$Iso)

# Extract posterior probabilities
post_probs <- bcp_result$posterior.prob

# Find changepoints using a threshold on posterior probabilities
threshold <- 0.5
bcp_cp_indices <- which(post_probs > threshold)

# If no changepoints meet the threshold, find top ones
if(length(bcp_cp_indices) < 2) {
  bcp_cp_indices <- order(post_probs, decreasing = TRUE)[1:5]
}

# Convert to time values
bcp_cp_times <- data$Time[bcp_cp_indices]

# Visualize Bayesian changepoint results
p_bcp <- ggplot() +
  geom_line(data = data.frame(Time = data$Time, Probability = post_probs),
            aes(x = Time, y = Probability)) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  geom_vline(xintercept = bcp_cp_times, linetype = "dotted", color = "blue") +
  labs(title = "Bayesian Changepoint Analysis",
       subtitle = paste("Posterior probabilities with threshold =", threshold),
       x = "Time (Distance from core)",
       y = "Posterior Probability") +
  theme_minimal()

print(p_bcp)

# Also visualize the isotope data with identified Bayesian changepoints
p_bcp_iso <- ggplot(data, aes(x = Time, y = Iso)) +
  geom_point(alpha = 0.3, color = "gray40") +
  geom_vline(xintercept = bcp_cp_times, color = "blue", linetype = "dashed") +
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(y = BaseSignal), color = "black", linetype = "solid")
  } +
  labs(title = "Bayesian Changepoint Detection on Isotope Data",
       subtitle = "Blue dashed lines show changepoints with high posterior probability",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

print(p_bcp_iso)

# 3. Structural Change Models (strucchange package)
# ================================================
# These approaches are common in econometrics for detecting
# structural changes in time series

# Method 1: CUSUM Process
# Cumulative sum of recursive residuals to detect structural changes
cusum_result <- efp(Iso ~ 1, data = data, type = "Rec-CUSUM")
cusum_bp <- breakpoints(cusum_result)

# Method 2: Breakpoints using BIC
# Find optimal breakpoints using information criteria
bp_result <- breakpoints(Iso ~ 1, data = data, h = 20)
summary(bp_result)

# Get optimal number of breaks using BIC
optimal_breaks <- bp_result$breakpoints[which.min(bp_result$RSS/log(bp_result$nobs))]
cat("Optimal number of breakpoints based on BIC:", length(optimal_breaks), "\n")

# Extract breakpoints
struc_bp_indices <- bp_result$breakpoints[1:length(optimal_breaks)]
struc_bp_times <- data$Time[struc_bp_indices]
cat("Structural change breakpoints at times:", struc_bp_times, "\n")

# Visualize structural breakpoints
p_struc <- ggplot(data, aes(x = Time, y = Iso)) +
  geom_point(alpha = 0.3, color = "gray40") +
  geom_vline(xintercept = struc_bp_times, color = "orange", linetype = "dashed") +
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(y = BaseSignal), color = "black", linetype = "solid")
  } +
  labs(title = "Structural Change Breakpoint Detection",
       subtitle = "Orange dashed lines show optimal breakpoints",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

print(p_struc)

# 4. Multivariate Changepoint Detection with MVCP
# ==============================================
# Multivariate changepoint detection using the mvcp package

# Prepare multivariate data matrix 
# Using smoothed isotope values and derivatives
data_matrix <- as.matrix(scaled_data[,c("ma", "derivative", "moving_var")])

# Run MVCP algorithm
set.seed(123)
mvcp_result <- mvcp(data_matrix, 
                    method = "BinSeg",  # Binary segmentation
                    penalty = "BIC",    # Bayesian Information Criterion
                    max.segments = 10)  # Maximum number of segments

# Extract changepoints
mvcp_cp_indices <- mvcp_result$breakpoints
mvcp_cp_times <- clean_data$Time[mvcp_cp_indices]
cat("MVCP multivariate changepoints at times:", mvcp_cp_times, "\n")

# Visualize MVCP results
p_mvcp <- ggplot(data, aes(x = Time, y = Iso)) +
  geom_point(alpha = 0.3, color = "gray40") +
  geom_vline(xintercept = mvcp_cp_times, color = "purple", linetype = "dashed") +
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(y = BaseSignal), color = "black", linetype = "solid")
  } +
  labs(title = "Multivariate Changepoint Detection (MVCP)",
       subtitle = "Purple dashed lines show multivariate changepoints",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

print(p_mvcp)

# 5. Combined Visualization of All Methods
# =======================================
# Create a comprehensive plot comparing all methods

# Combine all changepoints into one data frame for plotting
all_cp <- data.frame(
  Time = c(ecp_times[-1], ecp_multi_times[-1], bcp_cp_times, 
           struc_bp_times, mvcp_cp_times),
  Method = c(rep("Energy (Univariate)", length(ecp_times[-1])),
             rep("Energy (Multivariate)", length(ecp_multi_times[-1])),
             rep("Bayesian", length(bcp_cp_times)),
             rep("Structural", length(struc_bp_times)),
             rep("MVCP", length(mvcp_cp_times)))
)

# Add true state changes if available
if ("TrueState" %in% names(data)) {
  true_change_indices <- which(diff(floor(data$TrueState)) != 0)
  true_change_times <- data$Time[true_change_indices + 1]
  
  all_cp <- rbind(
    all_cp,
    data.frame(
      Time = true_change_times,
      Method = rep("True States", length(true_change_times))
    )
  )
}

# Create combined visualization
p_combined <- ggplot() +
  geom_point(data = data, aes(x = Time, y = Iso), alpha = 0.3, color = "gray50") +
  {if("BaseSignal" %in% names(data)) 
    geom_line(data = data, aes(x = Time, y = BaseSignal), color = "black", linetype = "solid")
  } +
  geom_vline(data = all_cp, aes(xintercept = Time, color = Method), 
             linetype = "dashed") +
  scale_color_brewer(name = "Method", palette = "Set1") +
  labs(title = "Comparison of All Behavioral Changepoint Methods",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_combined)

# 6. Accuracy Assessment
# =====================
# Calculate how well each method performed compared to true states
# (if available in the simulated data)

if ("TrueState" %in% names(data)) {
  # Function to calculate accuracy metrics for a set of changepoints
  assess_accuracy <- function(detected_times, true_times, time_tolerance) {
    # For each true changepoint, find the closest detected changepoint
    min_distances <- sapply(true_times, function(t) 
      min(abs(detected_times - t)))
    
    # Calculate true positives, false positives, false negatives
    true_positives <- sum(min_distances <= time_tolerance)
    false_negatives <- sum(min_distances > time_tolerance)
    
    # For each detected changepoint, find the closest true changepoint
    reverse_distances <- sapply(detected_times, function(t) 
      min(abs(t - true_times)))
    false_positives <- sum(reverse_distances > time_tolerance)
    
    # Calculate F1 score
    precision <- ifelse(true_positives + false_positives > 0, 
                        true_positives / (true_positives + false_positives), 0)
    recall <- true_positives / (true_positives + false_negatives)
    f1 <- ifelse(precision + recall > 0, 
                 2 * precision * recall / (precision + recall), 0)
    
    # Return metrics
    return(c(
      "True Positives" = true_positives,
      "False Positives" = false_positives,
      "False Negatives" = false_negatives,
      "Precision" = precision,
      "Recall" = recall,
      "F1 Score" = f1,
      "Mean Distance" = mean(min_distances)
    ))
  }
  
  # Define tolerance for matching (e.g., 5% of time range)
  time_range <- max(data$Time) - min(data$Time)
  tolerance <- time_range * 0.05
  
  # Calculate accuracy for each method
  accuracy_results <- data.frame(
    Method = c("Energy (Univariate)", "Energy (Multivariate)", 
               "Bayesian", "Structural", "MVCP"),
    stringsAsFactors = FALSE
  )
  
  # Add accuracy metrics for each method
  accuracy_results$TP <- NA
  accuracy_results$FP <- NA
  accuracy_results$FN <- NA
  accuracy_results$Precision <- NA
  accuracy_results$Recall <- NA
  accuracy_results$F1 <- NA
  accuracy_results$MeanDist <- NA
  
  # Energy (Univariate)
  metrics <- assess_accuracy(ecp_times[-1], true_change_times, tolerance)
  accuracy_results[1, 2:8] <- metrics
  
  # Energy (Multivariate)
  metrics <- assess_accuracy(ecp_multi_times[-1], true_change_times, tolerance)
  accuracy_results[2, 2:8] <- metrics
  
  # Bayesian
  metrics <- assess_accuracy(bcp_cp_times, true_change_times, tolerance)
  accuracy_results[3, 2:8] <- metrics
  
  # Structural
  metrics <- assess_accuracy(struc_bp_times, true_change_times, tolerance)
  accuracy_results[4, 2:8] <- metrics
  
  # MVCP
  metrics <- assess_accuracy(mvcp_cp_times, true_change_times, tolerance)
  accuracy_results[5, 2:8] <- metrics
  
  # Print accuracy results
  print("Accuracy assessment of changepoint methods compared to true states:")
  print(accuracy_results)
  
  # Visualize F1 scores
  p_f1 <- ggplot(accuracy_results, aes(x = Method, y = F1, fill = Method)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(F1, 2)), vjust = -0.5) +
    ylim(0, 1) +
    labs(title = "F1 Scores of Changepoint Detection Methods",
         subtitle = "Higher F1 scores indicate better alignment with true state changes",
         x = "Method",
         y = "F1 Score") +
    theme_minimal() +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_f1)
}

# 7. Advanced: Heterogeneous Auto-Segmentation (HDBSCAN + PCA)
# ===========================================================
# This approach combines dimension reduction and density-based clustering
# to identify segments with similar behavioral characteristics

# Load additional required libraries
library(dbscan)  # For density-based clustering
library(factoextra)  # For PCA visualization

# Create a feature matrix for clustering
# Using multiple derivatives and variance measures at different scales
features <- data.frame(
  isotope = data$Iso,
  ma_small = stats::filter(data$Iso, rep(1/3, 3), sides = 2),
  ma_large = stats::filter(data$Iso, rep(1/11, 11), sides = 2),
  deriv_small = c(NA, NA, diff(data$Iso, differences = 1, lag = 2) / diff(data$Time, lag = 2)),
  deriv_large = c(rep(NA, 5), diff(data$Iso, differences = 1, lag = 5) / diff(data$Time, lag = 5)),
  var_small = calculate_moving_var(data$Iso, 7),
  var_large = calculate_moving_var(data$Iso, 15)
)

# Remove rows with NAs
features_clean <- na.omit(features)
times_clean <- data$Time[complete.cases(features)]

# Scale the features
features_scaled <- scale(features_clean)

# Perform PCA for dimension reduction
pca_result <- prcomp(features_scaled, center = TRUE, scale. = TRUE)
summary(pca_result)

# Keep components that explain at least 90% of variance
cumulative_variance <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
n_components <- which(cumulative_variance >= 0.9)[1]
cat("Using", n_components, "principal components explaining", 
    round(cumulative_variance[n_components] * 100, 1), "% of variance\n")

# Extract PCA scores
pca_scores <- pca_result$x[, 1:n_components]

# Apply HDBSCAN clustering
# minPts controls the minimum cluster size
hdbscan_result <- hdbscan(pca_scores, minPts = 30)

# Get cluster assignments
clusters <- hdbscan_result$cluster

# Find changepoints as transitions between clusters
cluster_changes <- which(diff(clusters) != 0)
cluster_change_times <- times_clean[cluster_changes + 1]

# Visualize PCA and clustering results
pca_data <- data.frame(
  PC1 = pca_scores[, 1],
  PC2 = pca_scores[, 2],
  Cluster = factor(clusters)
)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(alpha = 0.7) +
  labs(title = "PCA Visualization of Behavioral Features",
       subtitle = "Colors represent HDBSCAN clusters",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()

print(p_pca)

# Visualize isotope data with cluster-based changepoints
data$cluster <- NA
indices <- which(complete.cases(features))
data$cluster[indices] <- clusters

# Create a clean plotting data frame
plot_data <- data.frame(
  Time = data$Time,
  Iso = data$Iso,
  Cluster = factor(data$cluster)
)

p_clusters <- ggplot(plot_data, aes(x = Time, y = Iso, color = Cluster)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = cluster_change_times, linetype = "dashed", color = "black") +
  {if("BaseSignal" %in% names(data)) 
    geom_line(aes(y = BaseSignal, color = NULL), color = "black", linetype = "solid")
  } +
  labs(title = "Behavioral Auto-Segmentation with HDBSCAN",
       subtitle = "Points colored by cluster, black dashed lines show transitions",
       x = "Time (Distance from core)",
       y = "87Sr/86Sr") +
  theme_minimal()

print(p_clusters)

# 8. Final Recommendations and Most Effective Method
# ================================================
# Based on comparison with true states (if available)

# Print final summary and recommendations
cat("\n===== BEHAVIORAL CHANGEPOINT ANALYSIS SUMMARY =====\n")

if ("TrueState" %in% names(data)) {
  # Find the best method based on F1 score
  best_method_idx <- which.max(accuracy_results$F1)
  best_method <- accuracy_results$Method[best_method_idx]
  
  cat("Best performing method based on F1 score:", best_method, "\n")
  cat("F1 score:", round(accuracy_results$F1[best_method_idx], 3), "\n\n")
  
  # Get changepoint times from the best method
  best_times <- switch(best_method,
                       "Energy (Univariate)" = ecp_times[-1],
                       "Energy (Multivariate)" = ecp_multi_times[-1],
                       "Bayesian" = bcp_cp_times,
                       "Structural" = struc_bp_times,
                       "MVCP" = mvcp_cp_times)
  
  # Create final recommended visualization
  p_final <- ggplot(data, aes(x = Time, y = Iso)) +
    geom_point(alpha = 0.3, color = "gray40") +
    geom_vline(xintercept = best_times, color = "blue", linetype = "dashed") +
    geom_vline(xintercept = true_change_times, color = "red", linetype = "dotted") +
    {if("BaseSignal" %in% names(data)) 
      geom_line(aes(y = BaseSignal), color = "black", linetype = "solid")
    } +
    labs(title = paste("Recommended Changepoint Detection:", best_method),
         subtitle = "Blue dashed: Detected changepoints, Red dotted: True state changes",
         x = "Time (Distance from core)",
         y = "87Sr/86Sr") +
    theme_minimal()
  
  print(p_final)
  
  # Recommendations
  cat("Recommendations:\n")
  cat("1. The", best_method, "method performed best for this dataset.\n")
  cat("2. ", accuracy_results$TP[best_method_idx], "of", 
      length(true_change_times), "true changepoints were correctly identified.\n")
  cat("3. There were", accuracy_results$FP[best_method_idx], 
      "false positive detections.\n")
} else {
  cat("True states not available for comparison.\n")
  cat("Recommendations based on method properties:\n")
  cat("1. Energy-based multivariate detection is robust for multiple metrics.\n")
  cat("2. Bayesian changepoint analysis provides uncertainty estimates.\n")
  cat("3. Structural change models work well for linear trends.\n")
  cat("4. HDBSCAN clustering captures complex behavioral patterns.\n")
  cat("5. Consider validating results with domain knowledge.\n")
}

# Optional: Save key plots
# ggsave("behavioral_cp_comparison.png", p_combined, width = 12, height = 8)
# ggsave("behavioral_cp_accuracy.png", p_f1, width = 8, height = 6)
# ggsave("behavioral_cp_clusters.png", p_clusters, width = 10, height = 6)