# Load required libraries
library(mgcv)
library(ggplot2)

# Function to run simplified GAM analysis and create facet plots by k value for each spline type
plot_gam_splines_by_method <- function(data, dataset_name) {
  # Create directories for results
  dir.create("Figures", showWarnings = FALSE)
  dir.create("Figures/GAM", showWarnings = FALSE)
  dir.create(paste0("Figures/GAM/", dataset_name), showWarnings = FALSE)
  
  # Define k values to test
  k_values <- c(10, 30, 75, 150)
  
  # Define spline types to test
  spline_types <- c("tp", "cr", "ps", "ad")
  
  # Nice labels for the spline types
  spline_labels <- c(
    tp = "Thin Plate Splines",
    cr = "Cubic Regression Splines",
    ps = "P-Splines",
    ad = "Adaptive Splines"
  )
  
  # Process each spline type
  for (spline in spline_types) {
    # Dataframe to store predictions for this spline type
    spline_predictions <- data.frame()
    
    # Process each k value for this spline type
    for (k in k_values) {
      # Print status
      cat("Dataset:", dataset_name, "- Fitting model with spline type:", spline, "and k =", k, "\n")
      
      # Try to fit the model
      tryCatch({
        # Fit the model based on spline type
        if (spline == "ad") {
          model <- gam(Iso ~ s(Time, k=k, bs=spline, m=2), data = data, method = "REML")
        } else {
          model <- gam(Iso ~ s(Time, k=k, bs=spline), data = data, method = "REML")
        }
        
        # Get predictions
        preds <- predict(model, se.fit = TRUE)
        
        # Calculate effective degrees of freedom (edf)
        edf <- sum(model$edf)
        
        # Add to predictions dataframe
        spline_predictions <- rbind(spline_predictions, data.frame(
          Time = data$Time,
          fitted = preds$fit,
          se = preds$se.fit,
          lower_ci = preds$fit - 1.96 * preds$se.fit,
          upper_ci = preds$fit + 1.96 * preds$se.fit,
          k = factor(k),
          edf = round(edf, 1)
        ))
      }, error = function(e) {
        cat("Error fitting model with spline type:", spline, "and k =", k, "\n")
        cat("Error message:", conditionMessage(e), "\n")
      })
    }
    
    # If we have predictions for this spline type, create the facet plot
    if (nrow(spline_predictions) > 0) {
      # Create the facet plot for this spline type
      p <- ggplot() +
        geom_point(data = data, aes(x = Time, y = Iso), alpha = 0.2, color = "gray50", size = 0.7) +
        geom_line(data = spline_predictions, aes(x = Time, y = fitted), color = "blue", size = 1) +
        facet_wrap(~ k, ncol = 2, labeller = labeller(k = function(k) paste0("k = ", k, " (edf = ", 
                                                                             subset(spline_predictions, k == k)$edf[1], ")"))) +
        labs(title = paste(dataset_name, "-", spline_labels[spline]),
             subtitle = "Panels show different k values (with effective degrees of freedom)",
             x = "Time", y = "Isotope Ratio") +
        theme_minimal() +
        theme(
          strip.background = element_rect(fill = "lightblue", color = "black"),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10)
        )
      
      # Add the true signal line if available (for simulated data)
      if ("BaseSignal" %in% colnames(data)) {
        p <- p + geom_line(data = data, aes(x = Time, y = BaseSignal), 
                           color = "red", linetype = "dashed", size = 1)
      }
      
      # Save the plot in the Figures/GAM directory
      output_file <- paste0("Figures/GAM/", dataset_name, "/", spline, "_faceted_by_k.png")
      ggsave(output_file, p, width = 10, height = 8, dpi = 300)
      
      cat("Created plot for", spline_labels[spline], "saved to", output_file, "\n")
    } else {
      cat("No predictions generated for", spline, "on", dataset_name, "\n")
    }
  }
  
  cat("Completed analysis for", dataset_name, "\n")
}

# Function to find and load data files safely
safe_read_csv <- function(file_patterns, name) {
  for (pattern in file_patterns) {
    # Find files matching the pattern
    files <- list.files(pattern = pattern, recursive = TRUE, full.names = TRUE)
    
    if (length(files) > 0) {
      cat("Found", name, "file:", files[1], "\n")
      return(read.csv(files[1]))
    }
  }
  
  cat("Could not find any", name, "files matching patterns:", paste(file_patterns, collapse=", "), "\n")
  return(NULL)
}

# Try to read the simulated data
simulated_data <- safe_read_csv(
  c("simulated_isotope_data\\.csv$", 
    "sim.*\\.csv$"),
  "simulated"
)

if (!is.null(simulated_data)) {
  plot_gam_splines_by_method(simulated_data, "simulated")
}

# Try to read the real data
real_data <- safe_read_csv(
  c("2016_yk_006\\.csv$"),
  "real otolith"
)

if (!is.null(real_data)) {
  plot_gam_splines_by_method(real_data, "2016_yk_006")
}

cat("Analysis complete. Facet plots saved in Figures/GAM directory.\n")
