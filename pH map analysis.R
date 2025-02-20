# Load libraries

library(tiff)
library(zoo)
library(ggplot2)
library(scales)
library(readxl)
library(reshape2)
library(gridExtra)
library(imager)

# Load image
img <- readTIFF("C:/Users/xxxxxxx/Ammonia_30degress.tif")

# Replace 0 values with 0.01 in both channels (ch2 and ch3)
img[,,2][img[,,2] == 0] <- 0.01
img[,,3][img[,,3] == 0] <- 0.01

#### Check the intensity signal to remove possible outliers (low intensity values)
{
  
  # Sum the signal from both channels (combined intensity signal)
  combined_signal <- img[,,2] + img[,,3]
  
  # Calculate the total number of pixels
  total_pixels <- length(combined_signal)
  
  # Prepare data for plotting
  combined_signal_df <- as.data.frame(as.table(combined_signal))
  colnames(combined_signal_df) <- c("X", "Y", "Intensity")
  
  # Plot the combined intensity signal as a heatmap
  
  ggplot(combined_signal_df, aes(x = X, y = Y, fill = Intensity)) +
    geom_raster() +
    scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"), na.value = "transparent") +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    ) +
    labs(fill = "Intensity")
  
  # Convert the combined signal matrix to a vector for plotting
  combined_signal_vector <- as.vector(combined_signal)
  
  # Create a density plot of the intensity distribution
  
  ggplot(data.frame(Intensity = combined_signal_vector), aes(x = Intensity)) +
    geom_density(fill = "green", alpha = 0.5) +
    theme_minimal() +
    labs(title = "Intensity Distribution Density Plot", x = "Intensity", y = "Density")
  
  # Create a density object for the combined signal
  density_signal <- density(combined_signal_vector)
  
  # Find the intensity value at the peak of the density plot
  peak_intensity <- density_signal$x[which.max(density_signal$y)]
  
  # Print the peak intensity
  cat("Peak intensity value: ", peak_intensity, "\n")
  
  # Count the number of pixels where the combined intensity signal < 255
  removed_pixels <- sum(combined_signal < 0.08, na.rm = TRUE)
  
  # Set pixels to NA where combined signal is less than 255
  img[,,2][combined_signal < 0.08] <- NA
  img[,,3][combined_signal < 0.08] <- NA
  
  # Now handle the image arrays:
  img[,,2][is.na(combined_signal)] <- NA
  img[,,3][is.na(combined_signal)] <- NA
}

# Calculate the ratio
ratio <- img[,,2] / (img[,,2] + img[,,3])

# Normalize (0-1) image data using the calibration range
#Calibration data
{calibration_data <- data.frame(
  pH = c(3.05, 3.05, 3.05, 4.00, 4.00, 4.96, 4.96, 4.96, 6.09, 6.09, 6.09, 
         6.51, 6.51, 6.97, 6.97, 6.97, 7.33, 7.33, 7.33, 7.63, 7.63, 7.63, 
         7.94, 7.94, 7.94, 8.17, 8.17, 8.17, 8.43, 8.43, 8.43, 8.97, 8.97, 
         8.97, 9.97, 9.97, 9.97),
  ratio = c(0.007, 0.007, 0.007, 0.007, 0.007, 0.010, 0.010, 0.010, 0.043, 
            0.040, 0.041, 0.077, 0.079, 0.155, 0.154, 0.157, 0.230, 0.229, 
            0.230, 0.300, 0.299, 0.296, 0.350, 0.349, 0.352, 0.381, 0.380, 
            0.379, 0.401, 0.400, 0.399, 0.415, 0.419, 0.422, 0.426, 0.424, 
            0.425))
}

# Define calibration range based on known calibration data
calibration_min <- min(calibration_data$ratio, na.rm = TRUE)
calibration_max <- max(calibration_data$ratio, na.rm = TRUE)

# Normalize image data using calibration range
normalized_image_ratio <- (ratio - calibration_min) / (calibration_max - calibration_min)
normalized_image_ratio[normalized_image_ratio > 1] <- 1
normalized_image_ratio[normalized_image_ratio < 0.01] <- 0.01

# Calculate pH image using the inversed Boltzmann function (Transform ratio into pH)
pH_image <- 7.262013 - 0.462603 * log((1.007084 - 0.002912) / (normalized_image_ratio - 0.002912) - 1)

#Plot the pH image

{
  # Convert matrix to data frame for plotting
  pH_df <- as.data.frame(as.table(pH_image))
  colnames(pH_df) <- c("X", "Y", "pH")
  aspect_ratio <- ncol(pH_image) / nrow(pH_image)
  
  # Custom pH levels and corresponding colours
  pH_levels <- c(5.0, 6.5, 6.9, 7.1, 7.2, 7.3, 7.4, 7.5, 7.7, 7.9, 8.1, 9, 10.0)
  custom_colors <- c("#30123BFF", "#3B1893FF", "#3F51D3FF", "#4681F7FF", "#28BBECFF", 
                     "#1AE4B6FF", "#56FA75FF", "#A3FA3DFF", "#E5D731FF", "#FCA738FF", 
                     "#E75418FF", "#BA1E02FF", "#931003FF", "#7A0403FF")
  
  # Create the pH plot using the original pH values
  pH_plot <- ggplot(pH_df, aes(x = X, y = Y, fill = pH)) +
    geom_raster() +
    scale_fill_gradientn(
      colours = custom_colors,  
      values = rescale(pH_levels, to = c(0, 1)),  # Scale colors from pH 5 to 10
      limits = c(5, 10),  # Ensure pH values stay in the correct range
      na.value = "transparent"
    ) +
    coord_fixed(ratio = aspect_ratio) +
    theme_void() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none")
  
  pH_plot
}

# Plot the Red (R) Channel

{#Extract TL image from the RGB file Red channel (R)

R_channel <- img[,,1]

# Define the dimensions for proper overlay
img_width <- ncol(pH_image)
img_height <- nrow(pH_image)

# Prepare R channel data for plotting
R_df <- as.data.frame(as.table(R_channel))
colnames(R_df) <- c("X", "Y", "R_value")

# Plot the Red (R) Channel
R_channel_plot <- ggplot(R_df, aes(x = X, y = Y, fill = R_value)) +
  geom_raster() +
  scale_fill_gradient(low = "black", high = "white", na.value = "transparent") +
  coord_fixed(ratio = aspect_ratio) +
  theme_void() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "pt"))+
  coord_fixed(ratio = 1, xlim = c(0, img_width), ylim = c(0, img_height), expand = FALSE)
R_channel_plot

}

# Create a custom color bar and normalized it

{
pH_levels <- c(5.0, 6.5, 6.9, 7.1, 7.2, 7.3, 7.4, 7.5, 7.7, 7.9, 8.1, 9, 10.0)
custom_colors <- c("#30123BFF", "#3B1893FF", "#3F51D3FF", "#4681F7FF", "#28BBECFF", 
                   "#1AE4B6FF", "#56FA75FF", "#A3FA3DFF", "#E5D731FF", "#FCA738FF", 
                   "#E75418FF", "#BA1E02FF", "#931003FF", "#7A0403FF")

# Normalize pH levels using the normal Boltzmann function
normalize_pH <- function(pH) {
  0.002912 + (1.007084 - 0.002912) / (1 + exp((pH - 7.262013) / -0.462603))
}

# Apply normalization to pH levels and rescale between 0 and 1
custom_values <- sapply(pH_levels, normalize_pH)
custom_values <- rescale(custom_values, to = c(0, 1))  # Rescale to range between 0 and 1

# Define which pH values should be displayed on the scale bar
displayed_pH_levels <- c(5.0, 6.0, 6.5, 7.0, 7.5, 8, 8.5, 10.0)

# Get corresponding normalized values for these pH levels
displayed_values <- rescale(sapply(displayed_pH_levels, normalize_pH), to = c(0, 1))

# Create the scale bar with selected labels

scale_bar <- ggplot(data.frame(x = 1, y = seq(0, 1, length.out = 500), pH = seq(5, 10, length.out = 500)),
                    aes(x = x, y = y, fill = pH)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = custom_colors,
    values = custom_values,       # Proportional colour distribution
    limits = c(5, 10)
  ) +
  scale_y_continuous(
    breaks = displayed_values,   # Only selected pH values will be labeled
    labels = sprintf("%.1f", displayed_pH_levels),  # Format labels with one decimal place
    expand = c(0, 0),
    position = "right"  # pH values on the right side
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 7, colour = "gray50",hjust = 0, margin = margin(l = 2, r = 5)),
    axis.ticks.y = element_line(size = 0.5, colour = "black"),
    axis.ticks.length.y = unit(-4, "pt"),
    axis.line.y = element_line(size = 0.5, colour = "black")  # Thicker border around the scale bar
  )
}

# combine plots

{final_plot <- R_channel_plot +
  annotation_custom(
    grob = ggplotGrob(pH_plot),
    xmin = 0, xmax = img_width,
    ymin = 0, ymax = img_height) +
  annotation_custom(
    grob = ggplotGrob(scale_bar),
    xmin = ncol(pH_image) * 1.05, xmax = ncol(pH_image) * 1.23,  # Adjust width based on image dimensions
    ymin = 50, ymax = 1000  # Adjust height to control the range of the scale bar
  )
}

# Save the final plot to a file

ggsave(filename = "pH map.png",
       plot = final_plot,
       width = 140, height = 50,
       units = "mm",
       dpi = 1000)