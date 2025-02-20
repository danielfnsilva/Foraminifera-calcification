# Necessary libraries

library(readxl)
library(ggplot2)
library(drc)
library(minpack.lm)

# Upload an excel file containing the measured seawater pH and the calculated ratios. 

data <- read_excel("")
head(data)

df <- data.frame(data)

# Normalize the ratio values between 0 and 1 to fit the Boltzmann function
df$normalized_ratio <- (df$ratio - min(df$ratio)) / (max(df$ratio) - min(df$ratio))

# Get a general plot of the data
plot(df$pH, df$normalized_ratio, main = "pH vs Normalized Ratio", xlab = "pH", ylab = "Normalized Ratio", col = "red", pch = 16)
plot(df$pH, df$ratio)
# Display normalized ratio and corresponding pH in the console
print(df[, c("pH", "normalized_ratio")])

# Fit the 4PL model to estimate initial parameters
# The 4PL model gives us four parameters: 
# A (minimum value of the curve), B (maximum value), C (inflection point, or midpoint), and D (slope).
# These parameters will serve as starting estimates for the Boltzmann function, as it is based on the same four parameters.
fm1 <- nls(normalized_ratio ~ SSfpl(pH, A, B, xmid, scal), data = df)
fm1
summary(fm1)

# The 4PL equation is: A + (B - A) / (1 + (x / C)^D)

# Now, we define the Boltzmann function. The parameters estimated from the 4PL are used for this:
# A1 (upper asymptote), A2 (lower asymptote), x0 (inflection point), dx (slope factor).
# The Boltzmann equation is: A2 + (A1 - A2) / (1 + exp((x0 - x) / dx))
# This function will be fitted to the normalized ratio data.

# Define the Boltzmann function (calculate the ratio based on the pH)
boltzmann <- function(x, A1, A2, x0, dx) {
  A2 + (A1 - A2) / (1 + exp((x0 - x) / dx))  # Boltzmann function formula
}

# Fit the Boltzmann model using the 4PL estimates as starting values
# The 4PL estimates are used here because both models share similar structures, with parameters A1 and A2 corresponding 
# to the upper and lower asymptotes, and x0 representing the inflection point.
model_boltzmann <- nls(normalized_ratio ~ boltzmann(pH, A1, A2, x0, dx), data = df, 
                       start = list(A1 = 1.007084, A2 = 0.002912, x0 = 7.2620127, dx = 0.4626025))

summary(model_boltzmann)

# Predicted values from the Boltzmann model
boltzmann_predictions <- predict(model_boltzmann)

# Add the predicted values to the data frame
df$boltzmann_fit <- boltzmann_predictions

# Plot the data and the fitted Boltzmann curve
df$predicted_ratio <- predict(model_boltzmann, newdata = df)

HPTSplot_boltzmann <- ggplot(df, aes(x = pH, y = normalized_ratio)) +
  geom_line(aes(y = predicted_ratio), color = "gray", size = 0.5) + 
  geom_point(size = 1, col = "black") +
  scale_x_continuous(
    minor_breaks = seq(2, 11, by = 0.2),
    breaks = seq(2, 11, by = 1), limits = c(2, 11),
    guide = "axis_minor"
  ) +
  scale_y_continuous(minor_breaks = seq(0, 1, by = 0.02),
                     breaks = seq(0, 1, by = 0.1), limits = c(0, 1.01),
                     guide = "axis_minor") +
  annotate("text", x = 6, y = 1, label = "0.002 + (1 - 0.002)/(1 + exp((7.26 - pH) / 0.46))", size = 2) +
  theme_bw() + 
  ylab(expression(Normalized~~λ^{"488"}/λ^{"488"}+~~λ^{"405"})) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 7, face = "bold"), axis.title.x = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))

# Display the plot
HPTSplot_boltzmann

# Save the plot to a file
ggsave(filename = "HPTS calibration 20X.png",
       plot = HPTSplot_boltzmann,
       width = 90, height = 90,
       units = "mm",
       dpi = 1000)