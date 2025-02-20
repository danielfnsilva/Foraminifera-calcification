# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 15:52:01 2023

@author: ldenooijer
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

## Load and Prepare the Data
## Reads the CSV file containing pH data and distance in micrometers (see example Ammonia_30degress.csv)
# Load Data
pH_data = pd.read_csv(".csv", 
                      sep=";", header=0, skipinitialspace=True, encoding='unicode_escape')

# Extract Independent (x) and Dependent (y) Variables
xdata = pH_data[["distance"]]  # Extracts the distance column
xdata = xdata.to_numpy()  # Converts the DataFrame into a NumPy array
xdata = xdata[:, 0]  # Ensures it's a 1D array

# Convert pH to Proton Concentration
ydata = 10**-(pH_data[["pH"]])  # Converts pH to H+ concentration using 10^(-pH)
ydata = ydata.to_numpy()
ydata = ydata[:, 0]

# Set the figure size for high-quality output (in inches)
plt.figure(figsize=(10, 6), dpi=300)

# Create the scatter plot
plt.scatter(xdata, ydata, s=30, c=ydata, cmap='viridis')

# Set labels and title with improved font size
plt.xlabel('Distance [m]', fontsize=12)
plt.ylabel('Proton Concentration [mol]', fontsize=12)
plt.title('Proton Concentration vs. Distance', fontsize=14)

# Add grid and colorbar
plt.grid(True)
plt.colorbar(label='Proton Concentration [mol]')

# Customize tick labels size
plt.tick_params(axis='both', which='major', labelsize=10)

# Define the Function for Curve Fitting
def func(x, alpha, k, H_inf):
    return (alpha * ((((2 * (x * k)) / 3.14) ** -0.5) * np.exp(-(x * k))) / ((k * x) ** -0.5)) + H_inf

# Fit the Model to Data
popt, pcov = curve_fit(func, xdata, ydata)  # Finds optimal parameters alpha, k, and H_inf
print(popt)  # Prints the fitted parameters

# Generate a Range of Distances for Plotting the Fitted Curve
dist = np.arange(1, 200, 1)  # Generates distance values from 1 to 199
plt.plot(dist, 
         (popt[0] * ((((2 * (dist * popt[1])) / 3.14) ** -0.5) * np.exp(-(dist * popt[1]))) / ((dist * popt[1]) ** -0.5)) + popt[2], 
         'r--')  # Plots the fitted curve as a red dashed line

# Show the plot (only once at the end)
plt.show()

# Compute Proton Concentration at Specific Distances
conc_1 = popt[0] * ((((2 * (0.05 * popt[1])) / 3.14) ** -0.5) * np.exp(-(0.05 * popt[1]))) / ((0.05 * popt[1]) ** -0.5) + popt[2]
conc_2 = popt[0] * ((((2 * (0.1 * popt[1])) / 3.14) ** -0.5) * np.exp(-(0.1 * popt[1]))) / ((0.1 * popt[1]) ** -0.5) + popt[2]

# Compute Proton Flux using Fick's First Law
D = 9.94e-5  # Diffusion constant for H+ in cm²/s at 30 degrees
J = -D * (conc_2 - conc_1) / (0.05 / 1e4)  # Converts µm to cm, calculates local flux

# Compute Proton Pumping Rate (Q)
R = 42 / 1e4  # Convert radius from µm to cm
Q = 4 * 3.14 * R**2 * J * 3600  # Converts flux from s⁻¹ to h⁻¹ assuming a spherical model
