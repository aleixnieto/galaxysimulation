import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

N = np.array([1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
execution_time = np.array([0.29000, 0.78000, 1.30000, 1.94000,
                            2.58000, 3.71000, 4.46000, 5.46000, 6.04000, 6.59000])

# Define complexity functions
def linear(N, a): return a * N
def nlogn(N, a): return a * N * np.log(N)
def quadratic(N, a): return a * N**2

# Fit models to data
popt_lin, _ = curve_fit(linear, N, execution_time)
popt_nlogn, _ = curve_fit(nlogn, N, execution_time)
popt_quad, _ = curve_fit(quadratic, N, execution_time)

# Generate fitted values
N_fit = np.linspace(min(N), max(N), 100)
time_lin = linear(N_fit, *popt_lin)
time_nlogn = nlogn(N_fit, *popt_nlogn)
time_quad = quadratic(N_fit, *popt_quad)

# Plot results
plt.figure(figsize=(8,6))
plt.scatter(N, execution_time, color="black", label="Measured data")
plt.plot(N_fit, time_lin, "r--", label=r"$O(N)$ Fit")
plt.plot(N_fit, time_nlogn, "g--", label=r"$O(N \log N)$ Fit")
plt.plot(N_fit, time_quad, "b--", label=r"$O(N^2)$ Fit")
plt.xlabel("N (Number of particles)")
plt.ylabel("Execution time (Seconds)")
plt.legend()
plt.grid()
plt.savefig("scaling_plot.png", dpi = 300)
