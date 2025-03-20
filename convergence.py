import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define second-order reference line function
def second_order_reference(dt_values, reference_error):
    return reference_error * (dt_values / dt_values[0]) ** 2

# Define first-order reference line function
def first_order_reference(dt_values, reference_error):
    return reference_error * (dt_values / dt_values[0])

# Load modified data for Velocity Verlet and Symplectic Euler
df_verlet = pd.DataFrame({
    "dt": [1.0e-03, 1.0e-04, 1.0e-05, 1.0e-06, 1.0e-07, 1.0e-08],
    "pos_maxdiff": [0.057047645, 0.003320566, 0.000305201, 0.000029980, 0.000002723, 0.000000000]
})

df_symplectic = pd.DataFrame({
    "dt": [1.0e-03, 1.0e-04, 1.0e-05, 1.0e-06, 1.0e-07, 1.0e-08],
    "pos_maxdiff": [0.068535808, 0.005896742, 0.000590297, 0.000058509, 0.000005319, 0.000000000]
})

# Remove rows where pos_maxdiff is exactly 0
df_verlet = df_verlet[df_verlet["pos_maxdiff"] > 0]
df_symplectic = df_symplectic[df_symplectic["pos_maxdiff"] > 0]

# Compute log values
log_dt_verlet = np.log(df_verlet["dt"])
log_error_verlet = np.log(df_verlet["pos_maxdiff"])
log_dt_symplectic = np.log(df_symplectic["dt"])
log_error_symplectic = np.log(df_symplectic["pos_maxdiff"])

# Generate second-order and first-order reference lines
dt_values = np.array(df_verlet["dt"])
second_order_errors = second_order_reference(dt_values, df_verlet["pos_maxdiff"].iloc[0])
first_order_errors = first_order_reference(dt_values, df_symplectic["pos_maxdiff"].iloc[0])

# Plot results
plt.figure(figsize=(8,6))
plt.plot(log_dt_verlet, log_error_verlet, 'o-', label="Velocity Verlet (O(Δt²))")
plt.plot(log_dt_symplectic, log_error_symplectic, 's-', label="Symplectic Euler (O(Δt))")
plt.plot(log_dt_verlet, np.log(second_order_errors), '--', label="Reference: O(Δt²)")
plt.plot(log_dt_verlet, np.log(first_order_errors), '--', label="Reference: O(Δt)")

plt.xlabel("log(Δt)")
plt.ylabel("log(pos_maxdiff)")
plt.title("Convergence Study: Velocity Verlet vs. Symplectic Euler")
plt.legend()
plt.grid()
plt.savefig("convergence_plot.png", dpi = 300)

