import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

DEGREES = [6, 9, 12, 15, 18, 24]

CSV_FILES = {
    6:  "degree6_sin.csv",
    9:  "degree9_sin.csv",
    12: "degree12_sin.csv",
    15: "degree15_sin.csv",
    18: "degree18_sin.csv",
    24: "degree24_sin.csv",
}

# use x values from one of the degree files and compute sin(5x)
ref_data = np.genfromtxt("degree6_sin.csv", delimiter=",", skip_header=1)
x_actual = ref_data[:, 0]
y_actual = np.sin(5*x_actual)

max_errors = []
rms_errors = []

for deg in DEGREES:
    data = np.genfromtxt(CSV_FILES[deg], delimiter=",", skip_header=1)
    y_approx = data[:, 2]
    err = np.abs(y_approx - y_actual)
    max_errors.append(err.max())
    rms_errors.append(np.sqrt(np.mean(err**2)))

# plot
fig, ax = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor('#F8F9FA')
ax.set_facecolor('#F8F9FA')

ax.semilogy(DEGREES, max_errors, "o-", color="red", linewidth=2, markersize=8, label="Max Error")
ax.semilogy(DEGREES, rms_errors, "s-", color="blue", linewidth=2, markersize=8, label="RMS Error")

for deg, err in zip(DEGREES, max_errors):
    ax.annotate(f"{err:.2e}", (deg, err), textcoords="offset points",
                xytext=(0, 10), ha="center", fontsize=8, color="red")

for deg, err in zip(DEGREES, rms_errors):
    ax.annotate(f"{err:.2e}", (deg, err), textcoords="offset points",
                xytext=(0, -16), ha="center", fontsize=8, color="blue")

ax.set_xlabel("Polynomial Degree", fontsize=11)
ax.set_ylabel("Error (log scale)", fontsize=11)
ax.set_title("Chebyshev — Max & RMS Error vs Polynomial Degree", fontsize=13, fontweight="bold")
ax.set_xticks(DEGREES)
ax.legend(fontsize=10)
ax.grid(True, which="both", linestyle="--", alpha=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig("chebyshev_error_vs_degree.png", dpi=150, bbox_inches="tight")
plt.show()