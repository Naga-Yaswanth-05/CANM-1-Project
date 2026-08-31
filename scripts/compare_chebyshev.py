import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# load the actual data
actual = np.genfromtxt("../data/data.csv", delimiter=",", skip_header=1)
x_actual = actual[:, 0]
y_actual = actual[:, 1]

# load the three approximation files
data1 = np.genfromtxt("../src/cheby_qr.csv", delimiter=",", skip_header=0)
data2 = np.genfromtxt("../src/cheby_normal.csv", delimiter=",", skip_header=0)
data3 = np.genfromtxt("../src/cheby_gradient.csv", delimiter=",", skip_header=0)

x1 = data1[:, 0]
y1 = data1[:, 2]

x2 = data2[:, 0]
y2 = data2[:, 2]

x3 = data3[:, 0]
y3 = data3[:, 2]

# calculate absolute errors
err1 = np.abs(y1 - y_actual)
err2 = np.abs(y2 - y_actual)
err3 = np.abs(y3 - y_actual)

# calculate rms errors
rms1 = np.sqrt(np.mean(err1**2))
rms2 = np.sqrt(np.mean(err2**2))
rms3 = np.sqrt(np.mean(err3**2))

# bar graph
methods = ['QR', 'Normal Equations', 'Gradient Descent']
rms_values = [rms1, rms2, rms3]
colors = ['#4C72B0', '#DD8452', '#55A868']

fig, ax = plt.subplots(figsize=(9, 6))
bars = ax.bar(methods, rms_values, color=colors, width=0.5, edgecolor='black', linewidth=0.8)

# write value on top of each bar
for bar, val in zip(bars, rms_values):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + (max(rms_values) - min(rms_values)) * 0.01,
            f'{val:.6e}', ha='center', va='bottom', fontsize=11, fontweight='bold')

# zoom in by setting y limits close to the actual values
padding = (max(rms_values) - min(rms_values)) * 0.3
ax.set_ylim(min(rms_values) - padding, max(rms_values) + padding * 4)

# format y axis to show enough precision
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.2e}'))

ax.set_ylabel('RMS Error', fontsize=12)
ax.set_title('RMS Error of Chebyshev Least Squares Approximations', fontsize=13, fontweight='bold')
ax.grid(True, axis='y', linestyle='--', alpha=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

plt.savefig("error_plot.png", dpi=150)
plt.show()