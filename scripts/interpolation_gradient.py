import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

#CONFIG
COMPARISON_FILE = "comparison.csv"   # columns: x, Newton, Chebyshev
GRAD_POINTS     = 600                # points for gradient error evaluation

C_NDD = "blue"
C_CHB = "red"


#Load comparison.csv 
def load_files():
    if not os.path.exists(COMPARISON_FILE):
        print(f"[ERROR] '{COMPARISON_FILE}' not found.")
        exit(1)
    df = pd.read_csv(COMPARISON_FILE)
    print(f"[✓] Loaded '{COMPARISON_FILE}' ({len(df)} points)")
    return df


#1a. Absolute Error
def plot_absolute_error(df):
    true_vals     = np.sin(5.0 * df["x"])
    newton_err    = np.abs(df["Newton"]    - true_vals)
    chebyshev_err = np.abs(df["Chebyshev"] - true_vals)

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle("Interpolation — Absolute Error  |P(x) − sin(5x)|",
                 fontsize=13, fontweight="bold")

    ax.semilogy(df["x"], newton_err,    color=C_NDD, linewidth=2,
                label="Newton Divided Difference")
    ax.semilogy(df["x"], chebyshev_err, color=C_CHB, linewidth=2,
                label="Chebyshev")
    ax.set_xlabel("x")
    ax.set_ylabel("Absolute Error (log scale)")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig("task4_1a_absolute_error.png", dpi=150, bbox_inches="tight")
    print("[✓] Saved task4_1a_absolute_error.png")
    plt.show()

    return newton_err, chebyshev_err


# ── 1b. Cumulative RMS Error 
def plot_rms_error(df, newton_err, chebyshev_err):
    n = np.arange(len(df)) + 1
    newton_rms    = np.sqrt(np.cumsum(newton_err**2)    / n)
    chebyshev_rms = np.sqrt(np.cumsum(chebyshev_err**2) / n)

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle("Interpolation — Cumulative RMS Error",
                 fontsize=13, fontweight="bold")

    ax.semilogy(df["x"], newton_rms,    color=C_NDD, linewidth=2,
                label="Newton Divided Difference")
    ax.semilogy(df["x"], chebyshev_rms, color=C_CHB, linewidth=2,
                label="Chebyshev")
    ax.set_xlabel("x")
    ax.set_ylabel("RMS Error (log scale)")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig("task4_1b_rms_error.png", dpi=150, bbox_inches="tight")
    print("[✓] Saved task4_1b_rms_error.png")
    plt.show()

    return newton_rms, chebyshev_rms


# ── 1c. Max & RMS Bar Chart
def plot_max_rms_comparison(newton_err, chebyshev_err, newton_rms, chebyshev_rms):
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.suptitle("Interpolation — Max & RMS Error Comparison",
                 fontsize=13, fontweight="bold")

    methods  = ["Newton DD", "Chebyshev"]
    max_errs = [newton_err.max(),          chebyshev_err.max()]
    rms_errs = [float(newton_rms.iloc[-1]), float(chebyshev_rms.iloc[-1])]

    x_pos = np.arange(len(methods))
    w = 0.35
    bars1 = ax.bar(x_pos - w/2, max_errs, w, label="Max Error",
                   color=[C_NDD, C_CHB], alpha=0.85)
    bars2 = ax.bar(x_pos + w/2, rms_errs, w, label="RMS Error",
                   color=[C_NDD, C_CHB], alpha=0.45, hatch="//")
    ax.set_yscale("log")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(methods)
    ax.set_ylabel("Error (log scale)")
    ax.legend()
    ax.grid(True, axis="y", linestyle="--", alpha=0.5)

    for bar in list(bars1) + list(bars2):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.1,
                f"{bar.get_height():.2e}", ha="center", va="bottom", fontsize=8)

    plt.tight_layout()
    plt.savefig("task4_1c_max_rms_comparison.png", dpi=150, bbox_inches="tight")
    print("[✓] Saved task4_1c_max_rms_comparison.png")
    plt.show()


# ── 2. Gradient Error ─────────────────────────────────────────────────────────
# Load original coefficients (c[0] first, one per line),
# differentiate analytically, evaluate at GRAD_POINTS, compare to 5cos(5x)
NEWTON_COEFF_FILE    = "newton_coeffs.csv"
CHEBYSHEV_COEFF_FILE = "chebyshev_coeffs.csv"

def load_coeffs(filename):
    if not os.path.exists(filename):
        print(f"[ERROR] '{filename}' not found.")
        exit(1)
    coeffs = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line:
                coeffs.append(float(line))
    print(f"[✓] Loaded {len(coeffs)} coefficients from '{filename}'")
    return np.array(coeffs)   # c[0], c[1], ..., c[n]

def deriv_coeffs(c):
    # P(x)  = c[0] + c[1]*x + c[2]*x^2 + ... + c[n]*x^n
    # P'(x) = c[1] + 2*c[2]*x + ... + n*c[n]*x^(n-1)
    # Derivative: multiply each coefficient by its degree, drop the constant
    return [i * c[i] for i in range(1, len(c))]


def eval_poly(c, x):
    # c is low-degree-first: c[0] + c[1]*x + c[2]*x^2 + ...
    # Use Horner's method for numerical stability and efficiency:
    # c[0] + x*(c[1] + x*(c[2] + x*(...)))
    result = 0
    for coeff in reversed(c):
        result = result * x + coeff
    return result

def plot_gradient_error():
    newton_c    = load_coeffs(NEWTON_COEFF_FILE)
    chebyshev_c = load_coeffs(CHEBYSHEV_COEFF_FILE)

    newton_dc    = deriv_coeffs(newton_c)
    chebyshev_dc = deriv_coeffs(chebyshev_c)

    x_grad    = np.linspace(-1.0, 1.0, GRAD_POINTS)
    true_grad = 5.0 * np.cos(5.0 * x_grad)

    newton_grad_err    = np.abs(eval_poly(newton_dc,    x_grad) - true_grad)
    chebyshev_grad_err = np.abs(eval_poly(chebyshev_dc, x_grad) - true_grad)

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle("Gradient Error  |P′(x) − 5cos(5x)|",
                 fontsize=13, fontweight="bold")

    ax.semilogy(x_grad, newton_grad_err,    color=C_NDD, linewidth=2,
                label="Newton Divided Difference")
    ax.semilogy(x_grad, chebyshev_grad_err, color=C_CHB, linewidth=2,
                label="Chebyshev")
    ax.set_xlabel("x")
    ax.set_ylabel("|P′(x) − 5cos(5x)| (log scale)")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig("task4_2_gradient_error.png", dpi=150, bbox_inches="tight")
    print("[✓] Saved task4_2_gradient_error.png")
    plt.show()


# ── 1d. Absolute Error at 25 Original Nodes (scatter) ────────────────────────
NODES_FILE = "comparision_25.csv"

def plot_nodes_error():
    if not os.path.exists(NODES_FILE):
        print(f"[SKIP] '{NODES_FILE}' not found — skipping nodes scatter plot.")
        return

    df25 = pd.read_csv(NODES_FILE, header=None, names=["x", "Newton", "Chebyshev"])
    print(f"[✓] Loaded '{NODES_FILE}' ({len(df25)} nodes)")

    true_vals     = np.sin(5.0 * df25["x"])
    newton_err    = np.abs(df25["Newton"]    - true_vals)
    chebyshev_err = np.abs(df25["Chebyshev"] - true_vals)

    fig, ax = plt.subplots(figsize=(12, 4))
    fig.suptitle("Absolute Error at Original 25 Data Nodes  (log scale)",
                 fontsize=13, fontweight="bold")

    ax.scatter(df25["x"], newton_err,    color=C_NDD, marker="o", s=60, zorder=3,
               label=f"Newton    max = {newton_err.max():.3e}")
    ax.scatter(df25["x"], chebyshev_err, color=C_CHB, marker="s", s=60, zorder=3,
               label=f"Chebyshev max = {chebyshev_err.max():.3e}")

    ax.set_yscale("log")
    ax.set_xlabel("x")
    ax.set_ylabel("|error|")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", alpha=0.5)

    # Summary box below plot
    n_max = newton_err.max()
    c_max = chebyshev_err.max()
    n_rms = np.sqrt((newton_err**2).mean())
    c_rms = np.sqrt((chebyshev_err**2).mean())
    summary = f"Newton -> max|err| = {n_max:.4e}   RMS = {n_rms:.4e}\nChebyshev -> max|err| = {c_max:.4e}   RMS = {c_rms:.4e}"
    fig.text(0.5, -0.04, summary, ha="center", fontsize=9,
             bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray"))

    plt.tight_layout()
    plt.savefig("task4_1d_nodes_scatter.png", dpi=150, bbox_inches="tight")
    print("[✓] Saved task4_1d_nodes_scatter.png")
    plt.show()


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    df = load_files()

    newton_err, chebyshev_err             = plot_absolute_error(df)
    newton_rms, chebyshev_rms             = plot_rms_error(df, newton_err, chebyshev_err)
    plot_max_rms_comparison(newton_err, chebyshev_err, newton_rms, chebyshev_rms)
    plot_gradient_error()
    plot_nodes_error()