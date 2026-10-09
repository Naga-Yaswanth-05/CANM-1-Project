"""Plot interpolation, RMS, node, and derivative errors for the active dataset."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from project_utils import (
    OUTPUT_DIR,
    active_dataset_filename,
    dataset_label,
    exact_derivative,
    exact_function,
)

GRAD_POINTS = 600
C_NDD = "blue"
C_CHB = "red"


def plot_absolute_error(df, dataset_filename):
    true_vals = exact_function(dataset_filename, df["x"].to_numpy())
    newton_err = np.abs(df["Newton"].to_numpy() - true_vals)
    chebyshev_err = np.abs(df["Chebyshev"].to_numpy() - true_vals)

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle(f"Interpolation — Absolute Error  |P(x) − {dataset_label(dataset_filename)}|",
                 fontsize=13, fontweight="bold")
    ax.semilogy(df["x"], newton_err, color=C_NDD, linewidth=2,
                label="Newton Divided Difference")
    ax.semilogy(df["x"], chebyshev_err, color=C_CHB, linewidth=2,
                label="Chebyshev")
    ax.set_xlabel("x")
    ax.set_ylabel("Absolute Error (log scale)")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    path = OUTPUT_DIR / "task4_1a_absolute_error.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[✓] Saved {path}")
    plt.show()
    return newton_err, chebyshev_err


def plot_rms_error(df, newton_err, chebyshev_err):
    n = np.arange(len(df)) + 1
    newton_rms = np.sqrt(np.cumsum(newton_err**2) / n)
    chebyshev_rms = np.sqrt(np.cumsum(chebyshev_err**2) / n)

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle("Interpolation — Cumulative RMS Error", fontsize=13, fontweight="bold")
    ax.semilogy(df["x"], newton_rms, color=C_NDD, linewidth=2,
                label="Newton Divided Difference")
    ax.semilogy(df["x"], chebyshev_rms, color=C_CHB, linewidth=2, label="Chebyshev")
    ax.set_xlabel("x")
    ax.set_ylabel("RMS Error (log scale)")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    path = OUTPUT_DIR / "task4_1b_rms_error.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[✓] Saved {path}")
    plt.show()
    return newton_rms, chebyshev_rms


def plot_max_rms_comparison(newton_err, chebyshev_err, newton_rms, chebyshev_rms):
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.suptitle("Interpolation — Max & RMS Error Comparison", fontsize=13, fontweight="bold")
    methods = ["Newton DD", "Chebyshev"]
    max_errs = [float(np.max(newton_err)), float(np.max(chebyshev_err))]
    rms_errs = [float(newton_rms[-1]), float(chebyshev_rms[-1])]
    x_pos = np.arange(len(methods))
    width = 0.35
    bars1 = ax.bar(x_pos - width / 2, max_errs, width, label="Max Error",
                   color=[C_NDD, C_CHB], alpha=0.85)
    bars2 = ax.bar(x_pos + width / 2, rms_errs, width, label="RMS Error",
                   color=[C_NDD, C_CHB], alpha=0.45, hatch="//")
    ax.set_yscale("log")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(methods)
    ax.set_ylabel("Error (log scale)")
    ax.legend()
    ax.grid(True, axis="y", linestyle="--", alpha=0.5)
    for bar in list(bars1) + list(bars2):
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, height * 1.1,
                    f"{height:.2e}", ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    path = OUTPUT_DIR / "task4_1c_max_rms_comparison.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[✓] Saved {path}")
    plt.show()


def load_coeffs(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f"Coefficient file not found: {path}. Run the C++ main program first.")
    coefficients = []
    with path.open(encoding="utf-8") as file:
        for line in file:
            value = line.strip()
            if value:
                coefficients.append(float(value))
    if not coefficients:
        raise ValueError(f"Coefficient file is empty: {path}")
    print(f"[✓] Loaded {len(coefficients)} coefficients from '{path}'")
    return np.asarray(coefficients, dtype=float)


def deriv_coeffs(coefficients):
    return np.asarray([i * coefficients[i] for i in range(1, len(coefficients))], dtype=float)


def eval_poly(coefficients, x):
    result = 0.0
    for coefficient in reversed(coefficients):
        result = result * x + coefficient
    return result


def plot_gradient_error(dataset_filename, x_min, x_max):
    newton_c = load_coeffs(OUTPUT_DIR / "polynomial_coefficients.csv")
    chebyshev_c = load_coeffs(OUTPUT_DIR / "chebyshev_coeffs.csv")
    newton_dc = deriv_coeffs(newton_c)
    chebyshev_dc = deriv_coeffs(chebyshev_c)

    x_grad = np.linspace(x_min, x_max, GRAD_POINTS)
    true_grad = exact_derivative(dataset_filename, x_grad)
    newton_grad_err = np.abs(eval_poly(newton_dc, x_grad) - true_grad)
    chebyshev_grad_err = np.abs(eval_poly(chebyshev_dc, x_grad) - true_grad)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_title(f"Derivative Error  |P′(x) − f′(x)| for {dataset_label(dataset_filename)}",
                 fontsize=13, fontweight="bold")
    ax.semilogy(x_grad, newton_grad_err, color=C_NDD, linewidth=2,
                label="Newton Divided Difference")
    ax.semilogy(x_grad, chebyshev_grad_err, color=C_CHB, linewidth=2, label="Chebyshev")
    ax.set_xlabel("x")
    ax.set_ylabel("Derivative absolute error (log scale)")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    path = OUTPUT_DIR / "task4_2_gradient_error.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[✓] Saved {path}")
    plt.show()


def plot_nodes_error(dataset_filename):
    nodes_file = OUTPUT_DIR / "comparison_nodes.csv"
    if not nodes_file.is_file():
        raise FileNotFoundError(f"Node comparison file not found: {nodes_file}. Run the C++ main program first.")
    df_nodes = pd.read_csv(nodes_file, header=None, names=["x", "Newton", "Chebyshev"])
    if df_nodes.empty or df_nodes.isna().any().any():
        raise ValueError(f"Invalid node comparison CSV: {nodes_file}")

    true_vals = exact_function(dataset_filename, df_nodes["x"].to_numpy())
    newton_err = np.abs(df_nodes["Newton"].to_numpy() - true_vals)
    chebyshev_err = np.abs(df_nodes["Chebyshev"].to_numpy() - true_vals)

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.set_title(f"Absolute Error at Original {len(df_nodes)} Data Nodes (log scale)",
                 fontsize=13, fontweight="bold")
    ax.scatter(df_nodes["x"], newton_err, color=C_NDD, marker="o", s=60, zorder=3,
               label=f"Newton max = {newton_err.max():.3e}")
    ax.scatter(df_nodes["x"], chebyshev_err, color=C_CHB, marker="s", s=60, zorder=3,
               label=f"Chebyshev max = {chebyshev_err.max():.3e}")
    ax.set_yscale("log")
    ax.set_xlabel("x")
    ax.set_ylabel("|error|")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    newton_rms = np.sqrt(np.mean(newton_err**2))
    chebyshev_rms = np.sqrt(np.mean(chebyshev_err**2))
    summary = (f"Newton: max = {newton_err.max():.4e}, RMS = {newton_rms:.4e}\n"
               f"Chebyshev: max = {chebyshev_err.max():.4e}, RMS = {chebyshev_rms:.4e}")
    fig.text(0.5, -0.04, summary, ha="center", fontsize=9,
             bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray"))
    plt.tight_layout()
    path = OUTPUT_DIR / "task4_1d_nodes_scatter.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[✓] Saved {path}")
    plt.show()


def main():
    dataset_filename = active_dataset_filename()
    df = load_files()
    newton_err, chebyshev_err = plot_absolute_error(df, dataset_filename)
    newton_rms, chebyshev_rms = plot_rms_error(df, newton_err, chebyshev_err)
    plot_max_rms_comparison(newton_err, chebyshev_err, newton_rms, chebyshev_rms)
    plot_gradient_error(dataset_filename, float(df["x"].min()), float(df["x"].max()))
    plot_nodes_error(dataset_filename)


def load_files():
    comparison_file = OUTPUT_DIR / "comparison.csv"
    if not comparison_file.is_file():
        raise FileNotFoundError(
            f"{comparison_file} was not found. Run the C++ main program first."
        )
    df = pd.read_csv(comparison_file)
    required = {"x", "Newton", "Chebyshev"}
    if not required.issubset(df.columns) or df.empty:
        raise ValueError(
            f"Invalid comparison CSV: {comparison_file}. Expected columns x, Newton, Chebyshev."
        )
    print(f"[✓] Loaded '{comparison_file}' ({len(df)} points)")
    return df


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(f"[ERROR] {error}") from error
