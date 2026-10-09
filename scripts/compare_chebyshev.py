"""Compare the three Chebyshev solvers on the currently selected dataset."""
import matplotlib.pyplot as plt
import numpy as np

from project_utils import (
    OUTPUT_DIR,
    active_dataset_filename,
    dataset_key,
    dataset_label,
    exact_function,
)


def load_approximation(label, path):
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {label} comparison file: {path}. Run the C++ program once with "
            "--method 1, --method 2, and --method 3 using the same dataset."
        )
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    if data.ndim != 2 or data.shape[1] < 3 or data.shape[0] == 0:
        raise ValueError(f"Invalid comparison CSV format: {path}")
    if not np.isfinite(data[:, :3]).all():
        raise ValueError(f"Non-numeric or non-finite values found in {path}")
    return data[:, 0], data[:, 2]


def main():
    dataset_filename = active_dataset_filename()
    key = dataset_key(dataset_filename)
    method_files = {
        "QR": OUTPUT_DIR / f"cheby_qr_{key}.csv",
        "Normal Equations": OUTPUT_DIR / f"cheby_normal_{key}.csv",
        "Gradient Descent": OUTPUT_DIR / f"cheby_gradient_{key}.csv",
    }

    errors = {}
    reference_x = None
    for label, path in method_files.items():
        x, y_approx = load_approximation(label, path)
        if reference_x is None:
            reference_x = x
        elif (x.shape != reference_x.shape or
              not np.allclose(x, reference_x, rtol=0, atol=1e-12)):
            raise ValueError("The three method files use different x grids.")
        y_true = exact_function(dataset_filename, x)
        errors[label] = float(np.sqrt(np.mean((y_approx - y_true) ** 2)))

    methods = list(errors)
    rms_values = [errors[name] for name in methods]
    fig, ax = plt.subplots(figsize=(9, 6))
    bars = ax.bar(methods, rms_values, width=0.5, edgecolor="black", linewidth=0.8)
    for bar, value in zip(bars, rms_values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{value:.6e}", ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("RMS Error", fontsize=12)
    ax.set_title(f"Chebyshev Approximation RMS Error — {dataset_label(dataset_filename)}",
                 fontsize=13, fontweight="bold")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y:.2e}"))
    ax.grid(True, axis="y", linestyle="--", alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plot_path = OUTPUT_DIR / f"chebyshev_method_comparison_{key}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"Saved {plot_path}")
    plt.show()


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(f"[ERROR] {error}") from error
