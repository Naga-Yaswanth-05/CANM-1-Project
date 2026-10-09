"""Plot Chebyshev max and RMS errors versus degree for the active dataset."""
import matplotlib.pyplot as plt
import numpy as np

from project_utils import (
    OUTPUT_DIR,
    active_dataset_filename,
    dataset_key,
    dataset_label,
    exact_function,
)

DEGREES = [6, 9, 12, 15, 18, 24]


def main():
    dataset_filename = active_dataset_filename()
    key = dataset_key(dataset_filename)
    csv_files = {degree: OUTPUT_DIR / f"degree{degree}_{key}.csv" for degree in DEGREES}
    missing = [str(path) for path in csv_files.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing degree-sweep CSV file(s):\n  " + "\n  ".join(missing) +
            "\nRun the C++ program with --degree set to each of: " +
            ", ".join(map(str, DEGREES)) + f" and --dataset {dataset_filename}."
        )

    max_errors = []
    rms_errors = []
    reference_x = None
    for degree in DEGREES:
        data = np.genfromtxt(csv_files[degree], delimiter=",", skip_header=1)
        if data.ndim != 2 or data.shape[1] < 3 or data.shape[0] == 0:
            raise ValueError(f"Invalid comparison data: {csv_files[degree]}")
        if not np.isfinite(data[:, :3]).all():
            raise ValueError(f"Non-numeric or non-finite values found in {csv_files[degree]}")
        x = data[:, 0]
        if reference_x is None:
            reference_x = x
        elif x.shape != reference_x.shape or not np.allclose(x, reference_x, rtol=0, atol=1e-12):
            raise ValueError(f"The x grid in {csv_files[degree]} does not match the other degrees.")

        y_true = exact_function(dataset_filename, x)
        error = np.abs(data[:, 2] - y_true)
        max_errors.append(float(error.max()))
        rms_errors.append(float(np.sqrt(np.mean(error**2))))

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor("#F8F9FA")
    ax.set_facecolor("#F8F9FA")
    ax.semilogy(DEGREES, max_errors, "o-", color="red", linewidth=2,
                markersize=8, label="Max Error")
    ax.semilogy(DEGREES, rms_errors, "s-", color="blue", linewidth=2,
                markersize=8, label="RMS Error")
    for degree, error in zip(DEGREES, max_errors):
        ax.annotate(f"{error:.2e}", (degree, error), textcoords="offset points",
                    xytext=(0, 10), ha="center", fontsize=8, color="red")
    for degree, error in zip(DEGREES, rms_errors):
        ax.annotate(f"{error:.2e}", (degree, error), textcoords="offset points",
                    xytext=(0, -16), ha="center", fontsize=8, color="blue")
    ax.set_xlabel("Polynomial Degree", fontsize=11)
    ax.set_ylabel("Error (log scale)", fontsize=11)
    ax.set_title(f"Chebyshev — Max & RMS Error vs Degree ({dataset_label(dataset_filename)})",
                 fontsize=13, fontweight="bold")
    ax.set_xticks(DEGREES)
    ax.legend(fontsize=10)
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plot_path = OUTPUT_DIR / f"chebyshev_error_vs_degree_{key}.png"
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"Saved {plot_path}")
    plt.show()


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(f"[ERROR] {error}") from error
