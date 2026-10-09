"""Plot the convergence-analysis CSVs produced by err_anal_secant_newton.cpp."""
import matplotlib.pyplot as plt
import pandas as pd

from project_utils import OUTPUT_DIR

FILES = {
    "Secant — Newton DD": OUTPUT_DIR / "secant_ndd_err_ord.csv",
    "Newton-Raphson — Newton DD": OUTPUT_DIR / "newton_rp_ndd_err_ord.csv",
    "Secant — Chebyshev": OUTPUT_DIR / "secant_chebyshev_err_ord.csv",
    "Newton-Raphson — Chebyshev": OUTPUT_DIR / "newton_rp_chebyshev_err_ord.csv",
}


def load_frame(label, path):
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {label} data: {path}. Compile and run "
            "convergence_analysis/err_anal_secant_newton.cpp after running the main program."
        )
    frame = pd.read_csv(path)
    required = {"iteration", "error", "order"}
    if not required.issubset(frame.columns) or frame.empty:
        raise ValueError(f"Invalid convergence CSV: {path}")
    return frame


def main():
    frames = {label: load_frame(label, path) for label, path in FILES.items()}
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Convergence Analysis: Secant vs Newton-Raphson",
                 fontsize=14, fontweight="bold")

    groups = [
        (axes[0][0], axes[0][1], frames["Secant — Newton DD"],
         frames["Newton-Raphson — Newton DD"], "Newton Divided Difference P′(x)"),
        (axes[1][0], axes[1][1], frames["Secant — Chebyshev"],
         frames["Newton-Raphson — Chebyshev"], "Chebyshev P′(x)"),
    ]

    for error_ax, order_ax, secant, newton, title in groups:
        error_ax.semilogy(secant["iteration"], secant["error"], "o-",
                          color="blue", label="Secant", linewidth=2, markersize=4)
        error_ax.semilogy(newton["iteration"], newton["error"], "s-",
                          color="magenta", label="Newton-Raphson", linewidth=2, markersize=4)
        error_ax.set_title(f"{title}\nSuccessive Step Error")
        error_ax.set_xlabel("Iteration")
        error_ax.set_ylabel("Error (log scale)")
        error_ax.legend()
        error_ax.grid(True, which="both", linestyle="--", alpha=0.5)

        secant_valid = secant[secant["order"] > 0]
        newton_valid = newton[newton["order"] > 0]
        order_ax.plot(secant_valid["iteration"], secant_valid["order"], "o-",
                      color="blue", label="Secant", linewidth=2, markersize=4)
        order_ax.plot(newton_valid["iteration"], newton_valid["order"], "s-",
                      color="magenta", label="Newton-Raphson", linewidth=2, markersize=4)
        order_ax.axhline(y=2.0, color="green", linestyle="--", alpha=0.6,
                         label="Order 2 (quadratic)")
        order_ax.axhline(y=1.618, color="orange", linestyle="--", alpha=0.6,
                         label="Order 1.618 (superlinear)")
        order_ax.set_title(f"{title}\nEstimated Convergence Order")
        order_ax.set_xlabel("Iteration")
        order_ax.set_ylabel("Order")
        order_ax.set_ylim(0, 3.5)
        order_ax.legend(fontsize=8, loc="lower right")
        order_ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plot_path = OUTPUT_DIR / "convergence_plot.png"
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"[✓] Saved {plot_path}")
    plt.show()


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(f"[ERROR] {error}") from error
