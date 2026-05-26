import pandas as pd
import matplotlib.pyplot as plt
import os


#  INPUT FILES
SEC_NDD_FILE = "secant_ndd_err_ord.csv"
NWT_NDD_FILE = "newton_rp_ndd_err_ord.csv"
SEC_CHB_FILE = "secant_chebyshev_err_ord.csv"
NWT_CHB_FILE = "newton_rp_chebyshev_err_ord.csv"

# Check all files exist
files = {
    "Secant — Newton DD":         SEC_NDD_FILE,
    "Newton-Raphson — Newton DD": NWT_NDD_FILE,
    "Secant — Chebyshev":         SEC_CHB_FILE,
    "Newton-Raphson — Chebyshev": NWT_CHB_FILE,
}
for label, f in files.items():
    if not os.path.exists(f):
        print(f"[ERROR] '{f}' not found.")
        exit(1)

# Load
df_sec_ndd = pd.read_csv(SEC_NDD_FILE)
df_nwt_ndd = pd.read_csv(NWT_NDD_FILE)
df_sec_chb = pd.read_csv(SEC_CHB_FILE)
df_nwt_chb = pd.read_csv(NWT_CHB_FILE)

#Plot
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("Convergence Analysis: Secant vs Newton-Raphson",
             fontsize=14, fontweight="bold")

datasets = [
    (axes[0][0], axes[0][1], df_sec_ndd, df_nwt_ndd, "Newton Divided Difference P′(x)"),
    (axes[1][0], axes[1][1], df_sec_chb, df_nwt_chb, "Chebyshev P′(x)"),
]

for ax_err, ax_ord, df_sec, df_nwt, title in datasets:
    
    #Successive Error (log scale)
    ax_err.semilogy(df_sec["iteration"], df_sec["error"], "o-",
                    color="blue", label="Secant", linewidth=2, markersize=4)
    ax_err.semilogy(df_nwt["iteration"], df_nwt["error"], "s-",
                    color="magenta", label="Newton-Raphson", linewidth=2, markersize=4)
    ax_err.set_title(f"{title}\nSuccessive Error $|x_{{n+1}} - x_n|$")
    ax_err.set_xlabel("Iteration")
    ax_err.set_ylabel("Error (Log Scale)")
    ax_err.legend()
    ax_err.grid(True, which="both", linestyle="--", alpha=0.5)

    # Convergence Order
    # Filter out rows where order is 0 (first row and precision-limit rows)
    df_sec_valid = df_sec[df_sec["order"] > 0]
    df_nwt_valid = df_nwt[df_nwt["order"] > 0]

    ax_ord.plot(df_sec_valid["iteration"], df_sec_valid["order"], "o-",
                color="blue", label="Secant", linewidth=2, markersize=4)
    ax_ord.plot(df_nwt_valid["iteration"], df_nwt_valid["order"], "s-",
                color="magenta", label="Newton-Raphson", linewidth=2, markersize=4)
    ax_ord.axhline(y=2.0,   color="green",  linestyle="--", alpha=0.6,
                   label="Order 2 (Quadratic)(Theoretical)")
    ax_ord.axhline(y=1.618, color="orange", linestyle="--", alpha=0.6,
                   label="Order 1.618 (Superlinear)(Theoretical)")
    ax_ord.set_title(f"{title}\nEstimated Convergence Order $p_n$")
    ax_ord.set_xlabel("Iteration")
    ax_ord.set_ylabel("Order")
    ax_ord.set_ylim(0, 3.5)
    ax_ord.legend(fontsize=8, loc="lower right")
    ax_ord.grid(True, linestyle="--", alpha=0.5)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("convergence_plot.png", dpi=150, bbox_inches="tight")
print("[✓] Saved convergence_plot.png")
plt.show()