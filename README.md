# Numerical Methods — Polynomial Interpolation, Approximation & Root Finding

Computer Aided Numerical Methods – I | **IIT Hyderabad**

> Full theory, derivations and plots → `docs/canm_project_report.pdf`

---

## What This Project Does

1. **Interpolation** — Newton Divided Difference + Newton Forward Difference
2. **Chebyshev Approximation** — via Normal Equations, Gradient Descent, and QR (Householder)
3. **Root Finding** — Newton-Raphson & Secant on the derivative P′(x)
4. **Convergence Analysis** — Successive error & convergence order for both root-finding methods

Dataset: 25 equally-spaced points of `sin(5x)` on `[-1, 1]`

---

## How to Run

### Step 1 — Main pipeline
```bash
cd src
g++ -O2 -std=c++17 main.cpp -o main
./main
```
Choose Chebyshev method when prompted:
- `1` → QR *(recommended)*
- `2` → Normal Equations
- `3` → Gradient Descent

Outputs CSVs and auto-launches `interpolation_gradient.py`.

---

### Step 2 — Chebyshev comparison plot
Run `./main` three times (once per method). After each run rename the output:

| Run | Rename to |
|---|---|
| QR | `cheby_qr.csv` |
| Normal Equations | `cheby_normal.csv` |
| Gradient Descent | `cheby_gradient.csv` |

```bash
cd scripts
python compare_chebyshev.py
```

---

### Step 3 — Degree vs RMS plot
Run `./main` with degree set to `6, 9, 12, 15, 18, 24` (edit `main.cpp`).
After each run rename `comparison.csv` → `degree{N}_sin.csv`.

```bash
cd scripts
python degree_vs_rms_plot.py
```

---

### Step 4 — Convergence analysis
```bash
cd convergence_analysis
g++ -O2 -std=c++17 err_anal_secant_newton.cpp -o err_anal
./err_anal
python ../scripts/conv_anal_secant_newton.py
```

---

## Data Flow

```
data.csv
   └─► main.cpp
           ├── comparison.csv ──────────────► interpolation_gradient.py
           ├── comparision_25.csv ──────────► interpolation_gradient.py
           ├── *_polynomial_coefficients.csv
           ├── Newton_Raphson_Roots.csv ─────┐
           └── Secant_Roots.csv ─────────────┤
                                             ▼
                              err_anal_secant_newton.cpp
                                             │
                              *_err_ord.csv ─► conv_anal_secant_newton.py

cheby_qr/normal/gradient.csv ──────────────► compare_chebyshev.py
degree*_sin.csv ────────────────────────────► degree_vs_rms_plot.py
```

---

## Project Structure

```
numerical-methods-project/
├── src/
│   ├── main.cpp                          # Entry point
│   ├── utils.h                           # Scaling & basis conversion utilities
│   ├── newton.h                          # Newton Divided Difference
│   ├── newton_forward_difference.cpp     # Forward Difference (equally spaced)
│   ├── NR_Secant.h                       # Newton-Raphson & Secant root finding
│   ├── Cheby_normal.h                    # Chebyshev via Normal Equations
│   ├── Cheby_gradient.h                  # Chebyshev via Gradient Descent
│   ├── qr.h                              # Chebyshev via QR (Householder)
│   └── evaluate.h                        # Clenshaw, Horner, Newton evaluation
├── scripts/
│   ├── interpolation_gradient.py         # Error & gradient accuracy plots
│   ├── compare_chebyshev.py              # RMS comparison of 3 Chebyshev solvers
│   ├── conv_anal_secant_newton.py        # Convergence order plots
│   └── degree_vs_rms_plot.py            # Degree vs RMS error plot
├── convergence_analysis/
│   ├── err_anal_secant_newton.cpp        # Successive error & convergence order
│   └── order_of_convergence_for_roots.cpp
├── data/
│   ├── data.csv                          # Group 2: sin(5x), 25 points
│   ├── tanh_data.csv                     # Custom: tanh(10x)
│   └── group1.csv                        # Group 1: 1/(1+25x²)
├── outputs/                              # Auto-generated (gitignored)
│   └── .gitkeep
├── docs/
│   ├── canm_project_report.pdf
│   ├── Contributions.xlsx
│   └── Instructions.docx
├── .gitignore
└── README.md
```

---

## Key Results

| Method | Max Error | RMS Error |
|---|---|---|
| Newton Divided Difference (degree 24) | 2.63e-12 | 3.35e-13 |
| Chebyshev QR (degree 12) | 5.16e-06 | 1.11e-06 |

All 3 Chebyshev solvers achieve identical RMS of **1.611747e-05** on training points.
QR is preferred — it avoids squaring the condition number, making it most stable.

| Root-finding Method | Convergence Order |
|---|---|
| Newton-Raphson | ~2 (Quadratic) |
| Secant | ~1.618 (Superlinear) |

---

## Dependencies

```bash
# C++ — standard library only, no external deps
g++ -O2 -std=c++17 main.cpp -o main

# Python
pip install numpy pandas matplotlib openpyxl
```

---

## Notes

- Input CSV format: header row, then `x,y` columns
- Input file is hardcoded as `data.csv` in `main.cpp` — change to switch datasets
- All polynomial coefficients are stored in **ascending power order**: `c[0] + c[1]x + c[2]x² ...`
- Chebyshev coefficients in output CSVs are already converted to the **original x domain** — evaluate directly with Horner's method
- Python scripts use `os.chdir` so run them from inside their own folder
