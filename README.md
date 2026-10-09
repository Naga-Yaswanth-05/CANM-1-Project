# Numerical Methods — Polynomial Interpolation, Approximation & Root Finding

**Course:** Computer Aided Numerical Methods – I

Full theory, derivations and plots are in [`doc/canm_project_report.pdf`](doc/canm_project_report.pdf).

## What this project does

1. **Interpolation:** Newton divided differences and Newton forward differences.
2. **Chebyshev approximation:** normal equations, gradient descent, and QR (Householder).
3. **Root finding:** Newton–Raphson and Secant methods applied to the derivative of the fitted polynomial.
4. **Convergence analysis:** successive step errors and estimated convergence orders.

## Datasets

All input datasets are stored in `data/`:

| Dataset argument | Function | File |
|---|---|---|
| `data.csv` | `sin(5x)` | [`data/data.csv`](data/data.csv) |
| `group1.csv` | `1/(1+25x²)` | [`data/group1.csv`](data/group1.csv) |
| `tanh_data.csv` | `tanh(10x)` | [`data/tanh_data.csv`](data/tanh_data.csv) |

Each file contains a header row followed by `x,y` rows. The C++ pipeline validates the dataset before processing it.

## Path handling

C++ paths are resolved to absolute paths from the project location discovered using the executable location. Python scripts resolve the project root from their own file locations. You **do not** need to edit dataset paths, use `cd` into specific folders, or hardcode your Windows username. Input files are read from `data/`; generated CSVs and plots are written to `outputs/`.

## Build and run the main pipeline

Run these commands from the project root.

### Windows PowerShell with g++

```powershell
g++ -O2 -std=c++17 src/main.cpp -o main.exe
.\main.exe --dataset data.csv --degree 12 --method 1
```

### Git Bash, Linux, or macOS with g++

```bash
g++ -O2 -std=c++17 src/main.cpp -o main
./main --dataset data.csv --degree 12 --method 1
```

Arguments:

- `--dataset`: `data.csv`, `group1.csv`, or `tanh_data.csv` (default: `data.csv`).
- `--degree`: Chebyshev polynomial degree (default: `12`; must be smaller than the number of input points).
- `--method`: `1` = QR, `2` = Normal Equations, `3` = Gradient Descent. If omitted, the program prompts for a method.
- `--no-plots`: generate the data files without automatically launching the interpolation plots.

For example, to use the `tanh(10x)` dataset with degree 12 and QR:

```bash
./main --dataset tanh_data.csv --degree 12 --method 1 --no-plots
```

On Windows PowerShell, use ` .\main.exe ` instead of ` ./main ` if the executable is named `main.exe`.

The program writes `comparison.csv`, `comparison_nodes.csv`, polynomial coefficients, root-finding results, a dataset marker, and a method-specific comparison CSV under `outputs/`. By default, it launches `scripts/interpolation_gradient.py` after generating the results. Add `--no-plots` to run without launching the plotting script.

> The remaining command examples use Git Bash/Linux/macOS syntax. In Windows PowerShell, use `.\main.exe` instead of `./main`, give compiled C++ executables an `.exe` suffix, and launch them with PowerShell syntax (for example, `.\convergence_analysis\err_anal.exe`).

## Compare the three Chebyshev methods

Run the main program once per method, using the **same dataset and degree** each time. Example for `data.csv`:

```bash
./main --dataset data.csv --degree 12 --method 1 --no-plots
./main --dataset data.csv --degree 12 --method 2 --no-plots
./main --dataset data.csv --degree 12 --method 3 --no-plots
python scripts/compare_chebyshev.py
```

The C++ program automatically saves the method comparisons as `outputs/cheby_qr_sin.csv`, `outputs/cheby_normal_sin.csv`, and `outputs/cheby_gradient_sin.csv`. There is no manual renaming. For the other datasets, the suffix is `group1` or `tanh`.

## Degree versus RMS error

Every QR run automatically saves a degree-specific comparison file. For `data.csv`, create the six files required by the degree plot using QR:

```bash
./main --dataset data.csv --degree 6 --method 1 --no-plots
./main --dataset data.csv --degree 9 --method 1 --no-plots
./main --dataset data.csv --degree 12 --method 1 --no-plots
./main --dataset data.csv --degree 15 --method 1 --no-plots
./main --dataset data.csv --degree 18 --method 1 --no-plots
./main --dataset data.csv --degree 24 --method 1 --no-plots
python scripts/degree_vs_rms_plot.py
```

The script reads the files from `outputs/`, computes errors against the correct analytical function for the selected dataset, and saves the plot there. Use the same `--dataset` argument for every degree in one sweep.

## Convergence analysis

First run the main pipeline for the dataset and Chebyshev method whose coefficients you want to analyze. Then, from the project root:

```bash
g++ -O2 -std=c++17 convergence_analysis/err_anal_secant_newton.cpp -o convergence_analysis/err_anal
./convergence_analysis/err_anal
python scripts/conv_anal_secant_newton.py
```

This produces four error/order CSV files in `outputs/` and the plot `outputs/convergence_plot.png`.

An additional root-iteration table is available from `order_of_convergence_for_roots.cpp`:

```bash
g++ -O2 -std=c++17 convergence_analysis/order_of_convergence_for_roots.cpp -o convergence_analysis/order_of_convergence
./convergence_analysis/order_of_convergence
```

It writes `secant_output.csv` and `newton_output.csv` under `outputs/`. The displayed true error uses the known stationary point for `sin(5x)` and `1/(1+25x²)`; it is reported as not applicable for `tanh(10x)`, whose derivative has no finite zero.

## Project structure

```text
CANM-1-Project/
├── convergence_analysis/
│   ├── err_anal_secant_newton.cpp
│   └── order_of_convergence_for_roots.cpp
├── data/
│   ├── data.csv
│   ├── group1.csv
│   └── tanh_data.csv
├── doc/
│   ├── canm_project_report.pdf
│   ├── Contributions.xlsx
│   └── Instructions.docx
├── outputs/                  # Generated CSVs and plots (gitignored)
│   └── .gitkeep
├── scripts/
│   ├── project_utils.py       # Shared absolute paths and analytical functions
│   ├── interpolation_gradient.py
│   ├── compare_chebyshev.py
│   ├── conv_anal_secant_newton.py
│   └── degree_vs_rms_plot.py
├── src/
│   ├── main.cpp              # Main pipeline
│   ├── project_paths.h       # C++ absolute-path resolution
│   ├── utils.h
│   ├── newton.h
│   ├── newton_forward_difference.cpp
│   ├── NR_Secant.h
│   ├── Cheby_normal.h
│   ├── Cheby_gradient.h
│   ├── qr.h
│   └── evaluate.h
├── .gitignore
├── main.cpp                  # Compatibility entry point forwarding to src/main.cpp
└── README.md
```

## Dependencies

- **C++:** C++17 compiler such as `g++`; the numerical code uses the standard library.
- **Python:** `numpy`, `pandas`, and `matplotlib`.

Install Python dependencies with:

```bash
python -m pip install numpy pandas matplotlib
```

## Notes

- Polynomial coefficients are written in ascending power order: `c[0] + c[1]x + c[2]x² + ...`.
- Chebyshev approximation coefficients saved in `chebyshev_coeffs.csv` have been converted to the original `x` domain.
- Generated files are kept in `outputs/` and excluded from Git; the source datasets in `data/` remain tracked.
- `main.cpp` in the project root is only a compatibility entry point. The maintained implementation is `src/main.cpp`.
