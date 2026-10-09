#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../src/project_paths.h"

using namespace std;
using Real = long double;

vector<Real> readCoefficients(const filesystem::path& path) {
    vector<Real> coefficients;
    ifstream file(path);
    if (!file.is_open()) {
        cerr << "[ERROR] Could not open coefficient file: " << path << '\n';
        return coefficients;
    }

    string line;
    while (getline(file, line)) {
        if (line.empty()) continue;
        try {
            size_t consumed = 0;
            Real value = stold(line, &consumed);
            if (consumed == line.size()) coefficients.push_back(value);
        } catch (const exception&) {
            // Ignore blank/non-numeric lines rather than breaking the load.
        }
    }
    return coefficients;
}

Real evaluatePolynomial(const vector<Real>& coefficients, Real x) {
    Real value = 0;
    for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
        value = value * x + *it;
    }
    return value;
}

vector<Real> derivativeCoefficients(const vector<Real>& coefficients) {
    vector<Real> derivative;
    if (coefficients.size() < 2) return derivative;
    derivative.reserve(coefficients.size() - 1);
    for (size_t i = 1; i < coefficients.size(); ++i) {
        derivative.push_back(static_cast<Real>(i) * coefficients[i]);
    }
    return derivative;
}

vector<Real> secantHistory(const vector<Real>& function, Real x0, Real x1) {
    vector<Real> history{x0, x1};
    constexpr int MAX_ITERATIONS = 50;
    constexpr Real TOLERANCE = 1e-12L;
    constexpr Real DENOMINATOR_TOLERANCE = 1e-18L;

    for (int iteration = 0; iteration < MAX_ITERATIONS; ++iteration) {
        const Real f0 = evaluatePolynomial(function, x0);
        const Real f1 = evaluatePolynomial(function, x1);
        const Real denominator = f1 - f0;
        if (!isfinite(f0) || !isfinite(f1) || fabsl(denominator) < DENOMINATOR_TOLERANCE) {
            break;
        }

        const Real next = x1 - f1 * (x1 - x0) / denominator;
        if (!isfinite(next)) break;
        history.push_back(next);
        if (fabsl(next - x1) < TOLERANCE) break;
        x0 = x1;
        x1 = next;
    }
    return history;
}

vector<Real> newtonHistory(const vector<Real>& function,
                           const vector<Real>& function_derivative,
                           Real initial_guess) {
    vector<Real> history{initial_guess};
    Real x = initial_guess;
    constexpr int MAX_ITERATIONS = 50;
    constexpr Real TOLERANCE = 1e-12L;
    constexpr Real DERIVATIVE_TOLERANCE = 1e-18L;

    for (int iteration = 0; iteration < MAX_ITERATIONS; ++iteration) {
        const Real fx = evaluatePolynomial(function, x);
        const Real dfx = evaluatePolynomial(function_derivative, x);
        if (!isfinite(fx) || !isfinite(dfx) || fabsl(dfx) < DERIVATIVE_TOLERANCE) {
            break;
        }

        const Real next = x - fx / dfx;
        if (!isfinite(next)) break;
        history.push_back(next);
        if (fabsl(next - x) < TOLERANCE) break;
        x = next;
    }
    return history;
}

vector<Real> successiveErrors(const vector<Real>& iterates) {
    vector<Real> errors;
    for (size_t i = 0; i + 1 < iterates.size(); ++i) {
        errors.push_back(fabsl(iterates[i + 1] - iterates[i]));
    }
    return errors;
}

// Estimate p from three successive step errors e[n-2], e[n-1], e[n]:
// p ~= log(e[n]/e[n-1]) / log(e[n-1]/e[n-2]).
vector<Real> convergenceOrders(const vector<Real>& errors) {
    vector<Real> orders(errors.size(), 0.0L);
    for (size_t i = 2; i < errors.size(); ++i) {
        const Real e0 = errors[i - 2];
        const Real e1 = errors[i - 1];
        const Real e2 = errors[i];
        if (e0 <= 0 || e1 <= 0 || e2 <= 0) continue;

        const Real denominator = logl(e1 / e0);
        const Real numerator = logl(e2 / e1);
        if (!isfinite(denominator) || !isfinite(numerator) || fabsl(denominator) < 1e-12L) {
            continue;
        }
        const Real order = numerator / denominator;
        // Discard nonphysical estimates produced by round-off after convergence.
        if (isfinite(order) && order > 0 && order < 10) orders[i] = order;
    }
    return orders;
}

bool saveVector(const filesystem::path& path, const vector<Real>& values) {
    ofstream file(path);
    if (!file.is_open()) {
        cerr << "[ERROR] Could not write: " << path << '\n';
        return false;
    }
    file << scientific << setprecision(18);
    for (Real value : values) file << value << '\n';
    return true;
}

bool saveErrorOrderCSV(const filesystem::path& path,
                       const vector<Real>& errors,
                       const vector<Real>& orders) {
    ofstream file(path);
    if (!file.is_open()) {
        cerr << "[ERROR] Could not write: " << path << '\n';
        return false;
    }
    file << "iteration,error,order\n" << scientific << setprecision(15);
    for (size_t i = 0; i < errors.size(); ++i) {
        file << i + 1 << ',' << errors[i] << ',' << orders[i] << '\n';
    }
    return true;
}

bool analyseMethod(const string& label,
                   const vector<Real>& polynomial_coefficients,
                   const string& root_file_name,
                   const string& error_file_name,
                   bool use_newton) {
    const vector<Real> first_derivative = derivativeCoefficients(polynomial_coefficients);
    const vector<Real> second_derivative = derivativeCoefficients(first_derivative);
    if (first_derivative.empty()) {
        cerr << "[ERROR] Polynomial has insufficient coefficients for root finding: "
             << label << '\n';
        return false;
    }

    vector<Real> iterates;
    if (use_newton) {
        iterates = newtonHistory(first_derivative, second_derivative, 0.15L);
    } else {
        iterates = secantHistory(first_derivative, 0.1L, 0.2L);
    }

    if (iterates.size() < 2) {
        cerr << "[ERROR] Root iteration did not advance for " << label << '\n';
        return false;
    }

    const vector<Real> errors = successiveErrors(iterates);
    const vector<Real> orders = convergenceOrders(errors);
    const filesystem::path root_path = project_paths::output_file(root_file_name);
    const filesystem::path error_path = project_paths::output_file(error_file_name);

    if (!saveVector(root_path, iterates) || !saveErrorOrderCSV(error_path, errors, orders)) {
        return false;
    }

    cout << "\n" << label << "\n";
    cout << left << setw(8) << "Iter" << setw(24) << "Root estimate"
         << setw(22) << "Step error" << setw(14) << "Order" << '\n';
    cout << string(68, '-') << '\n';
    for (size_t i = 0; i < iterates.size(); ++i) {
        cout << left << setw(8) << i + 1 << setw(24) << fixed << setprecision(14)
             << iterates[i];
        if (i == 0) {
            cout << setw(22) << "---" << setw(14) << "---" << '\n';
        } else {
            const Real error = errors[i - 1];
            const Real order = (i - 1 < orders.size()) ? orders[i - 1] : 0.0L;
            cout << scientific << setprecision(6) << setw(22) << error;
            if (order > 0) cout << fixed << setprecision(6) << setw(14) << order;
            else cout << setw(14) << "---";
            cout << '\n';
        }
    }
    cout << "Saved: " << root_path << '\n' << "Saved: " << error_path << '\n';
    return true;
}

int main() {
    try {
        const vector<Real> ndd_coefficients = readCoefficients(
            project_paths::output_file("polynomial_coefficients.csv"));
        const vector<Real> chebyshev_coefficients = readCoefficients(
            project_paths::output_file("chebyshev_coeffs.csv"));

        if (ndd_coefficients.empty() || chebyshev_coefficients.empty()) {
            cerr << "[ERROR] Missing polynomial coefficients. Run the main program first.\n"
                 << "Expected files:\n  "
                 << project_paths::output_file("polynomial_coefficients.csv") << "\n  "
                 << project_paths::output_file("chebyshev_coeffs.csv") << '\n';
            return 1;
        }

        bool ok = true;
        ok &= analyseMethod("Secant — Newton divided-difference P'(x)",
                            ndd_coefficients, "secant_ndd_roots.csv",
                            "secant_ndd_err_ord.csv", false);
        ok &= analyseMethod("Newton-Raphson — Newton divided-difference P'(x)",
                            ndd_coefficients, "newton_rp_ndd_roots.csv",
                            "newton_rp_ndd_err_ord.csv", true);
        ok &= analyseMethod("Secant — Chebyshev P'(x)",
                            chebyshev_coefficients, "secant_chebyshev_roots.csv",
                            "secant_chebyshev_err_ord.csv", false);
        ok &= analyseMethod("Newton-Raphson — Chebyshev P'(x)",
                            chebyshev_coefficients, "newton_rp_chebyshev_roots.csv",
                            "newton_rp_chebyshev_err_ord.csv", true);

        return ok ? 0 : 1;
    } catch (const exception& e) {
        cerr << "[ERROR] " << e.what() << '\n';
        return 1;
    }
}
