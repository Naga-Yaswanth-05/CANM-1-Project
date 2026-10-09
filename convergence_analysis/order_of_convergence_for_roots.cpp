#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "../src/project_paths.h"

using namespace std;
namespace fs = std::filesystem;

optional<string> readActiveDataset() {
    const fs::path marker_path = project_paths::output_file("active_dataset.txt");
    ifstream marker(marker_path);
    if (!marker.is_open()) {
        cerr << "Error: Could not open " << marker_path
             << ". Run the main program first.\n";
        return nullopt;
    }
    string filename;
    getline(marker, filename);
    if (filename.empty()) {
        cerr << "Error: Active-dataset marker is empty: " << marker_path << '\n';
        return nullopt;
    }
    return filename;
}

optional<double> analyticalStationaryPoint(const string& dataset) {
    if (dataset == "data.csv") return acos(-1.0) / 10.0; // sin(5x): f'(x)=0 at pi/10
    if (dataset == "group1.csv") return 0.0;              // 1/(1+25x^2): f'(x)=0 at x=0
    if (dataset == "tanh_data.csv") return nullopt;       // tanh(10x) has no finite stationary point
    return nullopt;
}

vector<double> readCoefficientsCSV(const fs::path& filename) {
    vector<double> coefficients;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open " << filename << ". Run the main program first.\n";
        return coefficients;
    }

    string line;
    while (getline(file, line)) {
        if (line.empty()) continue;
        try {
            size_t consumed = 0;
            const double value = stod(line, &consumed);
            if (consumed == line.size()) coefficients.push_back(value);
        } catch (const exception&) {
            // Coefficient files are one number per line; ignore nonnumeric lines.
        }
    }
    return coefficients;
}

double evaluatePoly(double x, const vector<double>& coefficients) {
    double result = 0.0;
    for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
        result = result * x + *it;
    }
    return result;
}

vector<double> derivativeCoefficients(const vector<double>& coefficients) {
    vector<double> derivative;
    if (coefficients.size() < 2) return derivative;
    derivative.reserve(coefficients.size() - 1);
    for (size_t i = 1; i < coefficients.size(); ++i) {
        derivative.push_back(static_cast<double>(i) * coefficients[i]);
    }
    return derivative;
}

bool saveToCSV(const string& filename, const vector<double>& history) {
    const fs::path path = project_paths::output_file(filename);
    ofstream file(path);
    if (!file.is_open()) {
        cerr << "Error: Could not write " << path << '\n';
        return false;
    }
    file << setprecision(17);
    for (double value : history) file << value << '\n';
    return static_cast<bool>(file);
}

void printHistory(const string& method,
                  const vector<double>& history,
                  const optional<double>& target) {
    cout << "\n--- " << method << " Method (solving P'(x)=0) ---\n";
    cout << left << setw(8) << "Iter" << setw(20) << "x value"
         << setw(20) << "Step Error" << setw(20) << "True Error" << '\n';
    cout << string(68, '-') << '\n';
    for (size_t i = 0; i < history.size(); ++i) {
        cout << left << setw(8) << i << setw(20) << fixed << setprecision(12) << history[i];
        if (i == 0) cout << setw(20) << "---";
        else cout << setw(20) << abs(history[i] - history[i - 1]);
        if (target.has_value()) cout << setw(20) << abs(*target - history[i]) << '\n';
        else cout << setw(20) << "N/A" << '\n';
    }
}

vector<double> secantHistory(double x0, double x1, const vector<double>& derivative) {
    vector<double> history{x0, x1};
    for (int iteration = 0; iteration < 30; ++iteration) {
        const double f0 = evaluatePoly(x0, derivative);
        const double f1 = evaluatePoly(x1, derivative);
        const double denominator = f1 - f0;
        if (!isfinite(f0) || !isfinite(f1) || abs(denominator) < 1e-18) break;

        const double next = x1 - f1 * (x1 - x0) / denominator;
        if (!isfinite(next)) break;
        history.push_back(next);
        if (abs(next - x1) < 1e-12) break;
        x0 = x1;
        x1 = next;
    }
    return history;
}

vector<double> newtonHistory(double initial_guess,
                             const vector<double>& derivative,
                             const vector<double>& second_derivative) {
    vector<double> history{initial_guess};
    double x = initial_guess;
    for (int iteration = 0; iteration < 30; ++iteration) {
        const double fx = evaluatePoly(x, derivative);
        const double dfx = evaluatePoly(x, second_derivative);
        if (!isfinite(fx) || !isfinite(dfx) || abs(dfx) < 1e-18) break;

        const double next = x - fx / dfx;
        if (!isfinite(next)) break;
        history.push_back(next);
        if (abs(next - x) < 1e-12) break;
        x = next;
    }
    return history;
}

int main() {
    try {
        const optional<string> active_dataset = readActiveDataset();
        if (!active_dataset) return 1;
        const optional<double> target = analyticalStationaryPoint(*active_dataset);
        const vector<double> polynomial = readCoefficientsCSV(
            project_paths::output_file("polynomial_coefficients.csv"));
        if (polynomial.size() < 2) {
            cerr << "Error: Polynomial coefficient file is missing or has too few values.\n";
            return 1;
        }

        const vector<double> derivative = derivativeCoefficients(polynomial);
        const vector<double> second_derivative = derivativeCoefficients(derivative);
        const vector<double> secant = secantHistory(0.1, 0.2, derivative);
        const vector<double> newton = newtonHistory(0.15, derivative, second_derivative);

        printHistory("Secant", secant, target);
        printHistory("Newton-Raphson", newton, target);

        if (!saveToCSV("secant_output.csv", secant) ||
            !saveToCSV("newton_output.csv", newton)) return 1;
        cout << "\nResults saved in: " << project_paths::output_dir() << '\n';
        if (!target.has_value()) {
            cout << "True error is N/A for tanh(10x), whose derivative has no finite zero.\n";
        }
        return 0;
    } catch (const exception& error) {
        cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
