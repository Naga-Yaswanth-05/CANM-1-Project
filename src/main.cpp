#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "project_paths.h"
#include "utils.h"
#include "newton.h"
#include "Cheby_normal.h"
#include "qr.h"
#include "evaluate.h"
#include "Cheby_gradient.h"
#include "NR_Secant.h"

using namespace std;
namespace fs = std::filesystem;

namespace {
constexpr int EVALUATION_POINTS = 600;

struct DatasetInfo {
    const char* filename;
    const char* key;
    const char* display_name;
};

const vector<DatasetInfo> DATASETS = {
    {"data.csv", "sin", "sin(5x)"},
    {"group1.csv", "group1", "1/(1+25x^2)"},
    {"tanh_data.csv", "tanh", "tanh(10x)"},
};

const DatasetInfo* find_dataset(const string& filename) {
    for (const auto& dataset : DATASETS) {
        if (filename == dataset.filename) return &dataset;
    }
    return nullptr;
}

void print_usage(const char* executable) {
    cout << "Usage: " << executable
         << " [--dataset data.csv|group1.csv|tanh_data.csv]"
         << " [--degree N] [--method 1|2|3] [--no-plots]\n\n"
         << "Methods: 1 = QR, 2 = Normal Equations, 3 = Gradient Descent\n"
         << "Defaults: --dataset data.csv --degree 12; method is prompted.\n";
}

bool load_dataset(const fs::path& path,
                  vector<long double>& x,
                  vector<long double>& y) {
    ifstream file(path);
    if (!file.is_open()) {
        cerr << "Error: Could not open dataset: " << path << '\n';
        return false;
    }

    string line;
    size_t line_number = 0;
    bool first_content_row = true;
    while (getline(file, line)) {
        ++line_number;
        if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
        if (first_content_row) {
            first_content_row = false;
            const string trimmed = line.substr(line.find_first_not_of(" \t\r\n"));
            if (trimmed.find("x") != string::npos && trimmed.find("y") != string::npos &&
                trimmed.find(',') != string::npos) {
                continue; // Header row: x,y
            }
        }

        stringstream row(line);
        long double x_value = 0.0L, y_value = 0.0L;
        char comma = '\0';
        if (!(row >> x_value >> comma >> y_value) || comma != ',' ||
            !isfinite(x_value) || !isfinite(y_value)) {
            cerr << "Error: Invalid CSV row at " << path << ':' << line_number
                 << ". Expected two finite numeric values separated by a comma.\n";
            return false;
        }
        row >> ws;
        if (!row.eof()) {
            cerr << "Error: Unexpected extra data at " << path << ':' << line_number << '\n';
            return false;
        }

        x.push_back(x_value);
        y.push_back(y_value);
    }

    if (x.size() < 2 || x.size() != y.size()) {
        cerr << "Error: Dataset must contain at least two valid x,y rows: " << path << '\n';
        return false;
    }
    vector<long double> sorted_x = x;
    sort(sorted_x.begin(), sorted_x.end());
    if (adjacent_find(sorted_x.begin(), sorted_x.end()) != sorted_x.end()) {
        cerr << "Error: Dataset contains duplicate x values, so interpolation is undefined: "
             << path << '\n';
        return false;
    }
    return true;
}

bool write_active_dataset(const string& filename) {
    const fs::path marker_path = project_paths::output_file("active_dataset.txt");
    ofstream marker(marker_path);
    if (!marker.is_open()) {
        cerr << "Error: Could not write active-dataset marker: " << marker_path << '\n';
        return false;
    }
    marker << filename << '\n';
    return static_cast<bool>(marker);
}

bool copy_file_overwrite(const fs::path& source, const fs::path& destination) {
    error_code ec;
    fs::copy_file(source, destination, fs::copy_options::overwrite_existing, ec);
    if (ec) {
        cerr << "Error: Could not save " << destination << " (" << ec.message() << ")\n";
        return false;
    }
    return true;
}
} // namespace

int main(int argc, char* argv[]) {
    try {
        string dataset_filename = "data.csv";
        int degree = 12;
        int method_choice = 0;
        bool launch_interpolation_plot = true;

        for (int i = 1; i < argc; ++i) {
            const string argument = argv[i];
            if (argument == "--help" || argument == "-h") {
                print_usage(argv[0]);
                return 0;
            } else if (argument == "--no-plots") {
                launch_interpolation_plot = false;
            } else if (argument == "--dataset" || argument == "--degree" || argument == "--method") {
                if (i + 1 >= argc) {
                    cerr << "Error: Missing value after " << argument << ".\n";
                    print_usage(argv[0]);
                    return 1;
                }
                const string value = argv[++i];
                try {
                    if (argument == "--dataset") {
                        dataset_filename = value;
                    } else if (argument == "--degree") {
                        size_t consumed = 0;
                        degree = stoi(value, &consumed);
                        if (consumed != value.size()) throw invalid_argument("invalid degree");
                    } else {
                        size_t consumed = 0;
                        method_choice = stoi(value, &consumed);
                        if (consumed != value.size()) throw invalid_argument("invalid method");
                    }
                } catch (const exception&) {
                    cerr << "Error: Invalid value for " << argument << ": " << value << '\n';
                    return 1;
                }
            } else {
                cerr << "Error: Unrecognized argument: " << argument << '\n';
                print_usage(argv[0]);
                return 1;
            }
        }

        const DatasetInfo* dataset_info = find_dataset(dataset_filename);
        if (dataset_info == nullptr) {
            cerr << "Error: Unsupported dataset '" << dataset_filename << "'. Choose one of:\n";
            for (const auto& dataset : DATASETS) cerr << "  " << dataset.filename << '\n';
            return 1;
        }

        const fs::path data_path = project_paths::data_file(dataset_filename);
        vector<long double> x, y;
        if (!load_dataset(data_path, x, y)) return 1;
        if (degree < 0 || degree >= static_cast<int>(x.size())) {
            cerr << "Error: Degree must be between 0 and " << x.size() - 1
                 << " for this dataset.\n";
            return 1;
        }
        if (method_choice != 0 && (method_choice < 1 || method_choice > 3)) {
            cerr << "Error: Method must be 1 (QR), 2 (Normal Equations), or 3 (Gradient Descent).\n";
            return 1;
        }
        cout << "Dataset: " << data_path << " (" << x.size() << " points; function "
             << dataset_info->display_name << ")\n";
        cout << "\n=============================\nNewton Polynomial\n=============================\n";
        vector<long double> newton_coeff = newton_divided_difference(x, y);
        if (newton_coeff.empty()) {
            cerr << "Error: Newton divided-difference coefficients could not be generated.\n";
            return 1;
        }

        cout << "\n=============================\nChebyshev Polynomial\n=============================\n";
        cout << "Polynomial degree: " << degree << '\n';

        if (method_choice == 0) {
            cout << "Choose Chebyshev Method:\n"
                 << "1. QR Method\n"
                 << "2. Normal Equations\n"
                 << "3. Gradient Descent\n"
                 << "Enter your choice: ";
            if (!(cin >> method_choice)) {
                cerr << "Error: Please enter 1, 2, or 3.\n";
                return 1;
            }
        }
        if (method_choice < 1 || method_choice > 3) {
            cerr << "Error: Please enter 1, 2, or 3.\n";
            return 1;
        }

        vector<long double> cheb_coeff;
        string method_suffix;
        if (method_choice == 1) {
            cheb_coeff = chebyshev_householder(x, y, degree);
            method_suffix = "qr";
        } else if (method_choice == 2) {
            cheb_coeff = Chebyshev_Using_Normal_Equation(x, y, degree);
            method_suffix = "normal";
        } else {
            cheb_coeff = Chebyshev_Using_Gradient_descent(x, y, degree);
            method_suffix = "gradient";
        }
        if (cheb_coeff.empty()) {
            cerr << "Error: Chebyshev coefficients could not be generated.\n";
            return 1;
        }

        const long double a = *min_element(x.begin(), x.end());
        const long double b = *max_element(x.begin(), x.end());
        const fs::path comparison_path = project_paths::output_file("comparison.csv");
        ofstream out(comparison_path);
        if (!out.is_open()) {
            cerr << "Error: Could not write comparison file: " << comparison_path << '\n';
            return 1;
        }

        out << fixed << setprecision(15) << "x,Newton,Chebyshev\n";
        for (int i = 0; i < EVALUATION_POINTS; ++i) {
            const long double xi = a + (b - a) * i / (EVALUATION_POINTS - 1);
            const long double y_newton = evaluate_newton(xi, x, newton_coeff);
            const long double y_cheb = evaluate_chebyshev(xi, cheb_coeff, a, b);
            out << xi << ',' << y_newton << ',' << y_cheb << '\n';
        }
        out.close();
        if (!out) {
            cerr << "Error: Failed while writing comparison file: " << comparison_path << '\n';
            return 1;
        }

        const fs::path method_comparison_path = project_paths::output_file(
            "cheby_" + method_suffix + "_" + dataset_info->key + ".csv");
        if (!copy_file_overwrite(comparison_path, method_comparison_path)) return 1;

        const fs::path nodes_path = project_paths::output_file("comparison_nodes.csv");
        ofstream nodes_out(nodes_path);
        if (!nodes_out.is_open()) {
            cerr << "Error: Could not write node comparison file: " << nodes_path << '\n';
            return 1;
        }
        nodes_out << fixed << setprecision(15);
        for (size_t i = 0; i < x.size(); ++i) {
            nodes_out << x[i] << ',' << evaluate_newton(x[i], x, newton_coeff) << ','
                      << evaluate_chebyshev(x[i], cheb_coeff, a, b) << '\n';
        }
        nodes_out.close();
        if (!nodes_out) {
            cerr << "Error: Failed while writing node comparison file: " << nodes_path << '\n';
            return 1;
        }

        cout << "\nComparison file generated: " << comparison_path << '\n'
             << "Method comparison file generated: " << method_comparison_path << '\n'
             << "Node comparison file generated: " << nodes_path << '\n';

        // Degree-versus-error plots use QR consistently so solver changes cannot
        // silently replace a degree-sweep result with another method's output.
        if (method_choice == 1) {
            const fs::path degree_path = project_paths::output_file(
                "degree" + to_string(degree) + "_" + dataset_info->key + ".csv");
            if (!copy_file_overwrite(comparison_path, degree_path)) return 1;
            cout << "Degree-sweep file generated: " << degree_path << '\n';
        }

        if (!write_active_dataset(dataset_filename)) return 1;

        if (launch_interpolation_plot) {
            cout << "\nOpening interpolation plots...\n";
            const fs::path script = project_paths::script_file("interpolation_gradient.py");
#ifdef _WIN32
            const string command = "python \"" + script.string() + "\"";
#else
            const string command = "python3 \"" + script.string() + "\"";
#endif
            const int script_status = system(command.c_str());
            if (script_status != 0) {
                cerr << "Warning: The plotting script did not complete successfully.\n"
                     << "Run it manually with Python: " << script << '\n';
            }
        }
        return 0;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
