#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>

using namespace std;

// This function reads raw numbers line-by-line from your input files.
vector<long double> readRootsFromCSV(const string& filename) {
    vector<long double> roots;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "[ERROR] Could not open " << filename << " for reading.\n";
        return roots;
    }

    string line;
    // We start reading from the very first line
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        stringstream ss(line);
        long double val;
        // Extracts the first numeric value found on the line
        if (ss >> val) {
            roots.push_back(val);
        }
    }

    file.close();
    return roots;
}

// Successive error: |x_{n+1} - x_n|
vector<long double> successiveErrors(const vector<long double>& roots) {
    vector<long double> errors;
    for (int i = 0; i + 1 < (int)roots.size(); i++)
        errors.push_back(fabsl(roots[i + 1] - roots[i]));
    return errors;
}

// Convergence order: p_n = log(e_{n+1}) / log(e_n)
vector<long double> convergenceOrders(const vector<long double>& errors) {
    vector<long double> orders;
    for (int i = 0; i + 1 < (int)errors.size(); i++) {
        long double e_curr = errors[i];
        long double e_next = errors[i + 1];
        // Ensure values are above machine noise (1e-19 for long double)
        if (e_curr > 1e-19L && e_next > 1e-19L)
            orders.push_back(logl(e_next) / logl(e_curr));
        else
            orders.push_back(0.0L); // Limit of precision reached
    }
    return orders;
}

// Print table to console
void printTable(const string& label,
                const vector<long double>& roots,
                const vector<long double>& errors,
                const vector<long double>& orders) {
    cout << "\n" << "\n";
    cout << "  " << label << "\n";
    cout << "\n";
    cout << setw(6)  << "Iter"
         << setw(20) << "Root Estimate"
         << setw(20) << "Succ. Error"
         << setw(16) << "Conv. Order" << "\n\n";

    if (roots.empty()) return;

    // Row 1 (Initial guess/first iteration)
    cout << setw(6)  << 1
         << setw(20) << fixed << setprecision(12) << roots[0]
         << setw(20) << "---"
         << setw(16) << "---" << "\n";

    for (int i = 0; i < (int)errors.size(); i++) {
        cout << setw(6)  << i + 2
             << setw(20) << fixed << setprecision(12) << roots[i + 1]
             << setw(20) << scientific << setprecision(6) << errors[i];

        if (i < (int)orders.size() && orders[i] != 0.0L)
            cout << setw(16) << fixed << setprecision(6) << orders[i];
        else
            cout << setw(16) << "---";

        cout << "\n";
    }
    cout << "\n";
}

// Export Function: Individual File Export
void exportIndividualCSV(const string& filename, 
                         const string& methodName,
                         const vector<long double>& errors, 
                         const vector<long double>& orders) {
    
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "[ERROR] Could not write to " << filename << "\n";
        return;
    }

    // Header specific to this method
    file << "iteration,error,order\n";

    for (int i = 0; i < (int)errors.size(); i++) {
        file << i + 1 << "," 
             << scientific << setprecision(15) << errors[i] << ",";
        
        if (i < (int)orders.size())
            file << fixed << setprecision(6) << orders[i] << "\n";
        else
            file << "0.000000\n"; // Last row where order can't be computed
    }

    file.close();
}

int main() {

    // Input files — roots logged from solvers
    string sec_ndd_file = "secant_ndd_roots.csv";
    string nwt_ndd_file = "newton_rp_ndd_roots.csv";
    string sec_chb_file = "secant_chebyshev_roots.csv";
    string nwt_chb_file = "newton_rp_chebyshev_roots.csv";

    // 1. Load
    vector<long double> sec_ndd_roots = readRootsFromCSV(sec_ndd_file);
    vector<long double> nwt_ndd_roots = readRootsFromCSV(nwt_ndd_file);
    vector<long double> sec_chb_roots = readRootsFromCSV(sec_chb_file);
    vector<long double> nwt_chb_roots = readRootsFromCSV(nwt_chb_file);

    // 2. Safety check
    if (sec_ndd_roots.empty() || nwt_ndd_roots.empty() ||
        sec_chb_roots.empty() || nwt_chb_roots.empty()) {
        cerr << "[!] Error: One or more root files are missing or empty.\n";
        return 1;
    }

    // 3. Analysis
    auto sec_ndd_err = successiveErrors(sec_ndd_roots);
    auto sec_ndd_ord = convergenceOrders(sec_ndd_err);

    auto nwt_ndd_err = successiveErrors(nwt_ndd_roots);
    auto nwt_ndd_ord = convergenceOrders(nwt_ndd_err);

    auto sec_chb_err = successiveErrors(sec_chb_roots);
    auto sec_chb_ord = convergenceOrders(sec_chb_err);

    auto nwt_chb_err = successiveErrors(nwt_chb_roots);
    auto nwt_chb_ord = convergenceOrders(nwt_chb_err);

    // 4. Console output
    printTable("SECANT - Newton DD P'(x)",         sec_ndd_roots, sec_ndd_err, sec_ndd_ord);
    printTable("NEWTON-RAPHSON - Newton DD P'(x)", nwt_ndd_roots, nwt_ndd_err, nwt_ndd_ord);
    printTable("SECANT - Chebyshev P'(x)",         sec_chb_roots, sec_chb_err, sec_chb_ord);
    printTable("NEWTON-RAPHSON - Chebyshev P'(x)", nwt_chb_roots, nwt_chb_err, nwt_chb_ord);

    // 5. Export
    exportIndividualCSV("secant_ndd_err_ord.csv",    "Secant",         sec_ndd_err, sec_ndd_ord);
    exportIndividualCSV("newton_rp_ndd_err_ord.csv", "Newton-Raphson", nwt_ndd_err, nwt_ndd_ord);
    exportIndividualCSV("secant_chebyshev_err_ord.csv",    "Secant",         sec_chb_err, sec_chb_ord);
    exportIndividualCSV("newton_rp_chebyshev_err_ord.csv", "Newton-Raphson", nwt_chb_err, nwt_chb_ord);

    return 0;
}