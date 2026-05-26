#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>

using namespace std;

// Analytical target for True Error calculation
const double ANALYTICAL_ROOT = M_PI / 10.0; 

// --- PART 1: DATA INGESTION ---
// Reads coefficients c0, c1, c2... from a single-column CSV
vector<double> readCoefficientsCSV(string filename) {
    vector<double> coeffs;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Error: Could not open " << filename << ". Ensure the file exists." << endl;
        return coeffs;
    }

    string line;
    while (getline(file, line)) {
        if (line.empty()) continue;
        try {
            coeffs.push_back(stod(line));
        } catch (...) {
            continue; // Skip headers or non-numeric lines
        }
    }
    file.close();
    return coeffs;
}

// --- PART 2: POLYNOMIAL TOOLS ---
double evaluatePoly(double x, const vector<double>& coeffs) {
    double result = 0;
    double x_pow = 1.0;
    for (double c : coeffs) {
        result += c * x_pow;
        x_pow *= x;
    }
    return result;
}

// --- PART 3: OUTPUT UTILITY ---
void saveToCSV(string filename, const vector<double>& history) {
    ofstream file(filename);
    for (double val : history) {
        file << fixed << setprecision(10) << val << endl;
    }
    file.close();
}

// --- PART 4: ROOT FINDING ---

void solveSecant(double x0, double x1, const vector<double>& derivCoeffs) {
    cout << "\n--- Secant Method (Solving P'(x)=0) ---" << endl;
    vector<double> x_history;
    x_history.push_back(x0);
    x_history.push_back(x1);

    double x_p = x0, x_c = x1;
    for (int i = 0; i < 30; i++) {
        double f_p = evaluatePoly(x_p, derivCoeffs);
        double f_c = evaluatePoly(x_c, derivCoeffs);
        if (abs(f_c - f_p) < 1e-18) break;

        double x_next = x_c - f_c * (x_c - x_p) / (f_c - f_p);
        x_history.push_back(x_next);
        
        if (abs(x_next - x_c) < 1e-12) break;
        x_p = x_c; x_c = x_next;
    }

    // Console output: 4 columns
    cout << left << setw(6) << "Iter" << setw(18) << "x value" << setw(18) << "Step Error" << setw(18) << "True Error" << endl;
    cout << string(60, '-') << endl;
    for (size_t i = 0; i < x_history.size(); i++) {
        cout << left << setw(6) << i << setw(18) << fixed << setprecision(10) << x_history[i];
        if (i == 0) cout << setw(18) << "---";
        else cout << setw(18) << abs(x_history[i] - x_history[i-1]);
        cout << setw(18) << abs(ANALYTICAL_ROOT - x_history[i]) << endl;
    }

    saveToCSV("secant_output.csv", x_history);
}

void solveNewtonRaphson(double x_guess, const vector<double>& derivCoeffs) {
    cout << "\n--- Newton-Raphson Method (Solving P'(x)=0) ---" << endl;
    vector<double> secondDeriv;
    for (size_t i = 1; i < derivCoeffs.size(); i++) {
        secondDeriv.push_back(i * derivCoeffs[i]);
    }

    vector<double> x_history;
    x_history.push_back(x_guess);
    double x = x_guess;
    for (int i = 0; i < 30; i++) {
        double fx = evaluatePoly(x, derivCoeffs);      
        double dfx = evaluatePoly(x, secondDeriv);    
        if (abs(dfx) < 1e-18) break;

        double x_next = x - fx / dfx; 
        x_history.push_back(x_next);

        if (abs(x_next - x) < 1e-12) break; 
        x = x_next;
    }

    // Console output: 4 columns
    cout << left << setw(6) << "Iter" << setw(18) << "x value" << setw(18) << "Step Error" << setw(18) << "True Error" << endl;
    cout << string(60, '-') << endl;
    for (size_t i = 0; i < x_history.size(); i++) {
        cout << left << setw(6) << i << setw(18) << fixed << setprecision(10) << x_history[i];
        if (i == 0) cout << setw(18) << "---";
        else cout << setw(18) << abs(x_history[i] - x_history[i-1]);
        cout << setw(18) << abs(ANALYTICAL_ROOT - x_history[i]) << endl;
    }

    saveToCSV("newton_output.csv", x_history);
}

int main() {
    // 1. Read p(x) coefficients from CSV
    vector<double> pCoeffs = readCoefficientsCSV("polynomial_coefficients.csv");
    if (pCoeffs.empty()) return 1;

    // 2. Compute p'(x) coefficients
    // If p(x) = c0 + c1*x + c2*x^2 ..., then p'(x) = c1 + 2*c2*x + 3*c3*x^2 ...
    vector<double> derivCoeffs;
    for (size_t i = 1; i < pCoeffs.size(); i++) {
        derivCoeffs.push_back(i * pCoeffs[i]);
    }

    // 3. Solve for p'(x) = 0
    solveSecant(0.1, 0.2, derivCoeffs);
    solveNewtonRaphson(0.15, derivCoeffs);

    return 0;
}