#ifndef QR_H
#define QR_H

#include<iostream>
#include<vector>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<algorithm>

#include "utils.h"

using namespace std;

long double T(int k, long double x){
    if(k == 0) return 1;
    if(k == 1) return x;

    long double T0 = 1, T1 = x, T2;
    for(int i = 2; i <= k; i++){
        T2 = 2*x*T1 - T0;
        T0 = T1;
        T1 = T2;
    }
    return T1;
}

void householderQR(vector<vector<long double>>& A, vector<long double>& b, int m, int n){

    for(int k = 0; k < n; k++){
        long double norm = 0;
        for(int i = k; i < m; i++){
            norm += A[i][k]*A[i][k];
        }
        norm = sqrtl(norm);

        if(norm == 0) continue;

        long double alpha = (A[k][k] >= 0) ? -norm : norm;

        vector<long double> v(m, 0);
        v[k] = A[k][k] - alpha;

        for(int i = k+1; i < m; i++){
            v[i] = A[i][k];
        }

        long double beta = 0;
        for(int i = k; i < m; i++){
            beta += v[i]*v[i];
        }

        if(beta == 0) continue;

        for(int j = k; j < n; j++){
            long double dot = 0;
            for(int i = k; i < m; i++){
                dot += v[i]*A[i][j];
            }
            dot *= (2.0L / beta);

            for(int i = k; i < m; i++){
                A[i][j] -= dot * v[i];
            }
        }

        long double dot_b = 0;
        for(int i = k; i < m; i++){
            dot_b += v[i]*b[i];
        }
        dot_b *= (2.0L / beta);

        for(int i = k; i < m; i++){
            b[i] -= dot_b * v[i];
        }

        A[k][k] = alpha;
        for(int i = k+1; i < m; i++){
            A[i][k] = 0;
        }
    }
}

vector<long double> solveR(const vector<vector<long double>>& R, const vector<long double>& b, int n){

    vector<long double> c(n, 0);

    for(int i = n-1; i >= 0; i--){
        long double sum = 0;
        for(int j = i+1; j < n; j++){
            sum += R[i][j] * c[j];
        }
        c[i] = (b[i] - sum) / R[i][i];
    }

    return c;
}

vector<long double> chebyshev_householder(vector<long double>& x, vector<long double>& y, int degree){

    int m = x.size();
    int n = degree + 1;

    vector<vector<long double>> A(m, vector<long double>(n));

    long double a = *min_element(x.begin(), x.end());
    long double b_val = *max_element(x.begin(), x.end());

    vector<long double> x_scaled(m);
    for(int i = 0; i < m; i++){
        x_scaled[i] = scale_to_chebyshev(x[i], a, b_val);
    }

    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            A[i][j] = T(j, x_scaled[i]);
        }
    }

    vector<long double> b = y;

    householderQR(A, b, m, n);

    vector<vector<long double>> R(n, vector<long double>(n));
    for(int i = 0; i < n; i++){
        for(int j = i; j < n; j++){
            R[i][j] = A[i][j];
        }
    }

    vector<long double> b1(n);
    for(int i=0;i<n;i++){
        b1[i] = b[i];
    }

    vector<long double> c = solveR(R, b1, n);

    vector<long double> standard = convert_to_standard(c, degree);

    // NEW STEP (convert x' → x)
    vector<long double> original = convert_to_original_x(standard, a, b_val);

    run_root_finding(original, x);

    cout << "\nResulting Polynomial:\nP(x) = ";

    long double EPS = 1e-15; // Filter out numerical noise
    bool first = true;

    // Loop backwards for descending order (i = degree down to 0)
    for (int i = original.size() - 1; i >= 0; i--) {
        
        // Skip terms that are effectively zero
        if (fabsl(original[i]) < EPS) continue;

        // Handle the sign between terms
        if (!first) {
            cout << (original[i] > 0 ? " + " : " - ");
        } else {
            // Special handling for the very first term
            if (original[i] < 0) cout << "-";
            first = false;
        }

        // Print the absolute value (magnitude) of the coefficient
        cout << fabsl(original[i]);

        // Print the power of x correctly
        if (i == 1) {
            cout << "x";
        } else if (i > 1) {
            cout << "x^" << i;
        }
    }
    cout << endl;

    ofstream out("qr_polynomial_coefficients.csv");
    out << fixed << setprecision(15);

    // out << a << "\n";
    // out << b_val << "\n";

    for(auto v : original){
        if(fabsl(v) < EPS) v = 0;
        out << v << "\n";
    }

    out.close();

    return c;
}

#endif