#ifndef EVALUATE_H
#define EVALUATE_H

#include <vector>

using namespace std;

// Evaluates the Chebyshev series sum(c_k * T_k(x_scaled)) using Clenshaw's Recurrence
long double evaluate_chebyshev(long double x, const vector<long double>& c, long double a, long double b) {
    // 1. Map original x back to [-1, 1]
    long double x_scaled = (2.0L * (x - a) / (b - a)) - 1.0L;

    // 2. Clenshaw's Algorithm
    int n = c.size() - 1;
    long double b2 = 0.0L;
    long double b1 = 0.0L;
    
    for (int k = n; k >= 1; k--) {
        long double b0 = 2.0L * x_scaled * b1 - b2 + c[k];
        b2 = b1;
        b1 = b0;
    }
    
    return x_scaled * b1 - b2 + c[0];
}

// Keep the old one for Newton (which is already in power basis)
long double evaluate_standard(long double x, const vector<long double>& coeff) {
    int n = coeff.size();
    long double result = coeff[n-1];
    for (int i = n-2; i >= 0; i--) {
        result = result * x + coeff[i];
    }
    return result;
}

long double evaluate_newton(long double xi, const vector<long double>& x_pts, const vector<long double>& coeffs) {
    int n = coeffs.size();
    long double result = coeffs[n - 1]; 
    for (int i = n - 2; i >= 0; i--) {
        result = result * (xi - x_pts[i]) + coeffs[i]; // Nested multiplication
    }
    return result;
}

#endif