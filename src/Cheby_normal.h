#ifndef CHEBY_NORMAL_H
#define CHEBY_NORMAL_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <functional>
#include <string>
#include<fstream>
#include<algorithm>
#include "NR_Secant.h"
#include "utils.h"

using namespace std;

void Chebyshev_Polynomial_Approximation(vector<long double>& x,vector<long double>& y,vector<long double>& c,int d){
    int N = x.size();
    vector<vector<long double>> A(N, vector<long double>(d+1, 0));
    for(int i=0;i<N;i++){
        for(int k=0;k<d+1;k++){
            A[i][k]=T_x_k(k,x[i]);
        }
    }

    vector<vector<long double>>AtA(d+1, vector<long double>(d+1, 0));
    for(int k=0;k<d+1;k++){
        for(int j=0;j<d+1;j++){
            for(int i=0;i<N;i++){
                AtA[k][j] += A[i][j]*A[i][k];
            }
        }
    }

    vector<vector<long double>>Aty(d+1,vector<long double>(1,0));
    for(int i=0;i<d+1;i++){
        for(int j=0;j<N;j++){
            Aty[i][0] += A[j][i]*y[j];
        }
    }

    int n = d + 1;
    vector<vector<long double>> M(n, vector<long double>(n+1));
    for (int i =0; i<n; i++) {
        for (int j=0; j<n; j++)
            M[i][j] = AtA[i][j];
        M[i][n] = Aty[i][0];
    }

    for (int col=0; col<n; col++){
        int pivot = col;
        for (int row=col+1; row<n; row++)
            if (abs(M[row][col]) > abs(M[pivot][col]))
                pivot = row;
        swap(M[col], M[pivot]);

        for (int row=col+1; row<n; row++) {
            long double factor = M[row][col]/M[col][col];
            for (int j=col; j<=n; j++)
                M[row][j] -= factor*M[col][j];
        } 
    }

    for (int i=n-1; i>=0; i--) {
        c[i] = M[i][n];
        for (int j=i+1; j<n; j++){
            c[i] -= M[i][j]*c[j];   
        }
        c[i] /= M[i][i];
    }
}

vector<long double> Chebyshev_Using_Normal_Equation(vector<long double>& x, vector<long double>& y, int degree){
    
    int N=x.size();
    int d=degree;
    long double a = *min_element(x.begin(), x.end());
    long double b_val = *max_element(x.begin(), x.end());

    vector<long double> x_scaled(N);
    for(int i=0; i<N;i++){
        x_scaled[i] = scale_to_chebyshev(x[i], a, b_val);
    }
    vector<long double>c(d+1,0);

    Chebyshev_Polynomial_Approximation(x_scaled,y,c,d);

    vector<long double> standard = convert_to_standard(c, d);
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

    ofstream out("Chebyshev_using_normal_equation_polynomial_coefficients.csv");
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

