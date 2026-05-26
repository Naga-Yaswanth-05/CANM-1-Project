#ifndef NEWTON_H
#define NEWTON_H

#include<iostream>
#include<vector>
#include<cmath>
#include <fstream>

#include "NR_Secant.h"

using namespace std;

vector<long double> newton_divided_difference(vector<long double>& x, vector<long double>& y){

    int n = x.size();

    vector<vector<long double>> table(n, vector<long double>(n));

    for(int i=0;i<n;i++){
        table[i][0] = y[i];
    }

    for(int j=1; j<n; j++){
        for(int i =0; i<n-j;i++){
            table [i][j] = (table[i+1][j-1] - table[i][j-1]) / (x[i+j] - x[i]);
        }
    }

    vector<long double> coeff(n);
    for(int j=0; j<n;j++){
        coeff[j] = table[0][j];
    }

    vector<long double> coeff_standard(n, 0);
    vector<long double> term = {1};

    for(int i=0;i<n;i++){
        for(int j=0;j<term.size() && j<n;j++){
            coeff_standard[j] += coeff[i] * term[j];
        }

        if(i == n-1){
            break;
        }

        vector<long double> new_term (term.size()+1, 0);

        for(int j=0;j<term.size(); j++){
            new_term[j] -= term[j] *x[i];
            new_term[j+1] += term[j];
        }

        term = new_term;
    }

    run_root_finding(coeff_standard, x);

    cout << "\n divided difference Polynomial in standard form"<< endl;
    cout << "\n P(x) = ";
    
    long double EPS = 1e-12;
    bool first = true;

    for(int i = 0; i < n; i++) {

        if(!first){
            if(coeff_standard[i] > 0){
                cout << " + ";
            }
            else{
                cout << " - ";
            }
        } 
        else{
            if(coeff_standard[i] < 0){
                cout << "-";
            }
            first = false;
        }

        cout << fabsl(coeff_standard[i]);

        if(i >= 1) cout << "x";
        if(i >= 2) cout << "^" << i;
    }

    cout << endl;

    ofstream out("polynomial_coefficients.csv");
    out << fixed << setprecision(15);

    for(int i = 0; i < n; i++){
        out << coeff_standard[i] << "\n";
    }

    out.close();
    return coeff;
}

#endif