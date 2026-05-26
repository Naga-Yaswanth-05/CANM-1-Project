#ifndef NR_SECANT_H
#define NR_SECANT_H

#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<algorithm>

using namespace std;

long double evaluate(long double c1, const vector<long double>& coeff) {
    long double total = 0;
    for (int i = 0; i < coeff.size(); i++) {
        total += coeff[i] * pow(c1, i);
    }
    return total;
}

// Secant Method
long double solveSecant(long double x0, long double x1, const vector<long double>& d_Coef) {

    long double x_prev = x0;
    long double x_curr = x1;

    for (int i = 0; i < 30; i++) {
        long double f_prev = evaluate(x_prev, d_Coef);
        long double f_curr = evaluate(x_curr, d_Coef);

        if (fabsl(f_curr - f_prev) < 1e-18) break;
        long double x = f_curr * (x_curr - x_prev) / (f_curr - f_prev);
        long double x_next = x_curr - x;

        if (fabsl(x_next - x_curr) < 1e-12)
            return x_next;

        x_prev = x_curr;
        x_curr = x_next;
    }

    return x_curr;
}

// Newton-Raphson Method
long double solveNewtonRaphson(long double x_guess, const vector<long double>& d_Coef) {

    vector<long double> secondDerivCoeffs;
    for (int i = 1; i < d_Coef.size(); i++) {
        secondDerivCoeffs.push_back(i*d_Coef[i]);
    }

    long double x = x_guess;

    for (int i = 0; i < 30; i++) {
        long double fx = evaluate(x, d_Coef);
        long double dfx = evaluate(x, secondDerivCoeffs);

        if (fabsl(dfx) < 1e-10) break;

        long double x_next = x - fx / dfx;

        if (fabsl(x_next - x) < 1e-12)
            return x_next;

        x = x_next;
    }

    return x;
}

//scan domain and find ALL roots
vector<long double> findAllRoots(long double a,long double b,int steps,const vector<long double>& d_Coef,string method) {
    vector<long double> roots;
    long double step = (b - a) / steps;

    for (int i = 0; i < steps; i++) {

        long double x1 = a + i * step;
        long double x2 = x1 + step;

        long double f1 = evaluate(x1, d_Coef);
        long double f2 = evaluate(x2, d_Coef);

        // sign change → root exists
        if (f1 * f2 <= 0) {

             long double root;

            if (method == "secant") {
                root = solveSecant(x1, x2, d_Coef);
            } else { // newton
                 long double mid = (x1 + x2) / 2;
                root = solveNewtonRaphson(mid, d_Coef);
            }

            // keep only roots inside domain
            if (root < a || root > b) continue;

            // avoid duplicates
            int c = 1;
            for (auto r : roots) {
                if (fabsl(r - root) < 1e-6) {
                    c = 0;
                    break;
                }
            }
            if (c) {
                roots.push_back(root);
            }
            
        }
    }

    return roots;
}

void run_root_finding(vector<long double>& coeff_standard, vector<long double>& x){

    vector<long double> derivCoeffs;
    for (int i = 1; i < coeff_standard.size(); i++) {
        derivCoeffs.push_back(i * coeff_standard[i]);
    }

    long double a = *min_element(x.begin(), x.end());
    long double b = *max_element(x.begin(), x.end());

    if (coeff_standard.size() > 2) {

        string method = "newton";

        vector<long double> roots = findAllRoots(a, b, 200, derivCoeffs, method);

        cout << fixed << setprecision(10);
        cout << "\nRoots using " << method << " method in [" << a << "," << b << "]:\n\n";

        for (int i = 0; i < roots.size(); i++) {
            cout << "Root " << i+1 << " = " << setprecision(16)<< roots[i] << endl;
        }

        ofstream out("Newton_Raphson_Roots.csv");

        for(int i = 0; i < roots.size(); i++){
            out << setprecision(16)<< roots[i] << "\n";
        }

        out.close();
    }

    if (coeff_standard.size() > 2) {

        string method = "secant";

        vector<long double> roots = findAllRoots(a, b, 200, derivCoeffs, method);

        cout << fixed << setprecision(10);
        cout << "\nRoots using " << method << " method in [" << a << "," << b << "]:\n\n";

        for (int i = 0; i < roots.size(); i++) {
            cout << "Root " << i+1 << " = " << setprecision(16)<< roots[i]<< endl;
        }

        ofstream out("Secant_Roots.csv");

        for(int i = 0; i < roots.size(); i++){
            out << setprecision(16)<< roots[i] << "\n";
        }

        out.close();
    }

    cout<<endl;
}

#endif