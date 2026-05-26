#ifndef CHEBY_GRADIENT_H
#define CHEBY_GRADIENT_H

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

    void Gradient_Descent(vector<long double>& x,vector<long double>& y,vector<long double>& c,int d,long double alpha,int max_iter,long double tol){
        int N = x.size();

        vector<vector<long double>> A(N, vector<long double>(d+1, 0));
        for(int i=0;i<N;i++){
            for(int k=0;k<d+1;k++){
                A[i][k]=T_x_k(k,x[i]);
            }
        }
        c.assign(d + 1, 0.0);

        cout << fixed << setprecision(8);
        cout << "\n--- Gradient Descent Log ---\n";
        cout << setw(8)  << "Iter"<< setw(18) << "Loss L(c)"<< setw(18) << "||gradient||" << "\n";
        cout << string(44, '-') << "\n";

        for (int iter = 0; iter <= max_iter; iter++) {

            vector<long double> residual(N,0);

            for(int i=0; i<N; i++) {
                long double s = 0.0;
                for(int j=0; j<d+1; j++){
                    s += c[j]*A[i][j];
                }
                residual[i] = s - y[i];
            }
            long double loss = 0;
            for (int i= 0; i<N; i++){
                loss += residual[i]*residual[i];
            }
            loss /= (2*N);

            vector<long double>grad(d+1,0);

            for(int i=0; i<d+1; i++){
                for(int j=0; j<N; j++){
                    grad[i] += A[j][i]*residual[j];
                }
            grad[i] /= N;
            }

            long double grad_norm = 0.0;
            for (int k = 0; k <= d; k++){
                grad_norm += grad[k] * grad[k];
            }
            grad_norm = sqrt(grad_norm);

            if (iter % 100 == 0)
                cout << setw(8)  << iter << setw(18) << loss << setw(18) << grad_norm << "\n";

            if (grad_norm<tol) {
                cout <<"\nConverged at iteration "<< iter<<"  (||grad|| = "<<grad_norm<<" < tol = "<<tol<< ")\n";
                break;
            }

            for (int k=0;k<=d;k++){
                c[k] -= alpha*grad[k];
            }
        }

        long double max_err=0, rms_err=0;
        for (int i=0; i<N; i++) {
            long double s = 0.0;
            for (int j=0; j<d+1; j++){
                s += c[j] * A[i][j];
            }
            long double e = fabsl(s-y[i]);
            max_err  = max(max_err, e);
            rms_err += e * e;
        }
        rms_err = sqrt(rms_err / N);
        
        cout<<endl;
        cout << "Max error : " << max_err  << "\n";
        cout << "RMS error : " << rms_err  << "\n";

    }

    vector<long double> Chebyshev_Using_Gradient_descent(vector<long double>& x, vector<long double>& y, int degree){
        int N=x.size();
        int d=degree;
        long double alpha = 0.1;
        int max_iter = 10000;
        long double tol = 1e-12;

        long double a = *min_element(x.begin(), x.end());
        long double b_val = *max_element(x.begin(), x.end());

        vector<long double> x_scaled(N);
        for(int i=0; i<N;i++){
            x_scaled[i] = scale_to_chebyshev(x[i], a, b_val);
        }
        vector<long double>c(d+1,0);

        Gradient_Descent(x_scaled,y,c,d,alpha,max_iter,tol);

        vector<long double> standard = convert_to_standard(c, d);
        vector<long double> original = convert_to_original_x(standard, a, b_val);


        run_root_finding(original, x);

        cout << "\nChebyshev Polynomial (Gadient Descent):\nP(x) = ";

        long double EPS = 1e-15;
        bool first = true;

        for(int i=original.size()-1; i>=0; i--){

            if(!first){
                if(original[i]>0) cout << " + ";
                else cout << " - ";
            } else {
                if(original[i]<0) cout << "-";
                first = false;
            }

            cout << fabsl(original[i]);
            if(i>=1) cout << "x";
            if(i>=2) cout << "^" << i;
        }

        cout << endl;

        ofstream out("Chebyshev_using_gradient_descent_polynomial_coefficients.csv");
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