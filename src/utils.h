#ifndef UTILS_H
#define UTILS_H

#include<iostream>
#include<vector>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<algorithm>

using namespace std;

long double scale_to_chebyshev(long double x, long double a, long double b){
    return (2.0*(x - a)/(b - a)) - 1.0;
}

vector<long double> convert_to_standard(const vector<long double>& c, int d){

    vector<vector<long double>> Tpoly(d+1);

    Tpoly[0] = {1};
    if(d >= 1) Tpoly[1] = {0,1};

    for(int k = 2; k <= d; k++){
        vector<long double> temp(k+1, 0);

        for(int i = 0; i < Tpoly[k-1].size(); i++){
            temp[i+1] += 2 * Tpoly[k-1][i];
        }

        for(int i = 0; i < Tpoly[k-2].size(); i++){
            temp[i] -= Tpoly[k-2][i];
        }

        Tpoly[k] = temp;
    }

    vector<long double> result(d+1, 0);

    for(int k = 0; k <= d; k++){
        for(int i = 0; i < Tpoly[k].size(); i++){
            result[i] += c[k] * Tpoly[k][i];
        }
    }

    return result;
}

vector<long double> convert_to_original_x(
    const vector<long double>& coeff,
    long double a,
    long double b
){
    int n = coeff.size();

    vector<long double> result(n, 0);

    long double m = 2.0/(b - a);
    long double c = - (b + a)/(b - a);

    vector<long double> power = {1};

    for(int k = 0; k < n; k++){

        for(int i = 0; i < power.size(); i++){
            result[i] += coeff[k] * power[i];
        }

        vector<long double> new_power(power.size()+1, 0);

        for(int i = 0; i < power.size(); i++){
            new_power[i] += c * power[i];
            new_power[i+1] += m * power[i];
        }

        power = new_power;
    }

    return result;
}

long double T_x_k(int k,long double x){
        if(k==0){
            return 1;
        }
        else if(k==1){
            return x;
        }
        long double T_n,T_n_1,T_n_2;
        T_n_2=1;
        T_n_1=x;
        for(int i=2;i<=k;i++){
            T_n = 2*x*T_n_1 - T_n_2;
            T_n_2 = T_n_1;
            T_n_1 = T_n;
        }
        return T_n;
    }

#endif