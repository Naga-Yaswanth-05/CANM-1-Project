#include<iostream>
#include<vector>
#include<cmath>
#include<iomanip>

using namespace std;

vector<long double> forward_difference(vector<long double>& x, vector<long double>& y){

    int n = x.size();

    long double h = x[1] - x[0];

    vector<vector<long double>> table(n, vector<long double>(n));

    //first column = y values
    for(int i=0; i<n;i++){
        table[i][0] = y[i];
    }

    //filling forward difference table
    for(int j=1; j<n;j++){
        for(int i=0; i<n-j;i++){
            table[i][j] = table[i+1][j-1] - table[i][j-1];
        }
    }

    //converting forward difference polynomial directly to standard polynomial
    vector<long double> coeff_standard(n, 0);
    vector<long double> term = {1}; // represents current product term

    for(int i=0;i<n;i++){

        long double fact = 1;
        for(int k=1;k<=i;k++){
            fact *= k;
        }

        long double coeff_i = table[0][i] / fact;

        long double scale = coeff_i / pow(h, i);

        for(int j=0;j<term.size() && j<n;j++){
            coeff_standard[j] += scale * term[j];
        }

        if(i == n-1){
            break;
        }

        //multiply term by (x - (x0 + i*h))
        vector<long double> new_term(term.size()+1, 0);

        for(int j=0;j<term.size();j++){
            new_term[j] -= term[j] * (x[0] + i*h);
            new_term[j+1] += term[j];
        }

        term = new_term;
    }

    //printing standard polynomial
    cout << "\n Polynomial in standard form"<< endl;
    cout << "\n P(x) = ";

    long double EPS = 1e-12;
    bool first = true;

    for(int i=0; i<n; i++){

        if(fabsl(coeff_standard[i]) < EPS){
            coeff_standard[i] = 0;
            continue;
        }

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

    return coeff_standard;
}