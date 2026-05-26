#include<iostream>
#include<vector>
#include<fstream>
#include<iomanip>
#include<algorithm>

#include "utils.h"
#include "newton.h"
#include "Cheby_normal.h"
#include "qr.h"
#include "evaluate.h"
#include "Cheby_gradient.h"
#include "NR_Secant.h"

using namespace std;

int main(){

    ifstream file("data.csv");

    vector<long double> x, y;
    long double x_val, y_val;
    char comma;

    file.ignore(1000, '\n');

    while(file >> x_val >> comma >> y_val){
        x.emplace_back(x_val);
        y.emplace_back(y_val);
    }

    file.close();

    cout << "\n=============================\n";
    cout << "Newton Polynomial\n";
    cout << "=============================\n";

    vector<long double> newton_coeff = newton_divided_difference(x,y);

    cout << "\n=============================\n";
    cout << "QR Chebyshev Polynomial\n";
    cout << "=============================\n";

    int degree = 12;

    int choice;
    cout << "Choose Chebyshev Method:\n";
    cout << "1. QR Method\n";
    cout << "2. Normal Equations\n";
    cout << "3. Gradient Descent\n";
    cout << "Enter Your choice: ";
    cin >> choice;

    vector<long double> cheb_coeff;

    if(choice == 1){
        cheb_coeff = chebyshev_householder(x, y, degree);
    }
    else if(choice == 2){
        cheb_coeff = Chebyshev_Using_Normal_Equation(x, y, degree);
    }
    else if(choice == 3){
        cheb_coeff = Chebyshev_Using_Gradient_descent(x, y, degree);
    }
    else{
        cout << "Invalid choice\n";
        return 1;
    }
    
    long double a = *min_element(x.begin(), x.end());
    long double b = *max_element(x.begin(), x.end());

    ofstream out("comparison.csv");

    out << fixed << setprecision(15);
    out << "x,Newton,Chebyshev\n";

    int N = 600;

    for(int i = 0; i <N; i++){
        long double xi = a + (b - a) * i / N;

        long double y_newton = evaluate_newton(xi, x, newton_coeff);
        long double y_cheb = evaluate_chebyshev(xi, cheb_coeff, a, b);

        out << xi << "," << y_newton << "," << y_cheb << "\n";
    }

    out.close();

    cout << "\nComparison file generated: comparison.csv\n";

    cout << "\nOpening plot...\n";


    ofstream out1("comparision_25.csv");
    out1 << fixed << setprecision(15);

    int N1 = 25;
    for(int i = 0; i < N1; i++){
        long double xi = x[i];

        long double y_newton_1 = evaluate_newton(xi, x, newton_coeff);
        long double y_cheb_1 = evaluate_chebyshev(xi, cheb_coeff, a, b);

        out1 << xi << "," << y_newton_1 << "," << y_cheb_1 << "\n";
    }

    out1.close();


    #ifdef _WIN32
        system("python interpolation_gradient.py");
    #else
        system("arch -arm64 python3 interpolation_gradient.py");
    #endif


    return 0;
}