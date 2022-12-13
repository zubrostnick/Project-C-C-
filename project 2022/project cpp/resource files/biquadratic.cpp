#include "biquadratic.h"

Biquadratic_eq::Biquadratic_eq(
    double _a2 = 0,
    double _a1 = 0,
    double _a0 = 0
):  Quadratic_eq(_a2, _a1, _a0) {}

void Biquadratic_eq::print(ostream& out = cout) {
    out <<""<< get_a2() <<"*x^4 + (" << get_a1() <<")*x^2 + (" << get_a0() <<") = 0 ";
}

void Biquadratic_eq::solve(ostream& out = cout){
    ostringstream res1;
    ostringstream o_res;
    int num = 0; //amount of roots 
    
    Quadratic_eq::solve(res1);

    istringstream i_res(res1.str());
    double val;
    while ((i_res >> val)) {
        if (val >= 0){
            num++;
            out << sqrt(val)<<endl << -sqrt(val) << endl;
        }                      
    }
    if (res1.str().find("nan") != string::npos || res1.str().find("inf") != string::npos) {
        out <<res1.str() << ends;
        // std::cout << "found inf or nan!" << '\n';
    } else {
        if (num == 0){  //(!num)
            out << NAN << endl;
        }  
    }  
    
}  
