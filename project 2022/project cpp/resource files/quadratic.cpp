#include "quadratic.h"
//quadratic equation class describtion

Quadratic_eq::Quadratic_eq(
    double _a2 = 0,
    double _a1 = 0,
    double _a0 = 0
):  Linear_eq(_a1, _a0), a2(_a2) {}

double Quadratic_eq::get_a2() {return a2;}
void Quadratic_eq::set_a2(double _a2) {this->a2 = _a2;}

void Quadratic_eq::input(istream& inp = cin) {
    inp >> a2;
    Linear_eq::input(inp);
}

void Quadratic_eq::print(ostream& out = cout) {
    out << a2 <<"*x^2 + ";
    Linear_eq::print(out);
}

//output of solvings in stream 
void Quadratic_eq::solve(ostream& out = cout){
    double res;
    if (a2 == 0){
        ostringstream check;
        Linear_eq :: solve (out);
        //Linear_eq :: solve(check);
        //cout << check.str() << endl;     
    }
    else if (get_a1() == 0) {
        if (-get_a0()/a2 < 0) {
            out << NAN << endl;
        }
        else out << sqrt(-get_a0()/a2) << endl << -sqrt(-get_a0()/a2) << endl;
    }
    else if (get_a0() == 0){
        out << 0 << endl;
        Linear_eq eq_q(a2, get_a1());
        eq_q.solve(out); // Але тут з відповіддю 0, може стояти inf, або nan 
    }
    else{
        
        double D = (get_a1()*get_a1() - 4*a2*get_a0());
        if (D < 0){
            out << NAN <<endl;
        } else if (D == 0){
            out << -get_a1()/2/a2 << endl;
        }
        else {
            out << (-get_a1() + sqrt(D))/2/a2 << endl << (-get_a1() - sqrt(D))/2/a2 <<endl;
        }
    }
    
} 
