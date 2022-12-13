//Linear describtion
#include "linear.h"

Linear_eq::Linear_eq(double _a1 = 0, double _a0 = 0): a0(_a0), a1(_a1) {}

double Linear_eq::get_a0() {return a0;}
void Linear_eq::set_a0(double a) {this->a0 = a;}
double Linear_eq::get_a1() {return a1;}
void Linear_eq::set_a1(double a) {this->a1 = a;}

//virtual: Буде викликано в залежності від поточного стану програми
void Linear_eq::input(istream& inp = cin) {
    inp >> a1 >> a0;
}

void Linear_eq::print(ostream& out = cout) {
    out <<"("<< a1 << ")*x + (" << a0 << ") = 0";
}

//Перевантаження оператору в загальному розумінні (як зовнішню функцію, а не внутрішній метод)
ostream& operator<<(ostream& out, Linear_eq& eq) {
    eq.print(out);
    return out;
}

istream& operator>>(istream& inp, Linear_eq& eq) {
    eq.input(inp);
    return inp;
}

//output of solvings in stream 
void Linear_eq::solve(ostream& out = cout){
    //cout << "check" << endl;
    if ((a1 == 0) && (a0 != 0)){
        //throw logic_error("Unsolvable");
        out << NAN << endl;
    } else if ((a1 == 0) && (a0 == 0)){
        out << INFINITY << endl;
    } else out << -a0/a1 << endl;
}
