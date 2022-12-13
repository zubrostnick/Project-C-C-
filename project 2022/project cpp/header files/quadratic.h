#pragma once
#include "linear.h"
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;

class Quadratic_eq: public Linear_eq {
    //a2*x^2 + a1 * x + a0 = 0

    double a2;

public:

    Quadratic_eq(
        double _a2,
        double _a1,
        double _a0
    );
    double get_a2();
    void set_a2(double _a2);

    virtual void input(istream& inp);

    virtual void print(ostream& out);

    //output of solvings in stream 
    void solve(ostream& out);

};