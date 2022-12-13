#pragma once
#include "linear.h"
#include "quadratic.h"

#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;


class Biquadratic_eq: public Quadratic_eq {

public:

    Biquadratic_eq(
        double _a2,
        double _a1,
        double _a0
    );

    virtual void print(ostream& out);

    void solve(ostream& out);

};

