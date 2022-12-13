#pragma once
#include "linear.h"
#include "quadratic.h"
#include "biquadratic.h"

#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;


class Equations_from_file{

    string filename;

public:

    void read_monomial(string monom, ostream& out);

    Equations_from_file(const string& filename);

    void exercise ();
};