#pragma once

#include <stdexcept>
#include <typeinfo>
#include <exception>
#include <iostream>
using namespace std;

//Exeption class for working on curves 
class ConicError: public logic_error {
    public: ConicError(const string& message);
};
