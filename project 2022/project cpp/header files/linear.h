#pragma once

#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;

class Linear_eq {
    //a1 * x + a0 = 0
    double a0;
    double a1;

public:

    Linear_eq(double _a1, double _a0);

    double get_a0();
    void set_a0(double a);     
    double get_a1();
    void set_a1(double a);

    //virtual: Буде викликано в залежності від поточного стану програми
    virtual void input(istream& inp);
    virtual void print(ostream& out);

    //Перевантаження оператору в загальному розумінні (як зовнішню функцію, а не внутрішній метод)
    friend ostream& operator<<(ostream& out, const Linear_eq& eq);

    friend istream& operator>>(istream& inp, Linear_eq& eq);

    //output of solvings in stream 
    void solve(ostream& out);
};