#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <cstring> //memmove, memcpy
#include <istream>
#include <ostream>
#include <sstream>
#include <cmath>

using namespace std;
#define _USE_MATH_DEFINES


struct Point{
/*Structure that creates 2d point + operans reassigment*/
    double x, y;
    
    // main constructor    
    Point(double _x, double _y);
   
    //default
    Point();

    virtual void input(istream& inp);
    
    virtual void print(ostream& out);

    friend ostream& operator<<(ostream& out,const Point& p);

    friend istream& operator>>(istream& inp, Point& p);

    //creates and returns a point 
    static Point create_input();
    /*Create and input point, writing in x, y coordinates of a point */
};
