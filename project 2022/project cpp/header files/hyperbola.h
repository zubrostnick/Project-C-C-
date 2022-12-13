#pragma once

#include "point.h"
#include "conic_error.h"
#include "linear.h"
#include "quadratic.h"
#include "biquadratic.h"
#include "from_file.h"
#include "parabola.h"

#include <iostream>
#include <string>
#include <fstream>
#include <cstring> //memmove, memcpy
#include <istream>
#include <ostream>
#include <sstream>
#include <cmath>

class Hyperbola {
/*using inverted standard formula for hyperbola: (y-y0)^2/b^2 - (x-x0)^2/a^2 = 1
Hyperbola axis is parallel to Oy-axis*/

    double *v; //vector, that contains a, b, x0, y0

public:

    Hyperbola(double a, double b, double x0, double y0); 
    /*creating parabola with given paremeters*/ 
        

    //default constructor
    Hyperbola();

    Hyperbola(double *ptr_v);
    /*Hyperbola from given array of length 4 with given parameters*/   
    
    ~Hyperbola();
    //destructor

    Hyperbola(const Hyperbola& a);
    /*copy-constructor: create Hyperbola from given*/
  

    static Hyperbola create_input();
    /*create and input Hyperbola */


    void input();
    /* input a, b, x0, y0 parameters of hyperbola (y-y0)^2/b^2 - (x-x0)^2/a^2 = 1 */
    

    //block of setters and getters
    double get_a();
    void set_a(double _a);

    double get_b();
    void set_b(double _b);

    double get_x0();
    void set_x0(double _x0);

    double get_y0();
    void set_y0(double _y0);  

    bool equal(const Hyperbola& other);
    /*perfomes direct comparison between two Ellipses*/ 


    virtual void print(ostream& out); 
    /*output of hyperbola*/

    friend ostream& operator<<(ostream& out, const Hyperbola& h);
    /*operator realoading for output*/


    Point get_center();
    /*getting center of hyperbola*/
  

    bool is_outside_hyperbola(Point p);
    /*checks whether point lies outside given hyperbola
    (above upper part, below lower part)*/

    
    void intersect_parabola(Parabola p, ostream& out);
    /*Intersection points of Hyperbola with given Parabola
    Intersection points (found numerically, rounded off)*/
        

    void intersect_hyperbola(Hyperbola other, ostream& out);
    /*Intersection points of Hyperbola with given Hyperbola
    Intersection points (found numerically, rounded off)*/
};