#pragma once

#include "point.h"
#include "conic_error.h"
#include "linear.h"
#include "quadratic.h"
#include "biquadratic.h"
#include "from_file.h"

#include <iostream>
#include <string>
#include <fstream>
#include <cstring> //memmove, memcpy
#include <istream>
#include <ostream>
#include <sstream>
#include <cmath>

class Parabola {
/*using unverted standard formula for parabola: (x-x0)^2 = 2p(y-y0)
Parabola axis is parallel to Oy-axis*/

    double *v; //vector, that contains p, x0, y0

public:

    Parabola( double p, double x0, double y0);
    /*creating parabola with given paremeters*/ 

    
    //default
    Parabola();

    Parabola(double *ptr_v);
    /*Parabola from given array of length 3 with given parameters*/
    
    //destructor
    ~Parabola();


    Parabola(const Parabola& a);
    /*copy-constructor: create Parabola from given*/


    static Parabola create_input();
    /*create and input Parabola */


    void input();
    /*Enter p, x0, y0 parameters of parabola (x-x0)^2 = 2p(y-y0) */
    

    //block of setters and getters
    double get_p();
    void set_p(double _p);

    double get_x0();
    void set_x0(double _x0);

    double get_y0();
    void set_y0(double _y0);  


    virtual void print(ostream& out);
    /*output of parabola*/

    friend ostream& operator<<(ostream& out, const Parabola& p);
    /*operator realoading for output*/


    Point get_center();
    /*getting center of parabola (focus)*/

    
    bool is_above_parabola(Point p);
    /*checks whether point lies above given parabola*/
    
    
    void intersect_line( double k, double b, ostream& out);
    /*Intersection points of Parabola with given line, given by equation y = kx + b*/


    double intersect_line_area( double k, double b);
    /*Intersection area of Parabola with given line, given by equation y = kx + b*/
    
    
    void intersect_parabola(Parabola p, ostream& out);
    /*Intersection points of Parabola with given Parabola*/


    double intersect_parabola_area(Parabola p);
    /*Intersection area of Parabola with given Parabola*/
};
