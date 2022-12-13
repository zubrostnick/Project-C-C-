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

//Задання еліпса у пдск, осі якого парлельні осям координат
class Ellipse {
/*canonical ellipse eqation (x-x0)^2/a^2 + (y-y0)^2/b^2 = 1 
wtih additional methods*/ 
    double a, b, x0, y0;
    
    void init(Point f1, Point f2, double r1, double r2);

public:

    //creating Ellipse with given paremeters using two radiuses and two focuses
    Ellipse(
        Point f1,
        Point f2,
        double r1,
        double r2
    );

    //default
    Ellipse();
    
    //destructor
    ~Ellipse();

    //copy-constructor: create Ellipse from given
    Ellipse(const Ellipse& e);

    //create and input Ellipse
    static Ellipse create_input();

    void input();
    
    //output of ellipse 
    virtual void print(ostream& out);


    friend ostream& operator<<(ostream& out, const Ellipse& e);

    //perfomes direct comparison between two Ellipses 
    bool equal(const Ellipse& other);
    bool operator==(const Ellipse& other);
    bool operator!=(const Ellipse& other);


    //block of setters and getters
    double get_a();
    void set_a(double _a);

    double get_b();
    void set_b(double _b);

    double get_x0();
    void set_x0(double _x0);

    double get_y0();
    void set_y0(double _y0);

    Point get_f1();
    Point get_f2();

    double area();

    //getting center of ellipse
    Point get_center();

    //checks whether point lies inside given ellipse
    bool is_inside_ellipse(Point p);
    
    //Intersection points of ellipse with given line, given by equation y = kx + h
    void intersect_line(double k, double h, ostream& out);
    
    //Intersection points of ellipses (found numerically, rounded off)
    void intersect_ellipse(Ellipse e, ostream& out);
};
