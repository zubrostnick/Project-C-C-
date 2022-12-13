#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>


#ifdef NAN
/* NAN is supported */
#endif
#ifdef INFINITY
/* INFINITY is supported */
#endif


#define _USE_MATH_DEFINES



typedef struct{
    //a2*x^2 + a1 * x + a0 = 0
    double a2;
    double a1;
    double a0;
}Quadratic_eq ;


Quadratic_eq init_quadratic(double _a2, double _a1, double _a0 );



double* solve(Quadratic_eq eq);


typedef struct{
/*Structure that creates 2d point*/
    double x, y;
} Point;

void print_point(Point p);
/*point output*/

//Задання еліпса у пдск, осі якого парлельні осям координат
typedef struct {
//canonical ellipse eqation (x-x0)^2/a^2 + (y-y0)^2/b^2 = 1
    double a, b, x0, y0;

}Ellipse;


//creating Ellipse with given paremeters задання еліпса координатами фокусів еліпса та 2 радіусами    
Ellipse init_ellipse(Point f1, Point f2, double r1, double r2);


//output of ellipse 
void print_ellipse(Ellipse e);

//perfomes direct comparison between two Ellipses 
int is_equal_ellipses(const Ellipse one, const Ellipse other);


//block of setters and getters
double ellipse_get_a(Ellipse e);
Ellipse ellipse_set_a(Ellipse e, double _a);

double ellipse_get_b(Ellipse e);
Ellipse ellipse_set_b(Ellipse e, double _b); 

double ellipse_get_x0(Ellipse e);
Ellipse ellipse_set_x0(Ellipse e, double _x0);

double ellipse_get_y0(Ellipse e);
Ellipse ellipse_set_y0(Ellipse e, double _y0);  

Point ellipse_get_f1(Ellipse e);

Point ellipse_get_f2(Ellipse e);


double ellipse_area(Ellipse e);

//getting center of ellipse
Point ellipse_get_center(Ellipse e);

//checks whether point lies inside given ellipse
bool is_inside_ellipse(Ellipse e, Point p);

//Intersection points of ellipse with given line, given by equation y = kx + h
void ellipse_intersect_line(Ellipse e, double k, double h);


//Intersection points of ellipses (found numerically, rounded off)
void intersect_ellipses(Ellipse e1, Ellipse e2);



typedef struct {
/*using unverted standard formula for parabola: (x-x0)^2 = 2p(y-y0)
Parabola axis is parallel to Oy-axis*/
    double v[3]; //vector, that contains p, x0, y0
} Parabola;


Parabola init_parabola( double p, double x0, double y0);
    /*creating parabola with given paremeters*/ 


Parabola parabola_from_array (double *array);
/*Parabola from given array of length 3 with given parameters*/

//block of setters and getters
double parabola_get_p(Parabola p);
Parabola parabola_set_p(Parabola p, double _p); 

double parabola_get_x0(Parabola p);
Parabola parabola_set_x0(Parabola p, double _x0);

double parabola_get_y0(Parabola p);
Parabola parabola_set_y0(Parabola p, double _y0);  


//output of parabola
void print_parabola(Parabola p);

Point parabola_get_center(Parabola p);

 
bool is_above_parabola(Parabola parab, Point p);
    /*checks whether point lies above given parabola*/
    

//Intersection points of Parabola with given line, given by equation y = kx + b
void parabola_intersect_line(Parabola p, double k, double b);

double parabola_intersect_line_area( Parabola p, double k, double b);
    

////Intersection points of Parabola with given Parabola
void parabola_intersect_parabola(Parabola p1, Parabola p2);


////Intersection area of Parabola with given Parabola
double parabola_intersect_parabola_area(Parabola p1, Parabola p2);

//////////////


typedef struct{
/*using inverted standard formula for hyperbola: (y-y0)^2/b^2 - (x-x0)^2/a^2 = 1
Hyperbola axis is parallel to Oy-axis*/
    double v[4]; //vector, that contains a, b, x0, y0
} Hyperbola;


Hyperbola init_hyperbola( double a, double b, double x0, double y0 );
    /*creating Hyperbola with given paremeters*/ 
    


Hyperbola hyperbola_from_array (double *array);
/*Hyperbola from given array of length 4 with given parameters*/
    

//block of setters and getters
double hyperbola_get_a(Hyperbola p);
Hyperbola hyperbola_set_a(Hyperbola p, double _a);

double hyperbola_get_b(Hyperbola p);
Hyperbola hyperbola_set_b(Hyperbola p, double _b);

double hyperbola_get_x0(Hyperbola p);
Hyperbola hyperbola_set_x0(Hyperbola p, double _x0);

double hyperbola_get_y0(Hyperbola p);
Hyperbola hyperbola_set_y0(Hyperbola p, double _y0);


bool is_equal_hyperbolas(Hyperbola one, Hyperbola other);
//performes direct comparison between hyperbolas

void print_hyperbola(Hyperbola h);
    /*output of hyperbola*/
    

Point hyperbola_get_center(Hyperbola h);
    /*getting center of hyperbola*/


 bool is_outside_hyperbola(Hyperbola h, Point p);
    /*checks whether point lies outside given hyperbola
    (above upper part, below lower part)*/
 


// ////Intersection points of Hyperbola with given Parabola
void hyperbola_intersect_parabola(Hyperbola h, Parabola p);
    
// ////Intersection points of Hyperbola with given Hyperbola
void hyperbola_intersect_hyperbola(Hyperbola one, Hyperbola other);
