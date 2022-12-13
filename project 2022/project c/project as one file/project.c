#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <setjmp.h>

#ifdef NAN
/* NAN is supported */
#endif
#ifdef INFINITY
/* INFINITY is supported */
#endif


// #include <iostream>
// #include <string>
// #include <fstream>
// #include <cstring> //memmove, memcpy
// #include <istream>
// #include <ostream>
// #include <sstream>
// #include <cmath>

#define _USE_MATH_DEFINES



typedef struct{
    //a2*x^2 + a1 * x + a0 = 0
    double a2;
    double a1;
    double a0;
}Quadratic_eq ;


Quadratic_eq init_quadratic(double _a2, double _a1, double _a0 ){
    Quadratic_eq eq;
    eq.a0 = _a0;
    eq.a1 = _a1;
    eq.a2 = _a2;
    return eq;
}



//output of solvings in stream 
double* solve(Quadratic_eq eq){
    static double res[2];
    if (eq.a2 == 0){
        if ((eq.a1 == 0) && (eq.a0 != 0)){
            res[0] = NAN;
            res[1] = NAN;
        } else if ((eq.a1 == 0) && (eq.a0 == 0)){
            res[0] = INFINITY;
            res[1] = INFINITY;
        } else {
            res[0] = -eq.a0/eq.a1;
            res[1] = res[0]; 
        }   
    }
    else if (eq.a1 == 0) {
        if (-eq.a0/eq.a2 < 0) {
            res[0] = NAN;
            res[1] = NAN;
        }
        else {
            res[0] = sqrt(-eq.a0/eq.a2);
            res[1] = -sqrt(-eq.a0/eq.a2);
        }
    }
    else if (eq.a0 == 0){
        res[0] = 0;
        if ((eq.a2 == 0) && (eq.a1 != 0)){
            res[1] = 0;
        } else if ((eq.a2 == 0) && (eq.a1 == 0)){
            res[0] = INFINITY;
            res[1] = INFINITY;
        } else res[1] = -eq.a1/eq.a2;
    }
    else{
        
        double D = (eq.a1*eq.a1 - 4*eq.a2*eq.a0);
        if (D < 0){
            res[0] = NAN;
            res[1] = NAN;
        } else if (D == 0){
            res[0] = -eq.a1/2/eq.a2;
            res[1] = -eq.a1/2/eq.a2;
        }
        else {
            res[0] = (-eq.a1 + sqrt(D))/2/eq.a2;
            res[1] = (-eq.a1 - sqrt(D))/2/eq.a2 ;
        }
    }
    return res;
} 



typedef struct{
/*Structure that creates 2d point*/
    double x, y;
} Point;

void print_point(Point p){
    printf("Point (%g,%g)\n", p.x, p.y);
}


//Задання еліпса у пдск, осі якого парлельні осям координат
typedef struct {
//canonical ellipse eqation (x-x0)^2/a^2 + (y-y0)^2/b^2 = 1
    double a, b, x0, y0;

}Ellipse;


//creating Ellipse with given paremeters задання еліпса координатами фокусів еліпса та 2 радіусами    
Ellipse init_ellipse(Point f1, Point f2, double r1, double r2) {
    Ellipse e;
    if (r1 > 0 && r2 > 0){
        if (f1.y == f2.y){
            double c = fabs(f1.x-f2.x)/2;
            e.a = (r1 + r2)/2; //using formula r1 = a + ex0; r2 = a-ex0 (M(x0,y0) - lies on ellipse)
            e.b = sqrt(e.a*e.a - c*c);
            e.x0 = fmax(f1.x, f2.x) - c;
            e.y0 = f1.y;
        }
        else if(f1.x == f2.x){
            double c = fabs(f1.y-f2.y)/2;
            e.b = (r1 + r2)/2; 
            e.a = sqrt(e.b*e.b - c*c);
            e.y0 = fmax(f1.y, f2.y) - c;  /// MAX ATTENTION
            e.x0 = f1.x;
        } else { 
            printf("Inappropriate Focuses input! (Not for standard equation), returning empty struct instance \n"); //ATTENTION
        }        
    } else {
        printf("Inappropriate Raduis input, returning empty struct Ellipse instance \n");
    }
    return e;
}


// 

// void input(Ellipse e){
//     double r1, r2;
//     Point f1, f2;
//     printf("Enter coordinates of F1: \n");                                           
//     scanf("%lf %lf", &f1.x, &f1.y);
//     printf("Enter coordinates of F2: \n");
//     scanf("%lf %lf", &f2.x, %f2.y);
//     printf("Enter 2 ellipse radiuses: \n");
//     scanf("%lf %lf", &r1, &r2);
//     init(f1, f2, r1, r2);
// }

//output of ellipse 
void print_ellipse(Ellipse e){
    printf("Ellipse (x+(%g))^2/(%g)^2 + ((y+(%g))^2/(%g)^2 = 1 \n", -e.x0, e.a, -e.y0, e.b ) ;
}

//perfomes direct comparison between two Ellipses 
int is_equal_ellipses(const Ellipse one, const Ellipse other){
    return one.a == other.a &&  one.b == other.b && one.x0 == other.x0 && one.y0 == other.y0;
}


//block of setters and getters
double ellipse_get_a(Ellipse e) {return e.a;}
Ellipse ellipse_set_a(Ellipse e, double _a) {e.a = _a; return e;}

double ellipse_get_b(Ellipse e) {return e.b;}
Ellipse ellipse_set_b(Ellipse e, double _b) {e.b = _b; return e;}

double ellipse_get_x0(Ellipse e) {return e.x0;}
Ellipse ellipse_set_x0(Ellipse e, double _x0) {e.x0 = _x0; return e;}

double ellipse_get_y0(Ellipse e) {return e.y0;}
Ellipse ellipse_set_y0(Ellipse e, double _y0) {e.y0 = _y0; return e;}  

Point ellipse_get_f1(Ellipse e){
    Point f1;
    double c;
    if (e.a >= e.b){
        c = sqrt(e.a*e.a - e.b*e.b);
        f1.x = e.x0-c;
        f1.y = e.y0;
    } else{
        c = sqrt(e.b*e.b - e.a*e.a);
        f1.x = e.x0;
        f1.y = e.y0-c;
    }
    return f1;
}

Point ellipse_get_f2(Ellipse e){
    Point f2;
    double c;
    if (e.a >= e.b){
        c = sqrt(e.a*e.a - e.b*e.b);
        f2.x = e.x0+c;  
        f2.y = e.y0;
    } else{
        c = sqrt(e.b*e.b - e.a*e.a);
        f2.x = e.x0;
        f2.y = e.y0+c;
    }
    return f2;
}


double ellipse_area(Ellipse e){
    return M_PI*e.a*e.b;                 //ATTENTION!
}

//getting center of ellipse
Point ellipse_get_center(Ellipse e){
    Point p = {e.x0, e.y0};
    return p;
}

//checks whether point lies inside given ellipse
bool is_inside_ellipse(Ellipse e, Point p){
    return pow(p.x +(-e.x0), 2)/(e.a*e.a) + pow(p.y+(-e.y0), 2)/(e.b*e.b) <= 1;
}

//Intersection points of ellipse with given line, given by equation y = kx + h
void ellipse_intersect_line(Ellipse e, double k, double h){
    /*x^2*(b^2 + a^2*k^2)+x(2a^2*k*h*y0 - 2b^2*x0) + (b^2*x0^2 + a^2*h^2*y0^2-a^2*b^2) = 0: result - intersection points*/
    
    double val;
    int num = 0;
    Quadratic_eq eq;
    eq.a2 = e.b*e.b + e.a*e.a*k*k;
    eq.a1 = 2*e.a*e.a*k*(h-e.y0) - 2*e.b*e.b*e.x0;
    eq.a0 = e.b*e.b*e.x0*e.x0 + e.a*e.a*(h-e.y0)*(h-e.y0) - e.a*e.a*e.b*e.b;

    Point p;

    double* sol = solve(eq);
     
    for (int i = 0; i < 2; i++) {
        if (sol[i]!=sol[i]){
            printf("NAN\n");
            break;
        } else {
            p.x = sol[i];
            p.y = k*p.x + h;
            print_point(p);
        }
    }
       
}


//Intersection points of ellipses (found numerically, rounded off)
void intersect_ellipses(Ellipse e1, Ellipse e2){

    Point p;
    double h = 0.001;
    double min_x = e1.x0 - e1.a;
    double max_x = e1.x0 + e1.a;
    int num = 0; 


    p.x = min_x - h;
    p.y = e1.y0;
    bool u_check_last = is_inside_ellipse(e2, p);  //setting check for upper part of the ellipse 
    bool l_check_last = u_check_last; //setting check for lower part of the ellipse

    bool u_check = u_check_last;
    bool l_check = l_check_last;


    //y = +-b/a*sqrt(a*a-(x-x0)*(x-x0)) + y0;
    if (is_equal_ellipses(e1, e2)){
        printf("INFINITY\n");

    } else {  
        for (float x = min_x; x <= max_x; x += h){
            p.x = x;
            p.y = e1.b/e1.a*sqrt(e1.a*e1.a-pow(x - e1.x0, 2)) + e1.y0; 
            u_check = is_inside_ellipse(e2, p);

            if (u_check != u_check_last){  
                print_point(p);
                num += 1;
            }
            u_check_last = u_check;
        
            p.y = -e1.b/e1.a*sqrt(e1.a*e1.a-pow(x - e1.x0, 2)) + e1.y0; 
            l_check = is_inside_ellipse(e2, p);

            if (l_check != l_check_last){  
                print_point(p);
                num += 1;
            }
            l_check_last = l_check;

        }

        if (!num){
            printf("NAN\n");
        }
    }
}



typedef struct {
/*using unverted standard formula for parabola: (x-x0)^2 = 2p(y-y0)
Parabola axis is parallel to Oy-axis*/
    double v[3]; //vector, that contains p, x0, y0
} Parabola;


Parabola init_parabola( double p, double x0, double y0) {
    /*creating parabola with given paremeters*/ 
    Parabola parab;
    parab.v[0] = p;
    parab.v[1] = x0;
    parab.v[2] = y0;
    return parab;
}


Parabola parabola_from_array (double *array){
/*Parabola from given array of length 3 with given parameters*/
    Parabola parab;
    memmove(parab.v, array, sizeof(*parab.v)*3);
    // for (int i = 0; i < 3; i++){
    //     parabola->equation[i] = array[i];
    // }
    return parab;
}



// void input(){
//     cout << "Enter p, x0, y0 parameters of parabola (x-x0)^2 = 2p(y-y0): "<< endl;
//     for(int i = 0; i < 3; i++){
//         cin >> v[i];
//     } 
// }

//block of setters and getters
double parabola_get_p(Parabola p) {return p.v[0];}
Parabola parabola_set_p(Parabola p, double _p) {p.v[0] = _p; return p;}

double parabola_get_x0(Parabola p) {return p.v[1];}
Parabola parabola_set_x0(Parabola p, double _x0) {p.v[1] = _x0; return p;}

double parabola_get_y0(Parabola p) {return p.v[2];}
Parabola parabola_set_y0(Parabola p, double _y0) {p.v[2] = _y0; return p;}  


//output of parabola
void print_parabola(Parabola p) {
    printf("Parabola (x+(%g))^2 = %g*2*(y+(%g)) \n", -p.v[1], p.v[0], -p.v[2]);
}


Point parabola_get_center(Parabola p){
    /*getting center of parabola (focus)*/
    Point center = {.x=p.v[1], .y=p.v[2]+p.v[0]/2};
    return center;
}

 
int is_above_parabola(Parabola parab, Point p){        //ATTENTION!
    /*checks whether point lies above given parabola*/
    double a = pow(p.x-parab.v[1], 2);
    double b = 2*parab.v[0]*(p.y-parab.v[2]);

    if (a < b) {

    }
    // printf("%g", bool(a<b));
    return a < b; //pow(p.x-parab.v[1], 2) < 2*parab.v[0]*(p.y-parab.v[2]);
}
    

//Intersection points of Parabola with given line, given by equation y = kx + b
void parabola_intersect_line(Parabola p, double k, double b){
    /*x^2 + x(-2x0 - 2pk) + x0*x0 + y0 - 2*b*p */
    double val;
    int num = 0;
    Quadratic_eq eq = {.a2=1, .a1=-2*p.v[1] - 2*p.v[0]*k, .a0= p.v[1]*p.v[1] + p.v[2] - 2*b*p.v[0]};
    Point point;

    double* sol = solve(eq);
     
    for (int i = 0; i < 2; i++) {
        if (sol[i]!=sol[i]){
            printf("NAN\n");
            break;
        } else {
            point.x = sol[i];
            point.y = k*point.x + b;
            print_point(point);
        }
    }   
}

double parabola_intersect_line_area( Parabola p, double k, double b){
    /*x^2 + x(-2x0 - 2pk) + x0*x0 + y0 - 2*b*p */

       /*x^2 + x(-2x0 - 2pk) + x0*x0 + y0 - 2*b*p */
    double val;
    Quadratic_eq eq = {.a2=1,.a1=-2*p.v[1] - 2*p.v[0]*k, .a0=  p.v[1]*p.v[1] + p.v[2] - 2*b*p.v[0]};
    Point point1, point2;

    double* sol = solve(eq);
     
    for (int i = 0; i < 2; i++) {
        if (sol[i]!=sol[i]){
            printf("NAN\n");
            return 0;
        }
    } 

    point1.x = sol[0];
    point1.y = k*point1.x + b;
    point2.x = sol[1];
    point2.y = k*point2.x + b;
        
    return fabs((pow(point1.x, 3)-pow(point2.x, 3))/(6*p.v[0])  +
     (point1.x*point1.x-point2.x*point2.x)/2*(-p.v[1]/p.v[0] - k) +
      (p.v[1]*p.v[1]/2/p.v[0] + p.v[2] - b)*(point1.x - point2.x));   
}
    

////Intersection points of Parabola with given Parabola
void parabola_intersect_parabola(Parabola p1, Parabola p2){
    
    Quadratic_eq eq = {.a2 = p2.v[0] - p1.v[0], .a1=  2*(p2.v[1]*p1.v[0] - p1.v[1]*p2.v[0]), .a0= p2.v[0]*p1.v[1]*p1.v[1] - p1.v[0]*p2.v[1]*p2.v[1] + 2*p2.v[0]*p1.v[0]*(p1.v[2] - p2.v[2])};
    Point point;

    double* sol = solve(eq);
     
    for (int i = 0; i < 2; i++) {
        if (sol[i]!=sol[i]){
            printf("NAN\n");
            break;
        } else if (isinf(sol[i])){
            printf("Equal parabolas \n");
            break;
        } else {
            point.x = sol[i];
            point.y = pow((point.x - p1.v[1]), 2)/2/p1.v[0] + p1.v[2];
            print_point(point);
        }
    } 
}

double parabola_intersect_parabola_area(Parabola p1, Parabola p2){
    
    Quadratic_eq eq = { .a2=p2.v[0] - p1.v[0], .a1= 2*(p2.v[1]*p1.v[0] - p1.v[1]*p2.v[0]), .a0=p2.v[0]*p1.v[1]*p1.v[1] - p1.v[0]*p2.v[1]*p2.v[1] + 2*p2.v[0]*p1.v[0]*(p1.v[2] - p2.v[2])};
    Point point;

    double* sol = solve(eq);
     
    for (int i = 0; i < 2; i++) {
        if (sol[i]!=sol[i]){
            return 0;
        } else if (isinf(sol[i])){
            return 0;
        }
    } 
    double x1 = sol[0];
    double x2 = sol[1];
    
    return fabs((pow(x1, 3) - pow(x2, 3))*(p2.v[0] - p1.v[0])/3 + (x1*x1-x2*x2)*(p2.v[1]*p1.v[0] - p1.v[1]*p2.v[0]) +  
    (p2.v[0]*p1.v[1]*p1.v[1] - p1.v[0]*p2.v[1]*p2.v[1] + 2*p2.v[0]*p1.v[0]*(p1.v[2] - p2.v[2]))*(x1 - x2)) ;
}


//////////////


typedef struct{
/*using inverted standard formula for hyperbola: (y-y0)^2/b^2 - (x-x0)^2/a^2 = 1
Hyperbola axis is parallel to Oy-axis*/
    double v[4]; //vector, that contains a, b, x0, y0
} Hyperbola;


Hyperbola init_hyperbola( double a, double b, double x0, double y0 ) {
    /*creating Hyperbola with given paremeters*/ 
    
    Hyperbola hyperbola;
    if (a<=0 || b<= 0){
        printf("Inappropriate a, b parameters input, returning empty struct instance \n");
        return hyperbola; //ATTENTION
    }
    hyperbola.v[0] = a;
    hyperbola.v[1] = b;
    hyperbola.v[2] = x0;
    hyperbola.v[3] = y0;
    return hyperbola;
}



Hyperbola hyperbola_from_array (double *array){
/*Hyperbola from given array of length 4 with given parameters*/
    Hyperbola hyperbola;
    if (array[0]<=0 || array[0]<= 0) {
        printf("Inappropriate a, b parameters input, returning empty struct instance \n");  //ATTENTION
        return hyperbola;
    }
    memmove(hyperbola.v, array, sizeof(*hyperbola.v)*4);
    // for (int i = 0; i < 3; i++){
    //     parabola->equation[i] = array[i];
    // }
    return hyperbola;
}


//block of setters and getters
double hyperbola_get_a(Hyperbola p) {return p.v[0];}
Hyperbola hyperbola_set_a(Hyperbola p, double _a) {p.v[0] = _a; return p;}

double hyperbola_get_b(Hyperbola p) {return p.v[1];}
Hyperbola hyperbola_set_b(Hyperbola p, double _b) {p.v[1] = _b; return p;}

double hyperbola_get_x0(Hyperbola p) {return p.v[2];}
Hyperbola hyperbola_set_x0(Hyperbola p, double _x0) {p.v[2] = _x0; return p;}

double hyperbola_get_y0(Hyperbola p) {return p.v[3];}
Hyperbola hyperbola_set_y0(Hyperbola p, double _y0) {p.v[3] = _y0; return p;}  



// void input(){
//     cout << "Enter a, b, x0, y0 parameters of hyperbola (y-y0)^2/b^2 - (x-x0)^2/a^2 = 1: "<< endl;
//     for(int i = 0; i < 4; i++){
//         cin >> v[i];
//     } 
// }


int is_equal_hyperbolas(Hyperbola one, Hyperbola other) {return one.v[0] == other.v[0] &&  one.v[1] == other.v[1]
    && one.v[2] == other.v[2] && one.v[3] == other.v[3];}

void print_hyperbola(Hyperbola h){
    /*output of hyperbola*/
    printf("Hyperbola (y-(%g))^2/(%g)^2 - (x-(%g))^2/(%g)^2  = 1 \n", h.v[3], h.v[1], h.v[2], h.v[0]);
}


Point hyperbola_get_center(Hyperbola h){
    /*getting center of hyperbola*/
    Point p ={.x = h.v[2], .y = h.v[3]};
    return p;
}


 bool is_outside_hyperbola(Hyperbola h, Point p){   //ATTENTION!
    /*checks whether point lies outside given hyperbola
    (above upper part, below lower part)*/
    return (pow(p.y - h.v[3], 2) / pow(h.v[1], 2) - pow(p.x - h.v[2], 2) / pow(h.v[0], 2) >= 1);
}


// ////Intersection points of Hyperbola with given Parabola
    //Intersection points of ellipses (found numerically, rounded off)

void hyperbola_intersect_parabola(Hyperbola h, Parabola p){
    
    Point point;
    double step = 0.001;
    double min_x = -100;// not universal method
    double max_x = 100;//
    int num = 0; 


    point.x = min_x - step;
    point.y = pow(point.x - parabola_get_x0(p), 2)/2/parabola_get_p(p) + parabola_get_y0(p);
    bool check_last = is_outside_hyperbola(h, point);
    bool check = check_last; 

    
    for (float x = min_x; x <= max_x; x += step){
        point.x = x;
        point.y = pow(point.x - parabola_get_x0(p), 2)/2/parabola_get_p(p) + parabola_get_y0(p); 
        check = is_outside_hyperbola(h, point);

        if (check != check_last){  
            print_point(point);
            num += 1;
        }
        check_last = check;
    }

    if (!num){
        printf("NAN \n");
    }
    
}
    

void hyperbola_intersect_hyperbola(Hyperbola one, Hyperbola other){
    
    Point p;
    double step = 0.001;
    double max_y = 100;
    double min_y = -100;
    int num = 0; 


    p.y = min_y;  // other.v[3]+other.v[1]  // y0 + b
    p.x = other.v[0]/other.v[1]*sqrt(pow(p.y - other.v[3], 2) - one.v[1]*one.v[1]) + other.v[2];
    bool u_check_last = is_outside_hyperbola(one, p); 
    bool l_check_last = u_check_last; //setting check 
    
    bool u_check = l_check_last;
    bool l_check = u_check;


    if (is_equal_hyperbolas(one, other)){
        printf("INFINITY \n");
    } else if ((one.v[1] + one.v[3]> other.v[1] + other.v[3]) && (-one.v[1] + one.v[3] < -other.v[1] + other.v[3]) &&
                one.v[0]/one.v[1] <= other.v[0]/other.v[1] ) {
        printf("NAN\n");
    }  else if (one.v[2] == other.v[2] && one.v[3] == other.v[3] && one.v[1] == other.v[1]) {
        p.x = one.v[2]; 
        p.y =  one.v[1] + one.v[3];
        print_point(p);
        p.y = -one.v[1] + one.v[3];
        print_point(p);
    } else {  
        for (float y = min_y; y <= max_y; y += step){
            if ((y < other.v[3]+other.v[1]) && (y > - other.v[3]+other.v[1])){ continue; }

            p.y = y;
            p.x = other.v[0]/other.v[1]*sqrt(pow(y - other.v[3], 2) - other.v[1]*other.v[1]) + other.v[2]; 
            u_check = is_outside_hyperbola(one, p);

            if (u_check != u_check_last){  
                print_point(p);
                num += 1;
            }
            u_check_last = u_check;
            
            p.x = -other.v[0]/other.v[1]*sqrt(pow(y - other.v[3], 2) - other.v[1]*other.v[1]) + other.v[2]; 

            l_check = is_outside_hyperbola(one, p);
            if (l_check != l_check_last){  
                print_point(p);
                num += 1;
            }
            l_check_last = l_check;
        }

        if (!num){
            printf("NAN\n");
        }
    }
}



int main(int argc, char* argv[]) {
    // Point p1{1, 1}; //
    // Point p2{-5, 1};
    // Ellipse e;
    // e = init_ellipse(p1,p2, 5, 5);
    // print_ellipse(e);
    // //printf("Ellipse area: %g\n", ellipse_area(e));
    // // Quadratic_eq eq{1, -3 , 2};  
    // // double* res = solve(eq);
    // // for (int i = 0; i < 2; i++){
    // //     printf("%g \n", res[i]);
    // // }

    // // //intersection of ellipse with line y = 5; (one point)
    // // ellipse_intersect_line(e, 0, 5); 
    // // //intersection of ellipse with line y = 3x+1; (2 points)
    // // ellipse_intersect_line(e, 3, 1);
    // // //intersection of ellipse with line y = 6; (no points)
    // // ellipse_intersect_line(e, 0, 6); 

    // //intersection of two ellipses
    // Ellipse e2 = e;
    // printf("Equal ellipses intersection: \n");
    // intersect_ellipses(e, e2);
    // e2 = ellipse_set_a(e2, 4);
    // e2 = ellipse_set_b(e2, 5);
    // print_ellipse(e2);
    // printf("Points of ellipses intersection: \n");
    // intersect_ellipses(e, e2);

    //// PARABOLA BLOCK

    // Point p3{0, 3}; 
    // Parabola par {0.5, 0 , 0}; //y = x^2
    // print_parabola(par);
    // printf("Center of Parabola: %g \n", parabola_get_center(par) );
    // printf("is point (0,3) above_parabola: %d \n", is_above_parabola(par, p3));     // wtf?
    // printf("Points of intersection of given parabola and line (y = 0*x + 4)\n");
    // parabola_intersect_line(par, 0, 4);
    // printf("Intersection area of given parabola and line (y = 0*x + 4): %g\n\n", parabola_intersect_line_area(par, 0, 4));
    // // 
 
    // printf("Intersection of two identical parabolas: \n");
    // parabola_intersect_parabola(par, par);
    // printf("Area of intersection: %g\n\n", parabola_intersect_parabola_area(par, par));


    // Parabola par2{-1, 1, 4};  // 
    // print_parabola(par2);
    // printf("Intersection of different parabolas (2 points of intersection): \n");
    // parabola_intersect_parabola(par, par2);
    // printf("Area of intersection: %g \n\n", parabola_intersect_parabola_area(par, par2));  // 7.64 (by hand calculations)
    

    // Parabola par3{1, 0, 0};  // 
    // print_parabola(par3);
    // printf("Intersection of different parabolas (1 points of intersection): \n");
    // parabola_intersect_parabola(par, par3);
    // printf("Area of intersection: %g \n\n", parabola_intersect_parabola_area(par, par3));
    

    // Parabola par4{0.5, 0, 4};  // 
    // print_parabola(par4);
    // printf("Intersection of different parabolas (0 points of intersection): \n");
    // parabola_intersect_parabola(par, par4);
    // printf("Area of intersection: %g\n\n", parabola_intersect_parabola_area(par, par4));


    //Hyperbola block

    Point p4 = {.x = 0, .y= -4};
    Point p5 = {0, -5};
    Hyperbola  hyp1 = init_hyperbola(3, 4, 0, 0);
    print_hyperbola(hyp1) ;
    printf("Center of Hyperbola: ");
    print_point(hyperbola_get_center(hyp1));
    printf("is point (0,-4) outside hyperbola: %d \n", is_outside_hyperbola(hyp1, p4));
    printf("is point (0,-5) outside hyperbola: %d \n\n", is_outside_hyperbola(hyp1, p5));

    printf("Intersection of Hyperbola with parabola: \n");
    Parabola p = {0.5, 0 , 0}; //y = x^2
    print_parabola(p);
    hyperbola_intersect_parabola(hyp1, p);



    printf("\nIntersection of two identical Hyperbolas: \n");
    hyperbola_intersect_hyperbola(hyp1, hyp1);
    

    Hyperbola hyp2 = {1, 4, 0, 0};   
    print_hyperbola(hyp2);
    printf("Intersection of different Hyperbolas (2 points of intersection): \n");
    hyperbola_intersect_hyperbola(hyp1, hyp2);
    
    printf("\n");

    Hyperbola hyp3 = {1, 3, 1, 1};   
    print_hyperbola(hyp3);
    printf("Intersection of different Hyperbolas (4 points of intersection): \n");
    hyperbola_intersect_hyperbola(hyp1, hyp3);
    
    printf("\n");

    Hyperbola hyp4 ={4, 3, 0, 0};   
    print_hyperbola(hyp4);
    printf("Intersection of different Hyperbolas (0 points of intersection): \n");
    hyperbola_intersect_hyperbola(hyp1, hyp4);
    
    printf("\n");

    // 4, 3, 0, 0 - NAN ; 1, 4, 0, 0 - 2 points; 1 3 0 1 - 4 points; 3 4 0 0 - inf (equal)

    return 0;
}
