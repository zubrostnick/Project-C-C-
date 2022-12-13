
#include "point.h"
#include "conic_error.h"
#include "linear.h"
#include "quadratic.h"
#include "biquadratic.h"
#include "from_file.h"
#include "parabola.h"
#include "ellipse.h"
#include "hyperbola.h"
#include <iostream>

void test_ellipse(){
    Point p1 = Point(1, 1); 
    Point p2 = Point(-5, 1);
    cout << "Creating ellipse with r1=r2=5 and f1(1,1), (-5,1)" << endl;
    Ellipse e(p1,p2, 5, 5);
    cout << e << endl;
    cout << "Ellipse area:" << e.area() << endl;
    
    cout << "Intersection of ellipse with line y = 5; (one point):" << endl;
    e.intersect_line(0, 4, cout); 
    cout << "Intersection of ellipse with line y = 3x+1; (2 points)" << endl;
    e.intersect_line(3, 1, cout);
    cout << "Intersection of ellipse with line y = 6; (no points)" << endl;
    e.intersect_line(0, 6, cout); 

    //intersection of two ellipses
    cout << "Testing intersection of two ellipses: "<<endl;
    Ellipse e2(e);
    e2.set_a(4);
    e2.set_b(5);
    cout << "setters and getters were tested"<<endl;
    cout << e2 << endl;
    
    // Ellipse e3 = Ellipse::create_input();  // this one obliges entering values from keyboard
    cout << "Points of ellipses intersection: "<< endl;
    e.intersect_ellipse(e2, cout);
    cout << "Ellipse testing is over" << endl;
}

void test_parabola(){
    //Parabola p = Parabola::create_input(); // this one obliges entering values from keyboard
    Point p3 = Point(0, 3); // 
    cout << "Creating parabola y = x^2:" << endl;
    Parabola par(0.5, 0 , 0); //y = x^2
    cout << par << endl; 
    cout << "Center of Parabola:" << endl;
    par.get_center().print(cout);
    cout << "is point (0,3) above_parabola: " << par.is_above_parabola(p3)<< endl;
    cout << "Points of intersection of given parabola and line (y = 0*x + 4): " <<  endl;
    par.intersect_line(0, 4, cout);
    cout << "Intersection area of given parabola and line (y = 0*x + 4): " << par.intersect_line_area(0, 4) << endl;
    cout << endl;

    cout << "Intersection of two identical parabolas: " << endl;
    par.intersect_parabola(par, cout);
    cout << "Area of intersection: " << par.intersect_parabola_area(par) << endl;  
    cout << endl;

    Parabola par2(-1, 1, 4);  // 
    cout << par2 << endl;
    cout << "Intersection of different parabolas (2 points of intersection): " << endl;
    par.intersect_parabola(par2, cout);
    cout << "Area of intersection: " << par.intersect_parabola_area(par2) << endl;  // 7.64 (by hand calculations)
    cout << endl;

    Parabola par3(1, 0, 0);  // 
    cout << par3 << endl;
    cout << "Intersection of different parabolas (1 points of intersection): " << endl;
    par.intersect_parabola(par3, cout);  
    cout << "Area of intersection: "<< par.intersect_parabola_area(par3) << endl;
    cout << endl;

    Parabola par4(0.5, 0, 4);  // 
    cout << par4 << endl;
    cout << "Intersection of different parabolas (0 points of intersection): " << endl;
    par.intersect_parabola(par4, cout);
    cout << "Area of intersection: " << par.intersect_parabola_area(par4) << endl;

    cout << "Parabola testing is over" << endl;
}


void test_hyperbola(){
    Point p4 = Point(0, -4);
    Point p5 = Point(0, -5);
    //Hyperbola h = Hyperbola::create_input();  // 3, 4, 0, 0 // obliges by hand input
    cout << "Creating hyperbola with a = 3, b = 4, (x,y0) = (0,0)" << endl;
    Hyperbola  h(3, 4, 0, 0);
    cout << h << endl;
    cout << "Center of Hyperbola: " << endl;
    h.get_center().print(cout);
    cout << "is point (0,-4) outside hyperbola: " << h.is_outside_hyperbola(p4) << endl;
    cout << "is point (0,-5) outside hyperbola: " << h.is_outside_hyperbola(p5) << endl;
    cout << endl;

    Parabola p(0.5, 0 , 0); //y = x^2
    cout << "Intersection of Hyperbola with \n" << p <<  ":" << endl;
    h.intersect_parabola(p, cout);
    cout << endl;

    Hyperbola  h2(1, 3, 0, 1);
    cout << "Intersection of Hyperbola with \n" << h2 <<  "(4 points of intersection) :" << endl;
    h.intersect_hyperbola(h2, cout);
    cout << endl;
    
    Hyperbola  h3(4, 3, 0, 0);// 4, 3, 0, 0 - NAN ; 1, 4, 0, 0 - 2 points; 1 3 0 1 - 4 points; 3 4 0 0 - inf (equal)
    cout << "Intersection of Hyperbola with \n" << h3 <<  "(0 points of intersection) :" << endl;
    h.intersect_hyperbola(h3, cout);
    cout << endl;

    Hyperbola  h4(1, 3, 0, 1);
    cout << "Intersection of Hyperbola with \n" << h4 <<  "(2 points of intersection) :" << endl;
    h.intersect_hyperbola(h4, cout);
    cout << endl;

    cout << "Intersection of Hyperbola with \n" << h4 <<  "(inf points of intersection) :" << endl;
    h.intersect_hyperbola(h4, cout);

    cout << "\nHyperbola testing is over" << endl;
}


int main(int argc, char* argv[]) {
    
    int mode = 0;
    while (true){
    
        cout << "Enter testing mode: 1 - Ellipse testing, 2 - Parbola testing, 3 - Hyeprbola testing, 4 - exit" << endl;
        cin >> mode;
        if (mode == 1){
            test_ellipse();
        } else if ( mode == 2 ) {
            test_parabola();
        } else if ( mode == 3 ){
            test_hyperbola();
        } else if ( mode == 4 ) break;
    }

    return 0;
}