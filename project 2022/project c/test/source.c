#include "curves.h"
#include "stdio.h"

void test_ellipse(){
    Point p1 = {1, 1}; //
    Point p2 = {-5, 1};
    Ellipse e;
    printf("Creating ellipse with r1=r2=5 and f1(1,1), (-5,1):\n");
    e = init_ellipse(p1,p2, 5, 5);
    print_ellipse(e);
    printf("Ellipse area: %g\n\n", ellipse_area(e));

    printf("Intersection of ellipse with line y = 5; (one point)\n");
    ellipse_intersect_line(e, 0, 5); 
    printf("Intersection of ellipse with line y = 3x+1; (2 points\n");
    ellipse_intersect_line(e, 3, 1);
    printf("Intersection of ellipse with line y = 6; (no points)\n");
    ellipse_intersect_line(e, 0, 6); 

    printf("\nTesting intersection of two ellipses:\n");
    Ellipse e2 = e;
    printf("Equal ellipses intersection: \n");
    intersect_ellipses(e, e2);
    printf("Testing setting new parameters for ellipse :\n");
    e2 = ellipse_set_a(e2, 4);
    e2 = ellipse_set_b(e2, 5);
    printf("setters and getters were tested\n");
    print_ellipse(e2);
    printf("Points of two different ellipses intersection: \n");
    intersect_ellipses(e, e2);
    printf("Ellipse testing is over\n");
}

void test_parabola(){
    Point p3 = {0, 3}; 
    printf("Creating parabola y = x^2:\n");
    Parabola par = {0.5, 0 , 0}; 
    print_parabola(par);
    printf("Center of Parabola: %g \n", parabola_get_center(par) );
    printf("is point (0,3) above_parabola: %d \n", is_above_parabola(par, p3));     
    printf("Points of intersection of given parabola and line (y = 0*x + 4)\n");
    parabola_intersect_line(par, 0, 4);
    printf("Intersection area of given parabola and line (y = 0*x + 4): %g\n\n", parabola_intersect_line_area(par, 0, 4));
    // 
 
    printf("Intersection of two identical parabolas: \n");
    parabola_intersect_parabola(par, par);
    printf("Area of intersection: %g\n\n", parabola_intersect_parabola_area(par, par));


    Parabola par2 = {-1, 1, 4};  // 
    print_parabola(par2);
    printf("Intersection of different parabolas (2 points of intersection): \n");
    parabola_intersect_parabola(par, par2);
    printf("Area of intersection: %g \n\n", parabola_intersect_parabola_area(par, par2));  // 7.64 (by hand calculations)
    

    Parabola par3 = {1, 0, 0};  // 
    print_parabola(par3);
    printf("Intersection of different parabolas (1 points of intersection): \n");
    parabola_intersect_parabola(par, par3);
    printf("Area of intersection: %g \n\n", parabola_intersect_parabola_area(par, par3));
    

    Parabola par4 = {0.5, 0, 4};  // 
    print_parabola(par4);
    printf("Intersection of different parabolas (0 points of intersection): \n");
    parabola_intersect_parabola(par, par4);
    printf("Area of intersection: %g\n\n", parabola_intersect_parabola_area(par, par4));

    printf("Parabola testing is over\n");
}

void test_hyperbola(){
    //Hyperbola block

    Point p4 = {.x = 0, .y= -4};
    Point p5 = {0, -5};
    printf("Creating hyperbola with a = 3, b = 4, (x,y0) = (0,0)\n");
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
    printf("Hyperbola testing over\n");

}

int main(int argc, char* argv[]) {
    int mode = 0;
    while (1){
    
        printf("Enter testing mode: 1 - Ellipse testing, 2 - Parbola testing, 3 - Hyeprbola testing, 4 - exit\n");
        scanf_s("%i", &mode);
        if (mode == 1){
            test_ellipse();
        } else if ( mode == 2) {
            test_parabola();
        } else if ( mode == 3 ){
            test_hyperbola();
        } else if ( mode == 4) break;
    }

    return 0;
}
