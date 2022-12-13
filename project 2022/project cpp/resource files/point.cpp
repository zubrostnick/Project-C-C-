#include "point.h"


//POINT DESCRIBTION 

//default
Point::Point() {}

// main constructor 
Point::Point(double _x, double _y): x(_x), y(_y) {}

void Point::input(istream& inp  = cin) {
    /*input of parameters x,y from stream (cin - by default)*/
    inp >> x >> y;
}

void Point::print(ostream& out = cout) {
    out << "Point(" << x << "," << y <<")"<<ends;
}

// TODO: ...
ostream& operator<<(ostream& out, Point& p){  
    //must be non-const, but somehow it doesn't funciton as it should 
    //(ellipse, hyperb., parab. realisations are similar); - works as a united file (project.cpp)
    // updated: worked in vs
    // it worked with const { g++ source.cpp -o source -fpermissive } -vscode 
    //
    p.print(out);
    //out << "Point(" << p.x << "," << p.y << ")" << ends;
    return out;
}

istream& operator>>(istream& inp, Point& p) {
    p.input(inp);
    return inp;
}

//creates and returns a point 
Point Point::create_input(){
    Point p;
    cout << "Input of point, please write in x, y coordinates of a point:" << endl;
    p.input();
    return p;
}
