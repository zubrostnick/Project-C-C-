#include "parabola.h"


/*using unverted standard formula for parabola: (x-x0)^2 = 2p(y-y0)
Parabola axis is parallel to Oy-axis*/
Parabola::Parabola( double p, double x0, double y0) {
    /*creating parabola with given paremeters*/ 
    
    v = new double[3];
    v[0] = p;
    v[1] = x0;
    v[2] = y0;
}

//default
Parabola::Parabola() {v = new double[3];}

Parabola::Parabola(double *ptr_v){
/*Parabola from given array of length 3 with given parameters*/
    
    v = new double[3];
    memmove(v, ptr_v, sizeof(*v)*3);
}


Parabola::~Parabola(){
//destructor
    delete[] v;
};


Parabola::Parabola(const Parabola& a){
/*copy-constructor: create Parabola from given*/
    v = new double[3];
    memmove(v, a.v, sizeof(*v)*3);
}

//input n-th coordinate of vector
//create and input vector of given size;
Parabola Parabola::create_input(){
    Parabola a;
    a.input();
    return a;
}

void Parabola::input(){
    cout << "Enter p, x0, y0 parameters of parabola (x-x0)^2 = 2p(y-y0): "<< endl;
    for(int i = 0; i < 3; i++){
        cin >> v[i];
    } 
}

//block of setters and getters
double Parabola::get_p() {return v[0];}
void Parabola::set_p(double _p) {this->v[0] = _p;}

double Parabola::get_x0() {return v[1];}
void Parabola::set_x0(double _x0) {this->v[1] = _x0;}

double Parabola::get_y0() {return v[2];}
void Parabola::set_y0(double _y0) {this->v[2] = _y0;}  


//output of parabola
void Parabola::print(ostream& out = cout) {
    out << "Parabola (x+("<< -v[1] <<"))^2 = " << v[0] << "*2*(y+("<< -v[2] <<"))"<<ends;
}


ostream& operator<<(ostream& out, Parabola& p) {
    p.print(out);
    return out;
}

Point Parabola::get_center(){
    /*getting center of parabola (focus)*/
    Point p(v[1], v[2]+v[0]/2);
    return p;
}


bool Parabola::is_above_parabola(Point p){
    /*checks whether point lies above given parabola*/
    return (pow(p.x-v[1], 2) < 2*v[0]*(p.y-v[2]));
}


//Intersection points of Parabola with given line, given by equation y = kx + b
void Parabola::intersect_line( double k, double b, ostream& out = cout){
    /*x^2 + x(-2x0 - 2pk) + x0*x0 + y0 - 2*b*p */
    
    ostringstream sout;
    double val;
    int num = 0;
    Quadratic_eq e(1,-2*v[1] - 2*v[0]*k,  v[1]*v[1] + v[2] - 2*b*v[0]);
    //cout << e << endl;
    Point p;

    e.solve(sout);
    istringstream i_res(sout.str());

    while ((i_res >> val)) {
        num++;
        p.x = val;
        p.y = k*p.x + b;
        out << p << endl;   
    }
        if (sout.str().find("nan") != string::npos) {  
        out <<sout.str() << ends;
    } else {
        if (num == 0){  //(!num)
            out << NAN << endl;
        }  
    }   
}

double Parabola::intersect_line_area( double k, double b){
    /*x^2 + x(-2x0 - 2pk) + x0*x0 + y0 - 2*b*p */
    
    ostringstream sout;
    intersect_line(k, b, sout);
    double x1, x2;
    int num = 0;
    Point p;
    string str;
    
    istringstream i_res(sout.str());
    
    while (i_res >> str){
        num ++;
        if (num == 1){
            x1 = stod( str.substr(str.find("(")+1, str.find(",")-str.find("(")-1));
            x2 = x1;
        } else if(num == 2){
            x2 = stod( str.substr(str.find("(")+1, str.find(",")-str.find("(")-1));
        }
    }   
    
    if (sout.str().find("nan") != string::npos) {  
        return 0;
    } else {
        if (num == 0){ 
            return 0;
        }
        else {
            return abs((pow(x1, 3)-pow(x2, 3))/(6*v[0])  + (x1*x1-x2*x2)/2*(-v[1]/v[0] - k) + (v[1]*v[1]/2/v[0] + v[2] - b)*(x1 - x2));
        }  
    }   
}

////Intersection points of Parabola with given Parabola
void Parabola::intersect_parabola(Parabola p, ostream& out = cout){
    
    ostringstream sout;
    double val;
    int num = 0;
    Quadratic_eq e(p.v[0] - v[0],  2*(p.v[1]*v[0] - v[1]*p.v[0]),  
    p.v[0]*v[1]*v[1] - v[0]*p.v[1]*p.v[1] + 2*p.v[0]*v[0]*(v[2] - p.v[2]));
    //cout << e << endl;
    Point point;

    e.solve(sout);
    istringstream i_res(sout.str());

    while ((i_res >> val)) {
        num++;
        point.x = val;
        point.y = pow((point.x - v[1]), 2)/2/v[0] + v[2];
        out << point << endl;   
    }
        if (sout.str().find("nan") != string::npos) {  
        out <<sout.str() << ends;
    } else {
        if (num == 0){  //(!num)
            out << NAN << endl;
        }  
    }   
}

double Parabola::intersect_parabola_area( Parabola p){
    
    ostringstream sout;
    intersect_parabola(p, sout);
    double x1, x2;
    int num = 0;
    Point point;
    string str;

    if (sout.str().find("nan") != string::npos) {  
        return 0;
    }
    
    istringstream i_res(sout.str());

    while (i_res >> str){
        num ++;
        if (num == 1){
            x1 = stod( str.substr(str.find("(")+1, str.find(",")-str.find("(")-1));
            x2 = x1;
        } else if(num == 2){
            x2 = stod( str.substr(str.find("(")+1, str.find(",")-str.find("(")-1));
        }
    }   

    
    if (num == 0){ 
        return 0;
    }
    else {
        return abs((pow(x1, 3) - pow(x2, 3))*(p.v[0] - v[0])/3 + (x1*x1-x2*x2)*(p.v[1]*v[0] - v[1]*p.v[0]) +  
    (p.v[0]*v[1]*v[1] - v[0]*p.v[1]*p.v[1] + 2*p.v[0]*v[0]*(v[2] - p.v[2]))*(x1 - x2)) ;
    }  
   
}
