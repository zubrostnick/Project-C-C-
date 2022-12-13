#include "hyperbola.h"

/*using inverted standard formula for hyperbola: (y-y0)^2/b^2 - (x-x0)^2/a^2 = 1
Hyperbola axis is parallel to Oy-axis*/

Hyperbola::Hyperbola(double a, double b, double x0, double y0) {
    /*creating parabola with given paremeters*/ 
    if (a<=0 || b<= 0)  throw ConicError("Inappropriate a, b parameters input");
    
    v = new double[4];
    v[0] = a;
    v[1] = b;
    v[2] = x0;
    v[3] = y0;
}

//default
Hyperbola::Hyperbola() {v = new double[4];}

Hyperbola::Hyperbola(double *ptr_v){
/*Hyperbola from given array of length 4 with given parameters*/
    
    v = new double[4];
    memmove(v, ptr_v, sizeof(*v)*4);
}


Hyperbola::~Hyperbola(){
//destructor
    delete[] v;
};


Hyperbola::Hyperbola::Hyperbola(const Hyperbola& a){
/*copy-constructor: create Hyperbola from given*/
    v = new double[4];
    memmove(v, a.v, sizeof(*v)*4);
}

//input n-th coordinate of vector
//create and input vector of given size;
Hyperbola Hyperbola::create_input(){
    Hyperbola a;
    a.input();
    return a;
}

void Hyperbola::input(){
    cout << "Enter a, b, x0, y0 parameters of hyperbola (y-y0)^2/b^2 - (x-x0)^2/a^2 = 1: "<< endl;
    for(int i = 0; i < 4; i++){
        cin >> v[i];
    } 
}

//block of setters and getters
double Hyperbola::get_a() {return v[0];}
void Hyperbola::set_a(double _a) {this->v[0] = _a;}

double Hyperbola::get_b() {return v[1];}
void Hyperbola::set_b(double _b) {this->v[1] = _b;}

double Hyperbola::get_x0() {return v[2];}
void Hyperbola::set_x0(double _x0) {this->v[2] = _x0;}

double Hyperbola::get_y0() {return v[3];}
void Hyperbola::set_y0(double _y0) {this->v[3] = _y0;}  

bool Hyperbola::equal(const Hyperbola& other) {return v[0] == other.v[0] &&  v[1] == other.v[1] && v[2] == other.v[2] && v[3] == other.v[3];}

void Hyperbola::print(ostream& out = cout) {
    /*output of hyperbola*/
    
    out << "Hyperbola (y-("<< -v[3] <<"))^2/("<< v[1] <<")^2 - (x-("<< -v[2] <<"))^2/("<< v[0] <<")^2  = 1"<<ends;
}

ostream& operator<<(ostream& out, Hyperbola& h) {
    h.print(out);
    return out;
}

Point Hyperbola::get_center(){
    /*getting center of hyperbola*/
    Point p(v[2], v[3]);
    return p;
}


bool Hyperbola::is_outside_hyperbola(Point p){
    /*checks whether point lies outside given hyperbola
    (above upper part, below lower part)*/
    return (pow(p.y - v[3], 2) / pow(v[1], 2) - pow(p.x - v[2], 2) / pow(v[0], 2) >= 1);
}


// ////Intersection points of Hyperbola with given Parabola
    //Intersection points of ellipses (found numerically, rounded off)
void Hyperbola::intersect_parabola(Parabola p, ostream& out = cout){
    
    Point point;
    double h = 0.001;
    double min_x = -100;// not universal method
    double max_x = 100;//
    int num = 0; 


    point.x = min_x - h;
    point.y = pow(point.x - p.get_x0(), 2)/2/p.get_p() + p.get_y0();
    bool check_last = this->is_outside_hyperbola(point);
    bool check = check_last; 

    
    for (float x = min_x; x <= max_x; x += h){
        point.x = x;
        point.y = pow(point.x - p.get_x0(), 2)/2/p.get_p() + p.get_y0(); 
        check = this->is_outside_hyperbola(point);

        if (check != check_last){  
            out << point << endl;
            num += 1;
        }
        check_last = check;
    }

    if (!num){
        out << NAN << endl;
    }
    
}
    

void Hyperbola::intersect_hyperbola(Hyperbola other, ostream& out = cout){
    
    Point p;
    double h = 0.001;
    double max_y = 100;
    double min_y = -100;
    int num = 0; 


    p.y = min_y;  // other.v[3]+other.v[1]  // y0 + b
    p.x = other.v[0]/other.v[1]*sqrt(pow(p.y - other.v[3], 2) - v[1]*v[1]) + other.v[2];
    bool u_check_last = this->is_outside_hyperbola(p); 
    bool l_check_last = u_check_last; //setting check 
    
    bool u_check = l_check_last;
    bool l_check = u_check;


    // cout << other.v[0] << " " <<other.v[1] << endl;
    // cout << v[1]*v[1]  << other.v[2]

    if (this->equal(other)){
        out << INFINITY << endl;
    } else if ((v[1] + v[3]> other.v[1] + other.v[3]) && (-v[1] + v[3] < -other.v[1] + other.v[3]) &&
                v[0]/v[1] <= other.v[0]/other.v[1] ) {
        out << NAN << endl;
    }  else if (v[2] == other.v[2] && v[3] == other.v[3] && v[1] == other.v[1]) {
        p.x = v[2]; 
        p.y =  v[1] + v[3];
        out << p << endl;
        p.y = -v[1] + v[3];
        out << p << endl;
    } else {  
        for (float y = min_y; y <= max_y; y += h){
            if ((y < other.v[3]+other.v[1]) && (y > - other.v[3]+other.v[1])){ continue; }

            p.y = y;
            p.x = other.v[0]/other.v[1]*sqrt(pow(y - other.v[3], 2) - other.v[1]*other.v[1]) + other.v[2]; 
            u_check = this->is_outside_hyperbola(p);

            if (u_check != u_check_last){  
                out << p << endl;
                num += 1;
            }
            u_check_last = u_check;
            
            p.x = -other.v[0]/other.v[1]*sqrt(pow(y - other.v[3], 2) - other.v[1]*other.v[1]) + other.v[2]; 

            l_check = this->is_outside_hyperbola(p);
            if (l_check != l_check_last){  
                out << p << endl;
                num += 1;
            }
            l_check_last = l_check;
        }

        if (!num){
            out << NAN << endl;
        }
    }
}
