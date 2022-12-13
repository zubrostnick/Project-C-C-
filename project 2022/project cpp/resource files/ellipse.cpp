#include "ellipse.h"


//Задання еліпса у пдск, осі якого парлельні осям координат
//canonical ellipse eqation (x-x0)^2/a^2 + (y-y0)^2/b^2 = 1
void Ellipse::init(Point f1, Point f2, double r1, double r2) {
    if (r1 > 0 && r2 > 0){
        if (f1.y == f2.y){
            double c = abs(f1.x-f2.x)/2;
            this->a = (r1 + r2)/2; //using formula r1 = a + ex0; r2 = a-ex0 (M(x0,y0) - lies on ellipse)
            this->b = sqrt(a*a - c*c);
            this->x0 = max(f1.x, f2.x) - c;
            this->y0 = f1.y;
        }
        else if(f1.x == f2.x){
            double c = abs(f1.y-f2.y)/2;
            this->b = (r1 + r2)/2; 
            this->a = sqrt(b*b - c*c);
            this->y0 = max(f1.y, f2.y) - c;
            this->x0 = f1.x;
        } else throw ConicError("Inappropriate Focuses input! (Not for standard equation)");
            
    } else throw ConicError("Inappropriate Raduis input");
}



//creating Ellipse with given paremeters задання еліпса координатами фокусів еліпса та 2 радіусами
Ellipse::Ellipse(
    Point f1,
    Point f2,
    double r1 = 5,
    double r2 = 5
) {init(f1, f2, r1, r2);}

//default
Ellipse::Ellipse() {}

//destructor
Ellipse::~Ellipse(){
    //delete[] v;
};

//copy-constructor: create Ellipse from given
Ellipse::Ellipse(const Ellipse& e){
    a = e.a;
    b = e.b;
    x0 = e.x0;
    y0 = e.y0;
}

//create and input Ellipse
Ellipse Ellipse::create_input(){
    cout << "Creating Ellipse"<<endl;
    Ellipse a;
    a.input();
    return a;
}

void Ellipse::input(){
    double r1, r2;
    Point f1, f2;
    cout << "Enter coordinates of F1: "<< endl;
    cin >> f1;
    cout << "Enter coordinates of F2: "<< endl;
    cin >> f2;
    cout << "Enter 2 ellipse radiuses: "<< endl;
    cin >> r1;
    cin >> r2; 
    init(f1, f2, r1, r2);
}

//output of ellipse 
void Ellipse::print(ostream& out = cout) {
    out << "Ellipse (x+("<< -x0 <<"))^2/(" << a << ")^2 + ((y+("<< -y0 <<"))^2/(" << b <<")^2 = 1"<<ends;
}


ostream& operator<<(ostream& out, Ellipse& e) {
    e.print(out);
    return out;
}

//perfomes direct comparison between two Ellipses 
bool Ellipse::equal(const Ellipse& other) {return a == other.a &&  b == other.b && x0 == other.x0 && y0 == other.y0;}
//bool operator==(const Ellipse& other) {return a == other.a &&  b == other.b && x0 == other.x0 && y0 == other.y0};//return this.equal(other);}
//bool operator!=(const Ellipse& other) {return a != other.a ||  b != other.b || x0 != other.x0 || y0 != other.y0;}


//block of setters and getters
double Ellipse::get_a() {return a;}
void Ellipse::set_a(double _a) {this->a = _a;}

double Ellipse::get_b() {return b;}
void Ellipse::set_b(double _b) {this->b = _b;}

double Ellipse::get_x0() {return x0;}
void Ellipse::set_x0(double _x0) {this->x0 = _x0;}

double Ellipse::get_y0() {return y0;}
void Ellipse::set_y0(double _y0) {this->y0 = _y0;}  

Point Ellipse::get_f1() {
    Point f1;
    double c;
    if (a >= b){
        c = sqrt(a*a - b*b);
        f1.x = x0-c;
        f1.y = y0;
    } else{
        c = sqrt(b*b - a*a);
        f1.x = x0;
        f1.y = y0-c;
    }
    return f1;
}

Point Ellipse::get_f2() {
    Point f2;
    double c;
    if (a >= b){
        c = sqrt(a*a - b*b);
        f2.x = x0+c;
        f2.y = y0;
    } else{
        c = sqrt(b*b - a*a);
        f2.x = x0;
        f2.y = y0+c;
    }
    return f2;
}

// void set_f1(double _x, double _y) {
//     this->a2 = _a2;
// }


double Ellipse::area(){
    return 3.141592*a*b;
}

//getting center of ellipse
Point Ellipse::get_center(){
    Point p(x0, y0);
    return p;
}

//checks whether point lies inside given ellipse
bool Ellipse::is_inside_ellipse(Point p){
    return pow(p.x +(-x0), 2)/(a*a) + pow(p.y+(-y0), 2)/(b*b) <= 1;
}

//Intersection points of ellipse with given line, given by equation y = kx + h
void Ellipse::intersect_line( double k = 1, double h = 0, ostream& out = cout){
    /*x^2*(b^2 + a^2*k^2)+x(2a^2*k*h*y0 - 2b^2*x0) + (b^2*x0^2 + a^2*h^2*y0^2-a^2*b^2) = 0: result - intersection points*/
    //static Point res[2];    
    
    ostringstream sout;
    double val;
    int num = 0;
    Quadratic_eq e(b*b + a*a*k*k, 2*a*a*k*(h-y0) - 2*b*b*x0, b*b*x0*x0 + a*a*(h-y0)*(h-y0) - a*a*b*b);
    //cout << e << endl;
    Point p;

    e.solve(sout);
    istringstream i_res(sout.str());

    while ((i_res >> val)) {
        num++;
        p.x = val;
        p.y = k*p.x + h;
        out << p << endl;   
    }
        if (sout.str().find("nan") != string::npos) {  //inf is impossible in this case 
        out <<sout.str() << ends;
        // std::cout << "found nan!" << '\n';
    } else {
        if (num == 0){  //(!num)
            out << NAN << endl;
        }  
    }   
}


//Intersection points of ellipses (found numerically, rounded off)
void Ellipse::intersect_ellipse(Ellipse e, ostream& out = cout){
    
    Point p;
    double h = 0.001;
    double min_x = e.x0 - e.a;
    double max_x = e.x0 + e.a;
    int num = 0; 


    p.x = min_x - h;
    p.y = e.y0;
    bool u_check_last = this->is_inside_ellipse(p);  //setting check for upper part of the ellipse 
    bool l_check_last = u_check_last; //setting check for lower part of the ellipse

    bool u_check = u_check_last;
    bool l_check = l_check_last;


    //y = +-b/a*sqrt(a*a-(x-x0)*(x-x0)) + y0;
    if (this->equal(e)){
        out << INFINITY << endl;
    
    } else {  
        for (float x = min_x; x <= max_x; x += h){
            p.x = x;
            p.y = e.b/e.a*sqrt(e.a*e.a-pow(x - e.x0, 2)) + e.y0; 
            u_check = this->is_inside_ellipse(p);

            if (u_check != u_check_last){  
                /*if (u_check == true){
                    p.x = x;
                } else p.x = x - h;*/
                //p.x = (x + x - h)/2;
                //p.y= e.b/e.a*sqrt(e.a*e.a-pow(p.x - e.x0, 2)) + e.y0;
                out << p << endl;
                num += 1;
            }
            u_check_last = u_check;
            
            p.y = -e.b/e.a*sqrt(e.a*e.a-pow(x - e.x0, 2)) + e.y0; 
            l_check = this->is_inside_ellipse(p);

            if (l_check != l_check_last){  
                out << p << endl;
                num += 1;
            l_check_last = l_check;

            }

        }

        if (!num){
            out << NAN << endl;
        }
    }
}