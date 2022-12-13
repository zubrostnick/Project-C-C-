#include <iostream>
#include <string>
#include <fstream>
#include <cstring> //memmove, memcpy
#include <istream>
#include <ostream>
#include <sstream>
#include <cmath>
//#include "Equations.h"
using namespace std;
#define _USE_MATH_DEFINES



class Linear_eq {
    //a1 * x + a0 = 0
    double a0;
    double a1;
     
public:

    Linear_eq(double _a1 = 0, double _a0 = 0): a0(_a0), a1(_a1) {}

    double get_a0() const {return a0;}
    void set_a0(double a) {this->a0 = a;}
    double get_a1() const {return a1;}
    void set_a1(double a) {this->a1 = a;}

    //virtual: Буде викликано в залежності від поточного стану програми
    virtual void input(istream& inp = cin) {
        inp >> a1 >> a0;
    }

    virtual void print(ostream& out = cout) const {
        out <<"("<< a1 << ")*x + (" << a0 << ") = 0";
    }

    //Перевантаження оператору в загальному розумінні (як зовнішню функцію, а не внутрішній метод)
    friend ostream& operator<<(ostream& out, const Linear_eq& eq) {
        eq.print(out);
        return out;
    }

    friend istream& operator>>(istream& inp, Linear_eq& eq) {
        eq.input(inp);
        return inp;
    }

    //output of solvings in stream 
    void solve(ostream& out = cout){
        //cout << "check" << endl;
        if ((a1 == 0) && (a0 != 0)){
            //throw logic_error("Unsolvable");
            out << NAN << endl;
        } else if ((a1 == 0) && (a0 == 0)){
            out << INFINITY << endl;
        } else out << -a0/a1 << endl;
    }
};


class Quadratic_eq: public Linear_eq {
    //a2*x^2 + a1 * x + a0 = 0

    double a2;

public:

    Quadratic_eq(
        double _a2 = 0,
        double _a1 = 0,
        double _a0 = 0
    ):  Linear_eq(_a1, _a0), a2(_a2) {}

    double get_a2() const {return a2;}
    void set_a2(double _a2) {this->a2 = _a2;}

    virtual void input(istream& inp = cin) {
        inp >> a2;
        Linear_eq::input(inp);
    }

    virtual void print(ostream& out = cout) const {
        out << a2 <<"*x^2 + ";
        Linear_eq::print(out);
    }

    //output of solvings in stream 
    void solve(ostream& out = cout){
        double res;
        if (a2 == 0){
            ostringstream check;
            Linear_eq :: solve (out);
            //Linear_eq :: solve(check);
            //cout << check.str() << endl;     
        }
        else if (get_a1() == 0) {
            if (-get_a0()/a2 < 0) {
                out << NAN << endl;
            }
            else out << sqrt(-get_a0()/a2) << endl << -sqrt(-get_a0()/a2) << endl;
        }
        else if (get_a0() == 0){
            out << 0 << endl;
            Linear_eq eq_q(a2, get_a1());
            eq_q.solve(out); // Але тут з відповіддю 0, може стояти inf, або nan 
        }
        else{
           
            double D = (get_a1()*get_a1() - 4*a2*get_a0());
            if (D < 0){
                out << NAN <<endl;
            } else if (D == 0){
                out << -get_a1()/2/a2 << endl;
            }
            else {
                out << (-get_a1() + sqrt(D))/2/a2 << endl << (-get_a1() - sqrt(D))/2/a2 <<endl;
            }
        }
        
    } 

};


struct Point{
/*Structure that creates 2d point + operans reassigment*/
    double x, y;

    Point(double _x, double _y): x(_x), y(_y) {}
   
    //default
    Point() {}

    virtual void input(istream& inp = cin) {
        inp >> x >> y;
    }
    
    virtual void print(ostream& out = cout) const {
        out << "Point(" << x << "," << y <<")"<<ends;
    }


    friend ostream& operator<<(ostream& out, const Point& p) {
        p.print(out);
        return out;
    }

    friend istream& operator>>(istream& inp, Point& p) {
        p.input(inp);
        return inp;
    }

    //creates and returns a point 
    static Point create_input(){
        Point p;
        cout << "Input of point, please write in x, y coordinates of a point:" << endl;
        p.input();
        return p;
    }
};


class ConicError: public logic_error {
    public: ConicError(const string& message): logic_error(message) {}
};




//Задання еліпса у пдск, осі якого парлельні осям координат
class Ellipse {
//canonical ellipse eqation (x-x0)^2/a^2 + (y-y0)^2/b^2 = 1
    double a, b, x0, y0;
    
    void init(Point f1, Point f2, double r1, double r2) {
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

public:

    //creating Ellipse with given paremeters задання еліпса координатами фокусів еліпса та 2 радіусами
    Ellipse(
        Point f1,
        Point f2,
        double r1 = 5,
        double r2 = 5
    ) {init(f1, f2, r1, r2);}

    //default
    Ellipse() {}
    
    //destructor
    ~Ellipse(){
        //delete[] v;
    };

    //copy-constructor: create Ellipse from given
    Ellipse(const Ellipse& e){
        a = e.a;
        b = e.b;
        x0 = e.x0;
        y0 = e.y0;
    }

    //create and input Ellipse
    static Ellipse create_input(){
        cout << "Creating Ellipse"<<endl;
        Ellipse a;
        a.input();
        return a;
    }

    void input(){
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
    virtual void print(ostream& out = cout) const {
        out << "Ellipse (x+("<< -x0 <<"))^2/(" << a << ")^2 + ((y+("<< -y0 <<"))^2/(" << b <<")^2 = 1"<<ends;
    }


    friend ostream& operator<<(ostream& out, const Ellipse& e) {
        e.print(out);
        return out;
    }

    //perfomes direct comparison between two Ellipses 
    bool equal(const Ellipse& other) {return a == other.a &&  b == other.b && x0 == other.x0 && y0 == other.y0;}
    //bool operator==(const Ellipse& other) const {return a == other.a &&  b == other.b && x0 == other.x0 && y0 == other.y0};//return this.equal(other);}
    //bool operator!=(const Ellipse& other) const {return a != other.a ||  b != other.b || x0 != other.x0 || y0 != other.y0;}


    //block of setters and getters
    double get_a() const {return a;}
    void set_a(double _a) {this->a = _a;}

    double get_b() const {return b;}
    void set_b(double _b) {this->b = _b;}

    double get_x0() const {return x0;}
    void set_x0(double _x0) {this->x0 = _x0;}

    double get_y0() const {return y0;}
    void set_y0(double _y0) {this->y0 = _y0;}  

    Point get_f1() const {
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

    Point get_f2() const {
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


    double area(){
        return M_PI*a*b;
    }

    //getting center of ellipse
    Point get_center(){
        Point p(x0, y0);
        return p;
    }

    //checks whether point lies inside given ellipse
    bool is_inside_ellipse(Point p){
        return pow(p.x +(-x0), 2)/(a*a) + pow(p.y+(-y0), 2)/(b*b) <= 1;
    }
    
    //Intersection points of ellipse with given line, given by equation y = kx + h
    void intersect_line( double k = 1, double h = 0, ostream& out = cout){
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
    void intersect_ellipse(Ellipse e, ostream& out = cout){
        
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
};


class Parabola {
/*using unverted standard formula for parabola: (x-x0)^2 = 2p(y-y0)
Parabola axis is parallel to Oy-axis*/

    double *v; //vector, that contains p, x0, y0

public:

    Parabola( double p, double x0, double y0) {
        /*creating parabola with given paremeters*/ 
        
        v = new double[3];
        v[0] = p;
        v[1] = x0;
        v[2] = y0;
    }

    //default
    Parabola() {v = new double[3];}

    Parabola(double *ptr_v){
    /*Parabola from given array of length 3 with given parameters*/
        
        v = new double[3];
        memmove(v, ptr_v, sizeof(*v)*3);
    }
    
    
    ~Parabola(){
    //destructor
        delete[] v;
    };


    Parabola(const Parabola& a){
    /*copy-constructor: create Parabola from given*/
        v = new double[3];
        memmove(v, a.v, sizeof(*v)*3);
    }

    //input n-th coordinate of vector
    //create and input vector of given size;
    static Parabola create_input(){
        Parabola a;
        a.input();
        return a;
    }

    void input(){
        cout << "Enter p, x0, y0 parameters of parabola (x-x0)^2 = 2p(y-y0): "<< endl;
        for(int i = 0; i < 3; i++){
            cin >> v[i];
        } 
    }
    
    //block of setters and getters
    double get_p() const {return v[0];}
    void set_p(double _p) {this->v[0] = _p;}

    double get_x0() const {return v[1];}
    void set_x0(double _x0) {this->v[1] = _x0;}

    double get_y0() const {return v[2];}
    void set_y0(double _y0) {this->v[2] = _y0;}  


    //output of parabola
    virtual void print(ostream& out = cout) const {
        out << "Parabola (x+("<< -v[1] <<"))^2 = " << v[0] << "*2*(y+("<< -v[2] <<"))"<<ends;
    }


    friend ostream& operator<<(ostream& out, const Parabola& p) {
        p.print(out);
        return out;
    }
    
    Point get_center(){
        /*getting center of parabola (focus)*/
        Point p(v[1], v[2]+v[0]/2);
        return p;
    }

    
    bool is_above_parabola(Point p){
        /*checks whether point lies above given parabola*/
        return (pow(p.x-v[1], 2) < 2*v[0]*(p.y-v[2]));
    }
    

    //Intersection points of Parabola with given line, given by equation y = kx + b
    void intersect_line(double k, double b, ostream& out = cout){
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

    double intersect_line_area( double k, double b){
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
    void intersect_parabola(Parabola p, ostream& out = cout){
        
        ostringstream sout;
        double val;
        int num = 0;
        Quadratic_eq e(p.v[0] - v[0],  2*(p.v[1]*v[0] - v[1]*p.v[0]),  
        p.v[0]*v[1]*v[1] - v[0]*p.v[1]*p.v[1] + 2*p.v[0]*v[0]*(v[2] - p.v[2]));
      
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

    double intersect_parabola_area( Parabola p){
        
        ostringstream sout;
        intersect_parabola(p, sout);
        double x1, x2;
        int num = 0;
        Point point;
        string str;

        istringstream i_res(sout.str());

        if (sout.str().find("nan") != string::npos) {  
           return 0;
        }
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
        } else {
            return abs((pow(x1, 3) - pow(x2, 3))*(p.v[0] - v[0])/3 + (x1*x1-x2*x2)*(p.v[1]*v[0] - v[1]*p.v[0]) +  
        (p.v[0]*v[1]*v[1] - v[0]*p.v[1]*p.v[1] + 2*p.v[0]*v[0]*(v[2] - p.v[2]))*(x1 - x2)) ;
        }  
         
    }
};


class Hyperbola {
/*using inverted standard formula for hyperbola: (y-y0)^2/b^2 - (x-x0)^2/a^2 = 1
Hyperbola axis is parallel to Oy-axis*/

    double *v; //vector, that contains a, b, x0, y0

public:

    Hyperbola(double a, double b, double x0, double y0) {
        /*creating parabola with given paremeters*/ 
        if (a<=0 || b<= 0)  throw ConicError("Inappropriate a, b parameters input");
        
        v = new double[4];
        v[0] = a;
        v[1] = b;
        v[2] = x0;
        v[3] = y0;
    }

    //default
    Hyperbola() {v = new double[4];}

    Hyperbola(double *ptr_v){
    /*Hyperbola from given array of length 4 with given parameters*/
        
        v = new double[4];
        memmove(v, ptr_v, sizeof(*v)*4);
    }
    
    
    ~Hyperbola(){
    //destructor
        delete[] v;
    };


    Hyperbola(const Hyperbola& a){
    /*copy-constructor: create Hyperbola from given*/
        v = new double[4];
        memmove(v, a.v, sizeof(*v)*4);
    }

    //input n-th coordinate of vector
    //create and input vector of given size;
    static Hyperbola create_input(){
        Hyperbola a;
        a.input();
        return a;
    }

    void input(){
        cout << "Enter a, b, x0, y0 parameters of hyperbola (y-y0)^2/b^2 - (x-x0)^2/a^2 = 1: "<< endl;
        for(int i = 0; i < 4; i++){
            cin >> v[i];
        } 
    }
    
    //block of setters and getters
    double get_a() const {return v[0];}
    void set_a(double _a) {this->v[0] = _a;}

    double get_b() const {return v[1];}
    void set_b(double _b) {this->v[1] = _b;}

    double get_x0() const {return v[2];}
    void set_x0(double _x0) {this->v[2] = _x0;}

    double get_y0() const {return v[3];}
    void set_y0(double _y0) {this->v[3] = _y0;}  

    bool equal(const Hyperbola& other) {return v[0] == other.v[0] &&  v[1] == other.v[1] && v[2] == other.v[2] && v[3] == other.v[3];}

    virtual void print(ostream& out = cout) const {
        /*output of hyperbola*/
        
        out << "Hyperbola (y-("<< -v[3] <<"))^2/("<< v[1] <<")^2 - (x-("<< -v[2] <<"))^2/("<< v[0] <<")^2  = 1"<<ends;
    }

    friend ostream& operator<<(ostream& out, const Hyperbola& h) {
        h.print(out);
        return out;
    }
    
    Point get_center(){
        /*getting center of hyperbola*/
        Point p(v[2], v[3]);
        return p;
    }

    
    bool is_outside_hyperbola(Point p){
        /*checks whether point lies outside given hyperbola
        (above upper part, below lower part)*/
        return (pow(p.y - v[3], 2) / pow(v[1], 2) - pow(p.x - v[2], 2) / pow(v[0], 2) >= 1);
    }
    
    
    // ////Intersection points of Hyperbola with given Parabola
        //Intersection points of ellipses (found numerically, rounded off)
    void intersect_parabola(Parabola p, ostream& out = cout){
        
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
        

    void intersect_hyperbola(Hyperbola other, ostream& out = cout){
        
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

};


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
    //Parabola p = Parabola::create_input(); // this one obliges entering values from keyboard
    Point p3 = Point(0, 3); // 
    cout << "Creating parabola y = x^2:" << endl;
    Parabola par(0.5, 0 , 0); //y = x^2
    cout << par << endl; 
    cout << "Center of Parabola:" << par.get_center() << endl;
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
    
        cout << "Enter testing mode: 1 - Ellipse testing, 2 - Parbola testing, 3 - Hyperbola testing, 4 - exit" << endl;
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