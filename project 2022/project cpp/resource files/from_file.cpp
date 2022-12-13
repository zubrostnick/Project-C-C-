#include "from_file.h"

void Equations_from_file::read_monomial(string monom, ostream& out = cout){
    string split("*x^");
    double val;
    int pos = monom.find(split);
    if (pos != -1){
        out << monom.substr(pos+split.length(), -1) <<endl;
        out << monom.substr(0, pos) <<endl;
    }
    else out << -1 << endl;
} 

Equations_from_file::Equations_from_file(const string& filename): filename(filename) {}

void Equations_from_file::exercise () {
    string biq_end("*x^4");    //for detecting type of equation
    string q_end("*x^2");
    string l_end("*x");
    string line;            //for reading one line (equation);
    double a;            //for transfering coef. to equations

    // block of values for statistics - fuck this shit, im out; 
    int n_nan = 0;
    int n_inf = 0;
    // no more spending time on it, it is worthless (only 3 points)

    ifstream finp(filename);
    if (finp.bad()) throw logic_error("");
    
    Linear_eq l_eq;
    Quadratic_eq q_eq;
    Biquadratic_eq biq_eq;
    int type; // for detecting the type of equation

    while (getline(finp, line)) {

        if (finp.fail()) break;
        //cout<< line << endl;
        istringstream line_s(line);

        if (line.find(biq_end) != -1){// if equation is biquadratic
            type = 2;
        } else if(line.find(q_end) != -1){ // if quation is quadratic
            type = 1;
        } else if(line.find(l_end) != -1){ // if quation is linear
            type = 0;
        } else continue; // not an equation

        string monomial; 
        //searching, setting of corresponding coefficients
        while ((line_s >> monomial)){                //monomial
            ostringstream read_m;
            read_monomial(monomial, read_m);
            istringstream coef_monom(read_m.str());
            coef_monom >> a;  // exponent

            if (a == 4){
                coef_monom >> a;
                biq_eq.set_a2(a);
            } else if (a == 2){
                coef_monom >> a;
                if (type == 2){  // if biquadratic
                    biq_eq.set_a1(a);
                } else {  // if quadratic
                    q_eq.set_a2(a);
                }
            } else if (a == 1){
                coef_monom >> a;
                if (type == 1){  // if quadratic
                    q_eq.set_a1(a);
                } else {  // if linear
                    l_eq.set_a1(a);
                }  
            } else if (a == 0){   
                coef_monom >> a;
                if (type == 2){    //if biquadratic
                    biq_eq.set_a0(a);
                } else if (type == 1){   // if quadratic
                    q_eq.set_a0(a);
                }  else {       // if linear
                    l_eq.set_a0(a);
                    //cout << a << endl;
                }
            }

        } 
        //solving, collecting data     
        if (type == 2){    //if biquadratic
            cout <<  biq_eq << endl;
            cout << "Solvings: "<< endl;
            biq_eq.solve();

            //setting to default after all the manipalations
            biq_eq.set_a2(0);
            biq_eq.set_a1(0);
            biq_eq.set_a0(0);
        } else if (type == 1){   // if quadratic
            cout <<  q_eq << endl;
            cout << "Solvings: "<< endl;
            q_eq.solve();

            q_eq.set_a2(0);
            q_eq.set_a1(0);
            q_eq.set_a0(0);
        }  else {       // if linear
            //cout << "check" << endl;
            cout <<  l_eq << endl;
            cout << "Solvings: "<< endl;
            l_eq.solve();
            
            l_eq.set_a1(0);
            l_eq.set_a0(0);
        }    
        cout << endl;
    }
    finp.close();
}