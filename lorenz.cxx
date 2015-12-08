#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

void func(double* y, const int a, const int b, const double c);
     
int main() {
    double y[3];   //x`=y(0), y`=y(1), z`=y(2)
    y[0] = 1; //Startwerte
    y[1] = 1;
    y[2] = 1;
    double tz = 1.0; //Zeitendwert
    int N = 10; //Anzahl der Schritte
    double dt = tz/N; //Schrittweite
    double k1[3];  //4k mit 3 Dimensionen
    double k2[3];
    double k3[3];
    double k4[3];
    double temp[3];
    const int a = 10;
    const int b = 28;
    const double c = 8./3;
    
    ofstream out("lorenzmodell.txt");
    out << 0 << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << endl;
    
    for (int i = 1; i < N; i++) {
           
        //k1
        temp[0] = y[0];
        temp[1] = y[1];
        temp[2] = y[2];
        func(temp, a, b, c);
        k1[0] = temp[0];
        k1[1] = temp[1];
        k1[2] = temp[2];
        cout << dt*i <<"\t"<< k1[0] << endl;
         cout << dt*i <<"\t"<< k1[1] << endl;
          cout << dt*i <<"\t"<< k1[2] << endl;
        //k2
        temp[0] = y[0] + (dt/2)*k1[0];
        temp[1] = y[1] + (dt/2)*k1[1];
        temp[2] = y[2] + (dt/2)*k1[2];
        func(temp, a, b, c);
        k2[0] = temp[0];
        k2[1] = temp[1];
        k2[2] = temp[2];
        
        //k3
        temp[0] = y[0] + (dt/2)*k2[0];
        temp[1] = y[1] + (dt/2)*k2[1];
        temp[2] = y[2] + (dt/2)*k2[2];
        func(temp, a, b, c);
        k3[0] = temp[0];
        k3[1] = temp[1];
        k3[2] = temp[2];
        
        //k4
        temp[0] = y[0] + dt*k3[0];
        temp[1] = y[1] + dt*k3[1];
        temp[2] = y[2] + dt*k3[2];
        func(temp, a, b, c);
        k4[0] = temp[0];
        k4[1] = temp[1];
        k4[2] = temp[2];
        
        for (int j = 0; j <3; j++){
            y[j] += dt*((1./6)*k1[j]+(1./3)*k2[j]+(1./3)*k3[j]+(1./6)*k4[j]);
        }
        
        out << dt*i <<"\t" << y[0] << "\t" << y[1] << "\t" << y[2]<< endl;
    
    }
    out.close();
    return 0;
}
    
    
    
//lorenzmodell
void func(double* y, const int a, const int b, const double c){
    double* t = y;
    y[0] = a*(t[1]-t[0]);
    y[1] = t[0]*(b-t[2])-t[1];
    y[2] = t[0]*t[1]-c*t[2];
}