#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

void func(double* y, double* k, const int a, const int b, const double c);
     
int main() {
    double y[3];   //x`=y(0), y`=y(1), z`=y(2)
    y[0] = 1; //Startwerte
    y[1] = 1;
    y[2] = 1;
    double tz = 100; //Zeitendwert
    int N = 10000; //Anzahl der Schritte
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
        func(temp, k1, a, b, c);
    
        //k2
        temp[0] = y[0] + (dt/2)*k1[0];
        temp[1] = y[1] + (dt/2)*k1[1];
        temp[2] = y[2] + (dt/2)*k1[2];
        func(temp, k2, a, b, c);

        
        //k3
        temp[0] = y[0] + (dt/2)*k2[0];
        temp[1] = y[1] + (dt/2)*k2[1];
        temp[2] = y[2] + (dt/2)*k2[2];
        func(temp, k3, a, b, c);
        
        //k4
        temp[0] = y[0] + dt*k3[0];
        temp[1] = y[1] + dt*k3[1];
        temp[2] = y[2] + dt*k3[2];
        func(temp, k4, a, b, c);
    
        //yn+1
        for (int j = 0; j <3; j++){
            y[j] += dt*((1./6)*k1[j]+(1./3)*k2[j]+(1./3)*k3[j]+(1./6)*k4[j]);
        }
        
        out << dt*i <<"\t" << y[0] << "\t" << y[1] << "\t" << y[2]<< endl;
    
    }
    out.close();
    return 0;
}
    
//lorenzmodell
void func(double* y, double* k, const int a, const int b, const double c){
    //k1 = f(y)...
    k[0] = a*(y[1]-y[0]);
    k[1] = y[0]*(b-y[2])-y[1];
    k[2] = y[0]*y[1]-c*y[2];
}