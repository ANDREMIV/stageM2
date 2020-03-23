#include <stdio.h>
#include <stdlib.h>

extern const double G;
extern const double H0;
extern const double Ol,Omo;
extern const double Oro;
extern const double Trado;
extern const double TI; //13.6eV/kb in Kelvin
extern const double DC; //Compton coupling constant
extern const double zeq;//+-26 plank 2018

void derivs(double t, double*y, double*dydt);
void desitter();
void einsteindesitter();
double expansion(double t, double a);
