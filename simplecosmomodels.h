#include <stdio.h>
#include <stdlib.h>

extern const double G;
extern const double H0;
extern const double hPl;
extern const double AR;
extern const double MH;
extern const double Me;
extern const double TS;
extern const double KB;
extern const double C;
extern const double Ol,Omo;
extern const double Oro;
extern const double Trado;
extern const double TI; //13.6eV/kb in Kelvin
extern const double zeq;//+-26 plank 2018

void derivs(double t, double*y, double*dydt);
void derivs2(double t, double*y, double*dydt);
void desitter();
void einsteindesitter();
double expansion(double t, double a);
double BoltzmanH2ROP(double T);
int expansion_calc(double **t, double **a, double **Trad,double **Tb);
