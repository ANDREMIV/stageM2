#include <stdio.h>
#include <stdlib.h>
#include "mathutils.h"
#include <math.h>

const double G=6.674e-11;
const double H0=2.42e-19;
const double Ol=0.6847,Omo=0.3153;
const double Oro=9.877e-7;
const double Trado=2.7255;
const double TI=1.579e5; //13.6eV/kb in Kelvin
const double DC=2.23625675e-6; //Compton coupling constant
const double zeq=3402;//+-26 plank 2018

double expansion(double t, double a)
{
    return  sqrt(Ol*a*a+Omo/a+Oro/a/a);
}

void derivs(double t, double*y, double*dydt)//n=3
//y0=a
//y1=da/dt'
//y2=Tbar
//y3=Trad
{
    dydt[0]=y[1];
    dydt[1]=-1/2*(Ol*y[0]-Omo/y[0]/y[0]-2*Oro/y[0]/y[0]/y[0]);
    dydt[2]=-2*y[1]/y[0]*y[2]\
            -DC*y[3]*y[3]*y[3]*y[3]*(y[3]-y[2])*XE(TI/y[2])*1;
    dydt[3]=-y[3]*y[1]/y[0];
}

void desitter()
{
    double expansion3(double t, double a)
    {
        return  a;
    }
    FILE* desit=fopen("expansionDS.txt","w");
    double a=1;
    double t0 = 0, t1 = -5, dt=-1e-4,t;
    for(t=t0; t>t1; t+=dt)
    {
        a=rk4(expansion3, dt, t, a);
        fprintf(desit,"%.13lf\t%.13lf\n",t+dt,a);
    }
    fclose(desit);
}
void einsteindesitter()
{
    double expansion3(double t, double a)
    {
        return  sqrt(1/a);
    }
    FILE* eindesit=fopen("expansionEDS.txt","w");
    double a=1;
    double t0 = 0, t1 = -0.96, dt=-1e-4,t;
    for(t=t0; t>t1; t+=dt)
    {
        a=rk4(expansion3, dt, t, a);
        fprintf(eindesit,"%.13lf\t%.13lf\n",t+dt,a);
    }
    fclose(eindesit);
}