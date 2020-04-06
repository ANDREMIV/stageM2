#include <stdio.h>
#include <stdlib.h>
#include "mathutils.h"
#include <math.h>
///IN SI
const double C=2.99792458e8;
const double KB=1.380649e-23;
const double MH=1.67265e-27;
const double Me=9.1093837015e-31;
const double TS=6.6524587158e-29;
const double AR=7.5657e-16;
const double hPl=6.62606876e-34;
const double G=6.674e-11;
const double H0=6.732e1*1e3/1e6/3.08567758149e16;
const double Ol=0.6847,Omo=0.3153;
const double Oro=9.877e-7;
const double Trado=2.7255;
const double TI=1.579e5; //13.6eV/kb in Kelvin
const double zeq=3402;//+-26 plank 2018


double expansion(double t, double a)
{
    return  sqrt(Ol*a*a+Omo/a+Oro/a/a);
}

/*void derivs(double t, double*y, double*dydt)//n=4
//y0=a
//y1=Tbar
//y2=Trad
//y3=a'
{
    double DC=8/3*TS*AR/Me/C/H0; //Compton coupling constant
    dydt[0]=y[3];
    dydt[1]=-2*y[3]/y[0]*y[1]\
            +DC*y[2]*y[2]*y[2]*y[2]*(y[2]-y[1])*1;//*XE(TI/y[2]);
    dydt[2]=-y[2]*y[3]/y[0];
    dydt[3]=+1.0/2*(2*Ol*y[0]-Omo/y[0]/y[0]-2*Oro/y[0]/y[0]/y[0]);
}*/
void derivs2(double t, double*y, double*dydt)//n=3
//y0=a
//y1=Tbar
//y2=Trad
{
    double DC=8/3*TS*AR/Me/C/H0; //Compton coupling constant
    dydt[0]=expansion(0,y[0]);
    /*double x= (-2*dydt[0]/y[0]*y[1]\
            +DC*y[2]*y[2]*y[2]*y[2]*(y[2]-y[1])*1e-4);
    double xx= -y[2]*dydt[0]/y[0];
    dydt[1]=x > xx ? xx : x; ///Overcompensation protection
    if(dydt[1]==-y[2]*dydt[0]/y[0])
        {
            int qsdsq;
    qsdsq+=5;
            ;}//debug purposes*/
    dydt[1]=-2*dydt[0]/y[0]*y[1]\
            +DC*y[2]*y[2]*y[2]*y[2]*(y[2]-y[1])*4e-4;
    dydt[2]=-y[2]*dydt[0]/y[0];
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



int expansion_calc(double **t, double **a, double **Tr,double **Tb, double (*xe)(double Tbar))
{
    FILE* results=fopen("expansion.txt","w");
    int i, n = MEMORYBLOC, m=NBMEMBLOC;
    double t0 = 0; //t0 = 0 when we are in zeq
    double dt=1e-17;
    **a=1.0/(zeq+1);
    **Tr=Trado*(zeq+1);
    **Tb=**Tr;
    **t=t0;
    double DT=1e-3;
    fprintf(results,"t'\ta(t')\tTrad(t')\tTb(t')\tz(t')\ti\n------------\n");

    i=0;
    double yo[3]; //initial condition for RK
    yo[0]=D(a,i);
        yo[1]=D(Tb,i);
        yo[2]=D(Tr,i);
        double tc=D(t,i); //current time
        double Oy[3]; //output of RK

    for (i = 0; D(a,i)<1 && i<n*m-1; i++)
    {
        if(!((i+1)%n))
        {
            (*(t+(i+1)/n))=(double *)malloc(sizeof(double) * (n));
            (*(a+(i+1)/n))=(double *)malloc(sizeof(double) * (n));
            (*(Tr+(i+1)/n))=(double *)malloc(sizeof(double) * (n));
            (*(Tb+(i+1)/n))=(double *)malloc(sizeof(double) * (n));
        }

        int j;
        for(j=0;j<10000;j++){
        rk62(derivs2,3,dt,tc,yo,Oy);

        yo[0]=Oy[0];
        yo[1]=Oy[1];
        yo[2]=Oy[2];
        tc=tc+dt;}



        D(a,i+1)=Oy[0];
        //if(Oy[1]>Oy[2])Oy[1]=Oy[2];//Delete positive surcompensation
        D(Tb,i+1)=Oy[1];
        D(Tr,i+1)=Oy[2];
        D(t,i+1)=tc;//+dt;

        /*double yo[4]; //initial condition for RK
        double Oy[4]; //output of RK
        yo[0]=D(a,i);
        yo[1]=D(Tb,i);
        yo[2]=D(Tr,i);
        yo[3]=expansion(0,D(a,i));
        double tc=D(t,i); //current time

        rk62(derivs,4,dt,tc,yo,Oy);

        D(a,i+1)=Oy[0];
        D(Tb,i+1)=Oy[1];
        D(Tr,i+1)=Oy[2];
        D(t,i+1)=tc+dt;*/

        //Variable step
        double TVAR=D(a,i+1)/D(a,i);
        if(TVAR>1+10*DT)
            dt/=2;
        else if(TVAR<1+DT)
            dt*=2;

        fprintf(results,"%.7le\t%.7le\t%.7le\t%.7le\t%.7le\t%d\n",\
                D(t,i), D(a,i),D(Tr,i),D(Tb,i),1/D(a,i)-1,i);
    }
    int imax=i;

    fclose(results);
    return imax;
}
