#include <stdio.h>
#include <stdlib.h>
#include "mathutils.h"
#include <math.h>
///IN SI
const double C=2.99792458e8;
const double KB=1.380649e-23;
const double MH=1.67265e-27;
const double TS=6.6524587158e-29;
const double AR=7.5657e-16;
const double hPl=6.62606876e-34;
const double G=6.674e-11;
const double H0=2.42e-18;
const double Ol=0.6847,Omo=0.3153;
const double Oro=9.877e-7;
const double Trado=2.7255;
const double TI=1.579e5; //13.6eV/kb in Kelvin
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
    double DC=8/3*TS*AR/MH/C/H0; //Compton coupling constant
    dydt[0]=y[1];
    dydt[1]=-1/2*(Ol*y[0]-Omo/y[0]/y[0]-2*Oro/y[0]/y[0]/y[0]);
    dydt[2]=-2*y[1]/y[0]*y[2]\
            +DC*y[3]*y[3]*y[3]*y[3]*fabs(y[3]-y[2])*1e-4;//*XE(TI/y[2]);
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



int expansion_calc(double **t, double **a, double **Trad,double **Tb)
{
    FILE* results=fopen("expansion.txt","w");
    int i, n = MEMORYBLOC, m=NBMEMBLOC;
    double t0 = 0, t1 = -0.96, dt=-1e-4;
    double DT=fabs(dt);
    **a=1;
    **Trad=Trado;
    **t=t0;
    fprintf(results,"t'\ta(t')\tTrad(t')\tTb(t')\ti\n------------\n");
    for (i = 1; *(*(t+i/n)+(i%n)) > t1 && i<n*m; i++)
    {
        double tp=*(*(t+(i-1)/n)+((i-1)%n));
        double ap=*(*(a+(i-1)/n)+((i-1)%n));



        if(!(i%n))
        {
            (*(t+i/n))=(double *)malloc(sizeof(double) * (n));
            (*(a+i/n))=(double *)malloc(sizeof(double) * (n));
            (*(Trad+i/n))=(double *)malloc(sizeof(double) * (n));
            (*(Tb+i/n))=(double *)malloc(sizeof(double) * (n));
        }
        *(*(a+i/n)+(i%n)) = rk4(expansion, dt, tp, ap);

        double A=*(*(a+i/n)+(i%n));
        *(*(Trad+i/n)+(i%n)) = Trado/A;
        //Variable step
        double TVAR=fabs(ap)/A;
        if(TVAR>1+5*DT)
            dt/=2;
        else if(TVAR<1+0.01*DT && TVAR>1)
            dt*=1.5;
        *(*(t+i/n)+(i%n))=tp+dt;

        if(fabs(dt)<1e-12)
            break;
    }
    int imax=i;



    //Initial condition for Tb = Tr
    for (i = 0; i < imax; i++)
        if(1/(*(*(a+i/n)+(i%n)))-1>=zeq)
            break;
    int ieq=i;
    *(*(Tb+ieq/n)+(ieq%n))=*(*(Trad+ieq/n)+(ieq%n));



    for (i = ieq; i > 0; i--)
    {
        double TRAD=*(*(Trad+i/n)+(i%n));
        double A=*(*(a+i/n)+(i%n));
        double TB=D(Tb,i);
        double tr=*(*(t+i/n)+(i%n));
        double dt=(*(*(t+(i-1)/n)+((i-1)%n))-tr);
        double yo[4];
        double Oy[4];
        yo[0]=A;
        yo[1]=expansion(tr,A);
        yo[2]=TB;
        yo[3]=TRAD;
        rk62(derivs,4,dt,tr,yo,Oy);
        *(*(Tb+(i-1)/n)+((i-1)%n))=Oy[2];

    }
    for (i = ieq; i < imax-1; i++)
    {
        double TRAD=*(*(Trad+i/n)+(i%n));
        double A=*(*(a+i/n)+(i%n));
        double TB=*(*(Tb+i/n)+(i%n));
        double tr=*(*(t+i/n)+(i%n));
        double dt=(*(*(t+(i+1)/n)+((i+1)%n))-tr);
        double yo[4];
        double Oy[4];
        yo[0]=A;
        yo[1]=expansion(tr,A);
        yo[2]=TB;
        yo[3]=TRAD;
        rk62(derivs,4,dt,tr,yo,Oy);
        *(*(Tb+(i+1)/n)+((i+1)%n))=Oy[2];

    }
    for (i = 0; i < imax; i++)
        fprintf(results,"%.13lf\t%.13lf\t%.13lf\t%.13lf\t%d\n",\
                *(*(t+i/n)+(i%n)), *(*(a+i/n)+(i%n)),*(*(Trad+i/n)+(i%n)),*(*(Tb+i/n)+(i%n)),i);
    fclose(results);
    return imax;
}
