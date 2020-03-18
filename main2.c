#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define D(pointer,index) ((D(pointer,(index))))
///t is the cosmic time t multiplied by Hubble constant Ho
///a is the scale factor
///f(a,t)= da/dt

double rk4(double(*f)(double, double), double dt, double t, double a)
{
    double    k1 = dt * f(t, a),
        k2 = dt * f(t + dt / 2, a + k1 / 2),
        k3 = dt * f(t + dt / 2, a + k2 / 2),
        k4 = dt * f(t + dt, a + k3);
    return a + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

const double Ol=0.6847,Omo=0.3153;
const double Oro=9.877e-7;
const double Trado=2.7255;


double expansion(double t, double a)
{
    return  sqrt(Ol*a*a+Omo/a+Oro/a/a);
}

double XE(double Xi)
{
    return 1-erf(sqrt(Xi))+2*sqrt(Xi)*exp(-Xi)/sqrt(M_PI);
}
const double TI=1.579e5;
const double DC=2.23625675e-6; //double xe=0.001; //ionization fraction
const double zeq=3402;//+-26 plank 2018


FILE* results;

int main(void)
{
    #define MEMORYBLOC 1e5
    #define NBMEMBLOC 100
    results=fopen("res.txt","w");
    double **a,**Trad, **Tb, **t;
    int i, n = MEMORYBLOC, m=NBMEMBLOC;
    double t0 = 0, t1 = -0.96, dt=-1e-4;
    double DT=fabs(dt);


    t = (double **)malloc(sizeof(double*) * (m));
    *t = (double *)malloc(sizeof(double) * (n));
    a = (double **)malloc(sizeof(double*) * (m));
    *a = (double *)malloc(sizeof(double) * (n));
    Trad= (double **)malloc(sizeof(double) * (m));
    *Trad= (double *)malloc(sizeof(double) * (n));
    Tb= (double **)malloc(sizeof(double) * (m));
    *Tb= (double *)malloc(sizeof(double) * (n));
    **a=1;
    **Trad=Trado;
    **t=t0;
    fprintf(results,"t'\ta(t')\tTrad(t')\tTb(t')\n------------\n");
    for (i = 1; (D(t,i)) > t1 && i<n*m; i++){
            double tp=(D(t,(i-1)));
            double ap=(D(a,(i-1)));



        if(!(i%n)){
                (*(t+i/n))=(double *)malloc(sizeof(double) * (n));
                (*(a+i/n))=(double *)malloc(sizeof(double) * (n));
                (*(Trad+i/n))=(double *)malloc(sizeof(double) * (n));
                (*(Tb+i/n))=(double *)malloc(sizeof(double) * (n));
        }
        (D(a,i)) = rk4(expansion, dt, tp, ap);
        double A=(D(a,i));
        (D(Trad,i)) = Trado/A;
        double TVAR=fabs(ap)/A;
        if(TVAR>1+5*DT)
            dt/=2;
        else if(TVAR<1+0.1*DT && TVAR>1)
            dt*=1.5;
        (D(t,i))=tp+dt;
        if(fabs(dt)<1e-9)break;
        }
        int imax=i;



    for (i = 0; i < imax; i++)
        if(1/((D(a,i)))-1>=zeq)break;
    int ieq=i;
    (D(Tb+ieq/n)+(ieq%n))=*(*(Trad,ieq));

    double afunc(double T)
{   int j;
    for (j = 0; j < imax; j++)
           if(((D(t,j)))<=T)break;
        double ao=(D(a,j));
        double to=(D(t,j));
        double dt=T-to;
    return rk4(expansion, dt, to, ao);
}

double DTbar(double t,double Tbar)
{
    double A=afunc(t);
    double B=-2*(expansion(t,A))/A*Tbar;//sure
    double C=-DC*Trado*Trado*Trado*Trado*(Trado/A-Tbar)/A/A/A/A*XE(TI/Tbar);
    return B+C;
}

    for (i = ieq; i > 0; i--){
        double TB=D(Tb,i);
        double tr=(D(t,i));
        double dt=((D(t,(i-1)))-tr);
        (D(Tb,(i-1))) = rk4(DTbar, dt, tr, TB);
        }
    for (i = ieq; i < imax-1; i++){
        double TB=(D(Tb,i));
        double tr=(D(t,i));
        double dt=((D(t,(i+1)))-tr);
        (D(Tb,(i+1))) = rk4(DTbar, dt, tr, TB);
        }
    for (i = 0; i < imax; i++)
    fprintf(results,"%.13lf\t%.13lf\t%.13lf\t%.13lf\n",\
    (D(t+i/n)+(i%n)), *(*(a+i/n)+(i%n)),*(*(Trad+i/n)+(i%n)),*(*(Tb,i)));


    fclose(results);

    return 0;
}
