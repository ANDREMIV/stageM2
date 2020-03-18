#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define D(pointer,index) (*(*(pointer+(index)/n)+((index)%n)))
///t is the cosmic time t multiplied by Hubble constant Ho
///a is the scale factor
///f(a,t)= da/dt

const double Ol=0.6847,Omo=0.3153;
const double Oro=9.877e-7;
const double Trado=2.7255;
const double TI=1.579e5; //13.6eV/kb in Kelvin
const double DC=2.23625675e-6; //Compton coupling constant
const double zeq=3402;//+-26 plank 2018

void rk42(void (*derivs)(double, double*, double*), int n\
          ,double dt, double t, double*yo, double *Oy)
{
    double *Y = malloc(sizeof(double)*n);
    double *dydt = malloc(sizeof(double)*n);
    double *k1 = malloc(sizeof(double)*n);
    double *k2 = malloc(sizeof(double)*n);
    double *k3 = malloc(sizeof(double)*n);
    double *k4 = malloc(sizeof(double)*n);

    int i;
    derivs(t,yo,dydt);
    for(i=0;i<n;i++) //first step
    k1[i]=dt*dydt[i];
    for(i=0;i<n;i++)
        Y[i]=yo[i]+k1[i]/2;
    derivs(t+dt/2,Y,dydt);
    for(i=0;i<n;i++) //second step
    k2[i]=dt*dydt[i];
    for(i=0;i<n;i++)
        Y[i]=yo[i]+k2[i]/2;
    derivs(t+dt/2,Y,dydt);
    for(i=0;i<n;i++) //third step
    k3[i]=dt*dydt[i];
    for(i=0;i<n;i++)
        Y[i]=yo[i]+k3[i];
    derivs(t+dt,Y,dydt);
    for(i=0;i<n;i++) //fourth step
    k4[i]=dt*dydt[i];
    for(i=0;i<n;i++) //final step
    Oy[i]=yo[i]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;

    free(dydt);
    free(Y);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
}

void rk62(void (*derivs)(double, double*, double*), int n\
          ,double dt, double t, double*yo, double *Oy)
{
    double *Y = malloc(sizeof(double)*n);
    double *dydt = malloc(sizeof(double)*n);
    double *k1 = malloc(sizeof(double)*n);
    double *k2 = malloc(sizeof(double)*n);
    double *k3 = malloc(sizeof(double)*n);
    double *k4 = malloc(sizeof(double)*n);
    double *k5 = malloc(sizeof(double)*n);
    double *k6 = malloc(sizeof(double)*n);
    static float a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,b31=3.0/40.0,b32=9.0/40.0,\
    b41=0.3,b42 = -0.9,b43=1.2,b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,\
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,b64=44275.0/110592.0,\
    b65=253.0/4096.0,c1=37.0/378.0,c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0;

    int i;
    derivs(t,yo,k1);
    for(i=0;i<n;i++)
        Y[i]=yo[i]+dt*(b21*k1[i]);
    derivs(t+a2*dt,Y,k2);
    for(i=0;i<n;i++)
        Y[i]=yo[i]+dt*(b31*k1[i]+b32*k2[i]);
    derivs(t+a3*2,Y,k3);
    for(i=0;i<n;i++)
        Y[i]=yo[i]+dt*(b41*k1[i]+b42*k2[i]+b43*k3[i]);
    derivs(t+a4*dt,Y,k4);
    for(i=0;i<n;i++)
        Y[i]=yo[i]+dt*(b51*k1[i]+b52*k2[i]+b53*k3[i]+b54*k4[i]);
    derivs(t+a5*dt,Y,k5);
    for(i=0;i<n;i++)
        Y[i]=yo[i]+dt*(b61*k1[i]+b62*k2[i]+b63*k3[i]+b64*k4[i]+b65*k5[i]);
    derivs(t+a6*dt,Y,k6);

    for(i=0;i<n;i++) //final step
    Oy[i]=yo[i]+dt*(c1*k1[i]+c3*k3[i]+c4*k4[i]+c6*k6[i]);

    free(dydt);
    free(Y);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(k5);
    free(k6);
}

double XE(double Xi)
{
    return 1-erf(sqrt(Xi))+2*sqrt(Xi)*exp(-Xi)/sqrt(M_PI);
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

double rk4(double(*f)(double, double), double dt, double t, double a)
{
    double    k1 = dt * f(t, a),
        k2 = dt * f(t + dt / 2, a + k1 / 2),
        k3 = dt * f(t + dt / 2, a + k2 / 2),
        k4 = dt * f(t + dt, a + k3);
    return a + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}




double expansion(double t, double a)
{
    return  sqrt(Ol*a*a+Omo/a+Oro/a/a);
}



FILE* results;

void desitter()
{
    double expansion3(double t, double a)
{
    return  a;
}
    FILE* desit=fopen("expansionDS.txt","w");
    double a=1; double t0 = 0, t1 = -5, dt=-1e-4,t;
    for(t=t0;t>t1;t+=dt)
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
    double a=1; double t0 = 0, t1 = -0.96, dt=-1e-4,t;
    for(t=t0;t>t1;t+=dt)
    {
        a=rk4(expansion3, dt, t, a);
        fprintf(eindesit,"%.13lf\t%.13lf\n",t+dt,a);
    }
    fclose(eindesit);
}

int main(void)
{
    #define MEMORYBLOC 1e5
    #define NBMEMBLOC 100
    desitter();
    einsteindesitter();

    results=fopen("expansion.txt","w");
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
    //fprintf(results,"t'\ta(t')\n------------\n");
    fprintf(results,"t'\ta(t')\tTrad(t')\tTb(t')\ti\n------------\n");
    for (i = 1; *(*(t+i/n)+(i%n)) > t1 && i<n*m; i++){
            double tp=*(*(t+(i-1)/n)+((i-1)%n));
            double ap=*(*(a+(i-1)/n)+((i-1)%n));



        if(!(i%n)){
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
        //fprintf(results,"%.13lf\t%.13lf\n",tp+dt,A);
        if(fabs(dt)<1e-12)break;
        }
        int imax=i;
        //fclose(results);
        //results=fopen("Vector.txt","w");



    //Initial condition for Tb = Tr
    for (i = 0; i < imax; i++)
        if(1/(*(*(a+i/n)+(i%n)))-1>=zeq)break;
    int ieq=i;
    *(*(Tb+ieq/n)+(ieq%n))=*(*(Trad+ieq/n)+(ieq%n));

    /*for (i = 0; i < imax/n+1; i++){
        free(*(a+i));
        free(*(t+i));
        }
    free(a);
    free(t);*/

   /* double afunc(double T)
{   int j;
    for (j = 0; j < imax; j++)
           if((*(*(t+j/n)+(j%n)))<=T)break;
        double ao=*(*(a+j/n)+(j%n));
        double to=*(*(t+j/n)+(j%n));
        double dt=T-to;
    return rk4(expansion, dt, to, ao);
}

double DTbar(double t,double Tbar)
{
    double A=afunc(t);
    double B=-2*(expansion(t,A))/A*Tbar;//sure
    double C=-DC*Trado*Trado*Trado*Trado*(Trado/A-Tbar)/A/A/A/A*XE(TI/Tbar);
    return B+C;
}*/

    for (i = ieq; i > 0; i--){
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
        //*(*(Tb+(i-1)/n)+((i-1)%n)) = rk4(DTbar, dt, tr, TB);
        }
    for (i = ieq; i < imax-1; i++){
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
        //*(*(Tb+(i+1)/n)+((i+1)%n)) = rk4(DTbar, dt, tr, TB);
        }
    for (i = 0; i < imax; i++)
    fprintf(results,"%.13lf\t%.13lf\t%.13lf\t%.13lf\t%d\n",\
    *(*(t+i/n)+(i%n)), *(*(a+i/n)+(i%n)),*(*(Trad+i/n)+(i%n)),*(*(Tb+i/n)+(i%n)),i);



        FILE* XEE=fopen("XE.txt","w");
        for (i = 0; i < 1e4; i++)
        fprintf(XEE,"%.13lf\t%.13lf\n",i*0.01,XE(i*0.01));
        fclose(XEE);


    fclose(results);

    return 0;
}
