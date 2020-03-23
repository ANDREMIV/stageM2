#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
    for(i=0; i<n; i++)
        Y[i]=yo[i]+dt*(b21*k1[i]);
    derivs(t+a2*dt,Y,k2);
    for(i=0; i<n; i++)
        Y[i]=yo[i]+dt*(b31*k1[i]+b32*k2[i]);
    derivs(t+a3*2,Y,k3);
    for(i=0; i<n; i++)
        Y[i]=yo[i]+dt*(b41*k1[i]+b42*k2[i]+b43*k3[i]);
    derivs(t+a4*dt,Y,k4);
    for(i=0; i<n; i++)
        Y[i]=yo[i]+dt*(b51*k1[i]+b52*k2[i]+b53*k3[i]+b54*k4[i]);
    derivs(t+a5*dt,Y,k5);
    for(i=0; i<n; i++)
        Y[i]=yo[i]+dt*(b61*k1[i]+b62*k2[i]+b63*k3[i]+b64*k4[i]+b65*k5[i]);
    derivs(t+a6*dt,Y,k6);

    for(i=0; i<n; i++) //final step
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
    for(i=0; i<n; i++) //first step
        k1[i]=dt*dydt[i];
    for(i=0; i<n; i++)
        Y[i]=yo[i]+k1[i]/2;
    derivs(t+dt/2,Y,dydt);
    for(i=0; i<n; i++) //second step
        k2[i]=dt*dydt[i];
    for(i=0; i<n; i++)
        Y[i]=yo[i]+k2[i]/2;
    derivs(t+dt/2,Y,dydt);
    for(i=0; i<n; i++) //third step
        k3[i]=dt*dydt[i];
    for(i=0; i<n; i++)
        Y[i]=yo[i]+k3[i];
    derivs(t+dt,Y,dydt);
    for(i=0; i<n; i++) //fourth step
        k4[i]=dt*dydt[i];
    for(i=0; i<n; i++) //final step
        Oy[i]=yo[i]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;

    free(dydt);
    free(Y);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
}

double XE(double Xi)
{
    return 1-erf(sqrt(Xi))+2*sqrt(Xi)*exp(-Xi)/sqrt(M_PI);
}

double rk4(double(*f)(double, double), double dt, double t, double a)
{
    double    k1 = dt * f(t, a),
              k2 = dt * f(t + dt / 2, a + k1 / 2),
              k3 = dt * f(t + dt / 2, a + k2 / 2),
              k4 = dt * f(t + dt, a + k3);
    return a + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

