#include <stdio.h>
#include <stdlib.h>
#include "mathutils.h"
#include <math.h>
#include "parsing.h"
#include <string.h>
///IN SI
const double C=2.99792458e8;
const double KB=1.380649e-23;
const double MH=1.67265e-27;
const double Me=9.1093837015e-31;
const double TS=6.6524587158e-29;
const double AR=7.5657e-16;
const double hPl=6.62606876e-34;
const double G=6.674e-11;
const double H0=6.736e1*1e3/1e6/3.08567758149e16;
const double Ol=0.6847,Omo=0.3153;
const double Trado=2.7255;
const double TI=1.579e5; //13.6eV/kb in Kelvin
const double zeq=3402;//+-26 plank 2018
const double Oro=9.28e-5;

double reciprocal_coef(double Tkin, double Eul, double gl, double gu, double coef_ul)
{
    return coef_ul*gu/gl*exp(-Eul/KB/Tkin);
}

double expansion(double t, double a)
{
    return  sqrt(Ol*a*a+Omo/a+Oro/a/a);
}



void derivs2(double t, double*y, double*dydt)//n=3
//y0=a
//y1=Tbar
//y2=Trad
{
    double DC=8.0/3*TS*AR/Me/C/H0; //Compton coupling constant
    dydt[0]=expansion(0,y[0]);
    dydt[1]=-2*dydt[0]/y[0]*y[1]\
            +DC*y[2]*y[2]*y[2]*y[2]*(y[2]-y[1])*1e-4;
    dydt[2]=-y[2]*dydt[0]/y[0];
}

void desitter()
{
   // Oro=AR*Trado*Trado*Trado*Trado*8*M_PI*G/3/C/C/H0/H0;//=9.877e-7;
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
void deLCDM()
{
   // Oro=AR*Trado*Trado*Trado*Trado*8*M_PI*G/3/C/C/H0/H0;//=9.877e-7;
    /*double expansion3(double t, double a)
    {
        return  a;
    }*/
    FILE* desit=fopen("expansionLCDM.txt","w");
    double a=1;
    double t0 = 0, t1 = -1, dt=-1e-5,t;
    for(t=t0; t>t1; t+=dt)
    {
        a=rk4(expansion, dt, t, a);
        fprintf(desit,"%.13lf\t%.13lf\n",t+dt,a);
    }
    fclose(desit);
}
void einsteindesitter()
{
//    Oro=AR*Trado*Trado*Trado*Trado*8*M_PI*G/3/C/C/H0/H0;//=9.877e-7;
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



int expansion_calc(double **t, double **a, double **Tr,double **Tb)
{
//    Oro=AR*Trado*Trado*Trado*Trado*8*M_PI*G/3/C/C/H0/H0;//=9.877e-7;
    FILE* results=fopen("expansion.txt","w");
    int i, n = MEMORYBLOC, m=NBMEMBLOC;
    double t0 = 0; //t0 = 0 when we are in zeq
    double dt=1e-15;
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
        for(j=0;j<100;j++){

        rk62(derivs2,3,dt,tc,yo,Oy,NULL);

        yo[0]=Oy[0];
        yo[1]=Oy[1];
        yo[2]=Oy[2];
        tc=tc+dt;}



        D(a,i+1)=Oy[0];
        //if(Oy[1]>Oy[2])Oy[1]=Oy[2];//Delete positive surcompensation
        D(Tb,i+1)=Oy[1];
        D(Tr,i+1)=Oy[2];
        D(t,i+1)=tc;//+dt;

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

#define IVRS 3
#define NBDENS 2
int expansion_calc2(double **t, double **a, double **Tr,double **Tb,struct datfile* p)
{
//    Oro=AR*Trado*Trado*Trado*Trado*8*M_PI*G/3/C/C/H0/H0;//=9.877e-7;
    FILE* results=fopen("expansion.txt","w");
    int i, n = MEMORYBLOC, m=NBMEMBLOC;
    int N=IVRS+NBDENS+p->npop;
    double *yo=(double *)malloc(sizeof(double) * (N)); //initial condition for RK
    double *Oy=(double *)malloc(sizeof(double) * (N)); //output of RK
    double t0 = 0; //t0 = 0 when we are in zeq
    double dt=1e-15;
    **a=1.0/(zeq+1);
    **Tr=Trado*(zeq+1);
    **Tb=**Tr;
    **t=t0;
    double DT=1e-3;
    fprintf(results,"t'\ta(t')\tTrad(t')\tTb(t')\tz(t')\ti\n------------\n");

    i=0;

    yo[0]=D(a,i);
        yo[1]=D(Tb,i);
        yo[2]=D(Tr,i);
        yo[3]=0.05/yo[0]/yo[0]/yo[0]*H0*H0*3/8/M_PI/G/MH*1e-6;//densities must be in cm^3
        yo[4]=1e-4*yo[3];
        double tc=D(t,i); //current time

    DummyRadexOut(p->Col_nb);
    char com[64]= {0};
    char mm[64]= {0};
    memset(com, 0, sizeof(com));
    strcpy(mm,p->DAT_file_name);
    char* unptr=strtok(mm,".");
    sprintf(com, "Lavels%s.txt",unptr);
    FILE* LEVELS=fopen(com,"w");
    radexinp(yo[2],yo[1],p,&(yo[IVRS]));
        strcpy(mm,p->DAT_file_name);
        unptr=strtok(mm,".");
        strcat(mm,".inp");
        memset(com, 0, sizeof(com));
        strcat(com, "C:\\Radex\\bin\\radex.exe < ");
        strcat(com,mm);
        system(com);

        double* levels;
        levels=radexout(p);
        #ifdef NOSTATEQ
        {
            int i;
            for(i=0;i<N-IVRS-NBDENS;i++)
                yo[i+IVRS+NBDENS]=levels[i];
            /*for(i=0;i<N-3-2;i++)
                printf("%lf\t",levels[i]);*/
        }
        #else
        {
            int i;
            for(i=0;i<N-IVRS-NBDENS;i++)
                yo[i+IVRS+NBDENS]=1.0/(N-IVRS-NBDENS);

        }
        #endif
        free(levels);
    char *cond=(char *)malloc(sizeof(char) * (N));
    for (i = 0; D(a,i)<1 && i<n*m-1; i++)
    {
        if(!((i+1)%n))
        {
            (*(t+(i+1)/n))=(double *)malloc(sizeof(double) * (n));
            (*(a+(i+1)/n))=(double *)malloc(sizeof(double) * (n));
            (*(Tr+(i+1)/n))=(double *)malloc(sizeof(double) * (n));
            (*(Tb+(i+1)/n))=(double *)malloc(sizeof(double) * (n));
        }

        int j; int sub=2;
        for(j=0;j<sub;j++){
        rk42(deriv_pop_net,N,dt,tc,yo,Oy,(void*)(p));
        tc=tc+dt;
        if(j==sub-1)
        {
            int i;
            for(i=0;i<N;i++)
            {
                double TVAR=fabs(Oy[i])/fabs(yo[i]);
                if(TVAR>1+10*DT||TVAR<1-9*DT)
                cond[i]=0;
                else if(TVAR<1+DT&&TVAR>1-0.9*DT)
                cond[i]=2;
            else
                cond[i]=1;
            }
            for(i=0;i<N;i++)
            if(!cond[i]){dt/=2;break;}
            else if(cond[i]!=2)break;

            if(i==N)dt*=1.5;
        }

        {
            int i;
            for(i=0;i<N;i++)
                yo[i]=Oy[i];
        }

        }



        D(a,i+1)=Oy[0];
        //if(Oy[1]>Oy[2])Oy[1]=Oy[2];//Delete positive surcompensation
        D(Tb,i+1)=Oy[1];
        D(Tr,i+1)=Oy[2];
        D(t,i+1)=tc;//+dt;

        /*//Variable step
        double TVAR=D(a,i+1)/D(a,i);
        if(TVAR>1+10*DT)
            dt/=2;
        else if(TVAR<1+DT)
            dt*=2;*/

        fprintf(results,"%.7le\t%.7le\t%.7le\t%.7le\t%.7le\t§§§\t",\
                D(t,i), D(a,i),D(Tr,i),D(Tb,i),1/D(a,i)-1);
        {
        int i;
        double rfg=0;
        for (i = 0; i < p->npop; i++)
            rfg+=yo[i+IVRS+NBDENS];
        fprintf(results,"%le\t$$$\t",rfg);
        for (i = 0; i < p->npop; i++)
            fprintf(results,"%le\t",yo[i+IVRS+NBDENS]);
        fprintf(results,"\n");}
    }
    int imax=i;

    fclose(results);
    fclose(LEVELS);
    free(yo);
    free(Oy);
    free(cond);
    return imax;
}
#undef IVRS
#undef NBDENS

