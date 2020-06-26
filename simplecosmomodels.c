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

void n_lsoda_terminate(void);
typedef void    (*_lsoda_f) (double, double *, double *, void *);
void lsoda(_lsoda_f f, int neq, double *y, double *t, double tout, int itol, double *rtol, double *atol,
		   int itask, int *istate, int iopt, int jt,
		   int iwork1, int iwork2, int iwork5, int iwork6, int iwork7, int iwork8, int iwork9,
		   double rwork1, double rwork5, double rwork6, double rwork7, void *_data);
#define IVRS 3
#define NBDENS 2
#define TR 2
#define TB 1
int expansion_calc2(double **t, double **a, double **Tr,double **Tb,struct datfile* p)
{
    FILE* ROP=fopen("ROP.txt","w");
//    Oro=AR*Trado*Trado*Trado*Trado*8*M_PI*G/3/C/C/H0/H0;//=9.877e-7;
    #ifdef NOSTATEQ
    FILE* results=fopen("expansion2.txt","w");
    #else
    FILE* results=fopen("STATEQlevels.txt","w");
    #endif
    int i, n = MEMORYBLOC, m=NBMEMBLOC;
    int N=IVRS+NBDENS+p->npop;
    double *yo=(double *)malloc(sizeof(double) * (N)); //initial condition for RK
    double *Oy=(double *)malloc(sizeof(double) * (N)); //output of RK


        double *atol=(double *)malloc(sizeof(double) * (N+1));
    double *rtol=(double *)malloc(sizeof(double) * (N+1));
    double *yy=(double *)malloc(sizeof(double) * (N+1));
    rtol[0] = 0.0;
	atol[0] = 0.0;
    double          rwork1, rwork5, rwork6, rwork7;
	double          tt, tout;
	int             iwork1, iwork2, iwork5, iwork6, iwork7, iwork8, iwork9;
	int             neq = N;
	int             itol, itask, istate, iopt, jt;

	iwork1 = iwork2 = iwork5 = iwork6 = iwork7 = iwork8 = iwork9 = 0;
	rwork1 = rwork5 = rwork6 = rwork7 = 0.0;
	{
	    int i;
        for(i=0;i<IVRS+NBDENS;i++)
        atol[i+1]=0;
        for(i=IVRS+NBDENS;i<N;i++)
        atol[i+1]=1e-5;
        for(i=0;i<IVRS+NBDENS;i++)
        rtol[i+1]=1e-5;
        for(i=IVRS+NBDENS;i<N;i++)
        rtol[i+1]=1e-3;
	}
		itol = 2;
	itask = 1;
	istate = 1;
	iopt = 0;
	jt = 2;



    double t0 = 0; //t0 = 0 when we are in zeq
    const double sdt=1e-10;
    double dt=sdt;
    **a=1.0/(zeq+1);
    **Tr=Trado*(zeq+1);
    **Tb=**Tr;
    **t=t0;
    double DT=1e-3;
    fprintf(results,"t'\ta(t')\tTrad(t')\tTb(t')\tz(t')\ti\n------------\n");
        i=0;
    #ifdef FROMNEQ

    double TVEC[]={1.38e-003, 1.1329e-002, 1.0006e+002, 2.4057e+002, 0, 0, 1.72e-002, 1.51e-001, 8.07e-002, 3.12e-001, 1.23e-001, 2.58e-001, 1.18e-003, 1.50e-004, 5.98e-003, 3.20e-006, 2.68e-005, 1.34e-005, 4.74e-005, 4.34e-002, 1.78e-005, 5.10e-005, 2.55e-003, 1.95e-005, 2.58e-003, 4.37e-005, 1.92e-005, 1.45e-007, 1.23e-006, 1.16e-004, 6.19e-007, 2.21e-006, 3.62e-005, 8.08e-007, 1.67e-004, 2.30e-006, 1.22e-005, 7.70e-007, 8.00e-006, 1.42e-005, 1.70e-006, 5.66e-007, 1.43e-008, 1.24e-007, 3.40e-006, 6.52e-008, 1.45e-005, 2.50e-007, 1.07e-006, 9.73e-008, 1.84e-006, 3.11e-007, 5.34e-007, 2.89e-007, 1.08e-007, 3.91e-007, 5.81e-007, 3.11e-007, 1.19e-006, 1.01e-007};
    {int i;for(i=0;i<N;i++)
    yo[i]=TVEC[i+1];
    }
    double tc=TVEC[0]; //current time
    #else

    yo[0]=D(a,i);
        yo[1]=D(Tb,i);
        yo[2]=D(Tr,i);
        double tc=D(t,i); //current time
    #endif

        yo[3]=0.05/yo[0]/yo[0]/yo[0]*H0*H0*3/8/M_PI/G/MH*1e-6;//densities must be in cm^3
        yo[4]=1e-4*yo[3];


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
        for (i = 0; i < p->npop; i++)
        fprintf(results,"v=%d___j=%d\t",p->QNs[i][0],p->QNs[i][1]);
        fprintf(results,"\n");
        printf("\n");
        #ifdef NOSTATEQ

        #else
        {
        }
        fprintf(results,"%.7le\t%.7le\t%.7le\t%.7le\t%.7le\n",\
                D(t,i), D(a,i),D(Tr,i),D(Tb,i),1/D(a,i)-1);
        #endif
        #ifndef FROMNEQ
            {int i;
            for(i=0;i<N-IVRS-NBDENS;i++)
                yo[i+IVRS+NBDENS]=levels[i];
                }
            #endif // FROMNEQ
        free(levels);
    char *cond=(char *)malloc(sizeof(char) * (N));


    for (i = 0; D(a,i)<1/10.0 && i<n*m-1; i++)
    {
        if(!((i+1)%n))
        {
            (*(t+(i+1)/n))=(double *)malloc(sizeof(double) * (n));
            (*(a+(i+1)/n))=(double *)malloc(sizeof(double) * (n));
            (*(Tr+(i+1)/n))=(double *)malloc(sizeof(double) * (n));
            (*(Tb+(i+1)/n))=(double *)malloc(sizeof(double) * (n));
        }

        int j; int sub=10;
        for(j=0;j<sub;j++){
        //rk42(deriv_pop_net,N,dt,tc,yo,Oy,(void*)(p));
        {
    {
	    int i;
	    for(i=0;i<N;i++)
        yy[i+1]=yo[i];
		}

	tt = tc;
	tout = tc+dt;

		lsoda(deriv_pop_net, neq, yy, &tt, tout, itol, rtol, atol, itask, &istate, iopt, jt,
		      iwork1, iwork2, iwork5, iwork6, iwork7, iwork8, iwork9,
		      rwork1, rwork5, rwork6, rwork7, (void*)(p));
		if (istate <= 0) {
			printf("error istate = %d\n", istate);
			exit(0);
		}
		{
	    int i;
	    for(i=0;i<N;i++)
        Oy[i]=yy[i+1];
		}

        }
        tc=tc+dt;


        if(j==sub-1);///time step management
        ;;;
        {
            int i;
            for(i=0;i<IVRS+NBDENS;i++)
            {
                double TVAR=fabs(Oy[i])/fabs(yo[i]);
                if(TVAR>1+10*DT||TVAR<1-9*DT)
                {cond[i]=0;break;}//lessen the time and redo calc
                else if(TVAR<1+DT&&TVAR>1-0.9*DT)
                cond[i]=2;//augment time
            else
                cond[i]=1;//keep time
            }

            for(i=IVRS+NBDENS;i<N;i++)
            {
                if(Oy[i]<0||Oy[i]>1){
                            cond[i]=2;break;}
                            else cond[i]=2;
            }

            for(i=0;i<N;i++)
            if(!cond[i]){//lessen time and redo calc
                    {int i;for(i=0;i<N;i++)
                    Oy[i]=yo[i];}
                    tc-=dt;
                    j--;
                    dt/=2;//if(dt<sdt)dt*=2;
                    break;
            }

            if(i==N)
            for(i=0;i<N;i++)
            if(cond[i]!=2)break;

            if(i==N)
                dt*=2;
        }

        {
            int i;
            for(i=0;i<N;i++)
                yo[i]=Oy[i];
        }

        }



        D(a,i+1)=Oy[0];
        D(Tb,i+1)=Oy[1];
        D(Tr,i+1)=Oy[2];
        D(t,i+1)=tc;//+dt;

        #ifdef NOSTATEQ
        fprintf(results,"%.2le\t%.4le\t%.4le\t%.4le\t%.4le\t§§§\t",\
                tc, yo[0],yo[2],yo[1],1/yo[0]-1);
        #endif // NOSTATEQ
        {
        int i;
        double rfg=0,zo=0,zp=0;
        for (i = 0; i < p->npop; i++)
            rfg+=yo[i+IVRS+NBDENS];
        for (i = 0; i < p->npop; i++)
            if(!(p->QNs[i][1]%2))zp+=yo[i+IVRS+NBDENS];
            else zo+=yo[i+IVRS+NBDENS];
        double ROP=zo/zp;
        double LC =0, GH=0;
        double nH2=6.3e-7*yo[IVRS];
        Cooling_heating(&LC,&GH,p,&(yo[IVRS+NBDENS]),nH2,&(yo[IVRS]),yo[TB]);

        fprintf(results,"%.2le\t%.4le\t%.3le\t%.4le\t",tc,rfg,ROP,GH-LC);
        for (i = 0; i < p->npop; i++)
            fprintf(results,"%.2le\t",yo[i+IVRS+NBDENS]);
        fprintf(results,"\n");}

    }
    int imax=i;

    fclose(results);
    fclose(LEVELS);
    fclose(ROP);
    free(yo);
    free(Oy);
    free(cond);
    free(atol);
    free(rtol);
    free(yy);
        n_lsoda_terminate();
    return imax;
}
#undef IVRS
#undef NBDENS
#undef TR
#undef TB

