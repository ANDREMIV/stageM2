#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mathutils.h"
#include "simplecosmomodels.h"
#define D(pointer,index) (*(*(pointer+(index)/n)+((index)%n)))
///t is the cosmic time t multiplied by Hubble constant Ho
///a is the scale factor
///f(a,t)= da/dt
///DANGERS WITH THE GLOBAL VARS SEE SIMPLECOSMOMODELS.H BE CAREFULL !!!

void radexinp(double Trad, double Tbar, double densityH, double densityp);
double* radexout(int n);


FILE* results;


int main(void)
{
    radexinp(123,456,789,909);
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



    FILE* XEE=fopen("XE.txt","w");
    for (i = 0; i < 1e4; i++)
        fprintf(XEE,"%.13lf\t%.13lf\n",i*0.01,XE(i*0.01));
    fclose(XEE);

    FILE* RAD=fopen("radoutComp.txt","w");
    fclose(RAD);
    FILE* LEVELS=fopen("levels.txt","w");
    FILE* ROP=fopen("rop.txt","w"); //ratio ORTHO PARA
    fprintf(LEVELS,"POP UP\t\t\tPOP LOW\n");
    int nlines=7;
    {
        int z;
        int i;
        int ieq2;
        int zeq2=150;
        for(i=0; i<imax; i++)
            if(1/(*(*(a+i/n)+(i%n)))-1>=zeq2)
                break;
        ieq2=i;

        for(z=0; z<zeq2; z++)
        {

            double H= 0.05*(1+z)*(1+z)*(1+z)*H0*H0*3/8/M_PI/G/MH*1e-6;//*1e-6; in cm^3
            radexinp(Trado*(1+z), Trado*D(a,ieq2)*(1+z)*(1+z), H, H*4e-4);
            system("C:\\Radex\\bin\\radex.exe < h2o.inp");
            system("cat < radex.out >> radoutComp.txt");


            double* levels;
            levels=radexout(nlines);

            fprintf(LEVELS,"Trad: %le\tTbar: %le\tNH: %le\tNp: %le\n",\
                    Trado*(1+z), Trado*D(a,ieq2)*(1+z)*(1+z), H, H*1e-4);
            for (i = 0; i < nlines; i++)
                fprintf(LEVELS,"%le\t%le\n",*(levels+2*i+0%2),*(levels+2*i+1%2));
            double rortho=0,rpara=0; rpara+=*(levels+2*0+1%2); rortho+=*(levels+2*1+1%2);
            for (i = 0; i < nlines; i++)
                if(i%2)rortho+=*(levels+2*i+0%2); else rpara+=*(levels+2*i+0%2);
            fprintf(ROP,"%d\t%le\t%le\t%le\t\n",z,rortho/rpara,BoltzmanH2ROP(Trado*(1+z)),BoltzmanH2ROP(Trado*D(a,ieq2)*(1+z)*(1+z)));

            free(levels);
        }

    }


    fclose(results);
    fclose(LEVELS);
    fclose(ROP);

    return 0;
}
