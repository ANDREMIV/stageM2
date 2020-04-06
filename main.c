#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mathutils.h"
#include "simplecosmomodels.h"
///t is the cosmic time t multiplied by Hubble constant Ho
///a is the scale factor
///f(a,t)= da/dt
///DANGERS WITH THE GLOBAL VARS SEE SIMPLECOSMOMODELS.H BE CAREFULL !!!

void radexinpHp(double Trad, double Tbar, double densityH, double densityp);
void radexinpH(double Trad, double Tbar, double densityH);
double* radexoutHp(int n);
double* radexoutH(int n);

int main(void)
{
    desitter();
    einsteindesitter();

    double **a,**Tr, **Tb, **t;
    int n = MEMORYBLOC, m=NBMEMBLOC;


    t = (double **)malloc(sizeof(double*) * (m));
    *t = (double *)malloc(sizeof(double) * (n));
    a = (double **)malloc(sizeof(double*) * (m));
    *a = (double *)malloc(sizeof(double) * (n));
    Tr= (double **)malloc(sizeof(double) * (m));
    *Tr= (double *)malloc(sizeof(double) * (n));
    Tb= (double **)malloc(sizeof(double) * (m));
    *Tb= (double *)malloc(sizeof(double) * (n));

    double xe1e_4(){return 1e-4;}

    int imax = expansion_calc(t,a,Tr,Tb);

    XEplot();

    FILE* RAD=fopen("radoutComp.txt","w"); //forces new empty file for later cat >> commands
    fclose(RAD);



    /*{
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
            radexinp2(Trado*(1+z), Trado*D(a,ieq2)*(1+z)*(1+z), H);
            system("C:\\Radex\\bin\\radex.exe < h22o.inp");
            system("cat < radex.out >> radoutComp.txt");




            fprintf(LEVELS,"Trad: %le\tTbar: %le\tNH: %le\t",\
                    Trado*(1+z), Trado*D(a,ieq2)*(1+z)*(1+z), H);

            double* levels;
            levels=radexout(nlines);

            for (i = 0; i < nlines; i++)
                fprintf(LEVELS,"%le\t%le\n",*(levels+2*i+0%2),*(levels+2*i+1%2));
            double rortho=0,rpara=0; rpara+=*(levels+2*0+1%2); rortho+=*(levels+2*1+1%2);
            for (i = 0; i < nlines; i++)
                if(i%2)rortho+=*(levels+2*i+0%2); else rpara+=*(levels+2*i+0%2);
            fprintf(ROP,"%d\t%le\t%le\t%le\t\n",z,rortho/rpara,BoltzmanH2ROP(Trado*(1+z)),BoltzmanH2ROP(Trado*D(a,ieq2)*(1+z)*(1+z)));

            free(levels);
        }

    }*/
    void ROPTHp()
    {
        system("cp radexHH+.out radex.out");
         FILE* LEVELS=fopen("levelsT.txt","w");
    FILE* ROP=fopen("ropT.txt","w"); //ratio ORTHO PARA
    fprintf(LEVELS,"POP UP\t\t\tPOP LOW\n");
    int nlines=7;
        int z;
        int i;

        for(z=10; z<zeq; z++)
        {
            for(i=0; 1/D(a,i)-1>z; i++);int ic=i;

            double H= 0.05*(1+z)*(1+z)*(1+z)*H0*H0*3/8/M_PI/G/MH*1e-6;//*1e-6; in cm^3
            radexinpHp(D(Tr,ic), D(Tb,ic), H, H*4e-4);
            system("C:\\Radex\\bin\\radex.exe < h2o.inp");
            system("cat < radex.out >> radoutTComp.txt");




            fprintf(LEVELS,"Trad: %le\tTbar: %le\tNH: %le\tNp: %le\n",\
                    D(Tr,ic), D(Tb,ic), H, H*4e-4);

            double* levels;
            levels=radexoutHp(nlines);

            for (i = 0; i < nlines; i++)
                fprintf(LEVELS,"%le\t%le\n",*(levels+2*i+0%2),*(levels+2*i+1%2));
            double rortho=0,rpara=0; rpara+=*(levels+2*0+1%2); rortho+=*(levels+2*1+1%2);
            for (i = 0; i < nlines; i++)
                if(i%2)rortho+=*(levels+2*i+0%2); else rpara+=*(levels+2*i+0%2);
            fprintf(ROP,"%d\t%le\t%le\t%le\t\n",z,rortho/rpara,BoltzmanH2ROP(D(Tr,ic)),BoltzmanH2ROP(D(Tb,ic)));

            free(levels);
        }
        fclose(LEVELS);
    fclose(ROP);

    }




    void ROPHp()
    {
        system("cp radexHH+.out radex.out");
         FILE* LEVELS=fopen("levels.txt","w");
    FILE* ROP=fopen("rop.txt","w"); //ratio ORTHO PARA
    fprintf(LEVELS,"POP UP\t\t\tPOP LOW\n");
    int nlines=7;
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
            radexinpHp(Trado*(1+z), Trado*D(a,ieq2)*(1+z)*(1+z), H, H*4e-4);
            system("C:\\Radex\\bin\\radex.exe < h2o.inp");
            system("cat < radex.out >> radoutComp.txt");




            fprintf(LEVELS,"Trad: %le\tTbar: %le\tNH: %le\tNp: %le\n",\
                    Trado*(1+z), Trado*D(a,ieq2)*(1+z)*(1+z), H, H*4e-4);

            double* levels;
            levels=radexoutHp(nlines);

            for (i = 0; i < nlines; i++)
                fprintf(LEVELS,"%le\t%le\n",*(levels+2*i+0%2),*(levels+2*i+1%2));
            double rortho=0,rpara=0; rpara+=*(levels+2*0+1%2); rortho+=*(levels+2*1+1%2);
            for (i = 0; i < nlines; i++)
                if(i%2)rortho+=*(levels+2*i+0%2); else rpara+=*(levels+2*i+0%2);
            fprintf(ROP,"%d\t%le\t%le\t%le\t\n",z,rortho/rpara,BoltzmanH2ROP(Trado*(1+z)),BoltzmanH2ROP(Trado*D(a,ieq2)*(1+z)*(1+z)));

            free(levels);
        }
        fclose(LEVELS);
    fclose(ROP);

    }
    void ROPH()
    {
        system("cp radexH.out radex.out");
         FILE* LEVELS=fopen("levels2.txt","w");
    FILE* ROP=fopen("rop2.txt","w"); //ratio ORTHO PARA
    fprintf(LEVELS,"POP UP\t\t\tPOP LOW\n");
    int nlines=7;
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
            radexinpH(Trado*(1+z), Trado*D(a,ieq2)*(1+z)*(1+z), H);
            system("C:\\Radex\\bin\\radex.exe < h22o.inp");
            system("cat < radex.out >> radoutComp2.txt");




            fprintf(LEVELS,"Trad: %le\tTbar: %le\tNH: %le\n",\
                    Trado*(1+z), Trado*D(a,ieq2)*(1+z)*(1+z), H);

            double* levels;
            levels=radexoutH(nlines);

            for (i = 0; i < nlines; i++)
                fprintf(LEVELS,"%le\t%le\n",*(levels+2*i+0%2),*(levels+2*i+1%2));
            double rortho=0,rpara=0; rpara+=*(levels+2*0+1%2); rortho+=*(levels+2*1+1%2);
            for (i = 0; i < nlines; i++)
                if(i%2)rortho+=*(levels+2*i+0%2); else rpara+=*(levels+2*i+0%2);
            fprintf(ROP,"%d\t%le\t%le\t%le\t\n",z,rortho/rpara,BoltzmanH2ROP(Trado*(1+z)),BoltzmanH2ROP(Trado*D(a,ieq2)*(1+z)*(1+z)));

            free(levels);
        }
        fclose(LEVELS);
    fclose(ROP);

    }
    ROPTHp();
    //ROPH();
    //ROPHp();




    return 0;
}
