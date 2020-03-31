#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mathutils.h"
#include "simplecosmomodels.h"

void radexinp(double Trad, double Tbar, double densityH, double densityp);
double* radexout(int n);

void ROP1(int imax)
    {
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
            radexinp(Trado*(1+z), Trado*D(a,ieq2)*(1+z)*(1+z), H, H*4e-4);
            system("C:\\Radex\\bin\\radex.exe < h2o.inp");
            system("cat < radex.out >> radoutComp.txt");




            fprintf(LEVELS,"Trad: %le\tTbar: %le\tNH: %le\tNp: %le\n",\
                    Trado*(1+z), Trado*D(a,ieq2)*(1+z)*(1+z), H, H*4e-4);

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
        fclose(LEVELS);
    fclose(ROP);

    }
