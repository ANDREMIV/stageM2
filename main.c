#include <stdio.h>
#include <stdlib.h>
#include "mathutils.h"

#include "parsing.h"
#include "simplecosmomodels.h"
#include <string.h>
///t is the cosmic time t multiplied by Hubble constant Ho
///a is the scale factor

///DANGERS WITH THE GLOBAL VARS SEE SIMPLECOSMOMODELS.H BE CAREFULL !!!



struct Vector_Data GD;

int expansion_calc2(double **t, double **a, double **Trad,double **Tb,struct datfile* p);

int main(void)
{
    //deLCDM();
    //desitter();
    //einsteindesitter();

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

    GD.a=a;GD.t=t;GD.Tb=Tb;GD.Tr=Tr;
    struct datfile h2_h,hd_h_rot,h2_h_p_rot;
    DATinit("h2-h.dat",&h2_h);
    DATinit("hd-h-rot.dat",&hd_h_rot);
    DATinit("h2-h-H+-rot.dat",&h2_h_p_rot);
    //int imax = expansion_calc(t,a,Tr,Tb);
    int imax2 = expansion_calc2(t,a,Tr,Tb,(&h2_h_p_rot));
    //XEplot();

    /*FILE* RAD=fopen("radoutComp.txt","w"); //forces new empty file for later cat >> commands
    fclose(RAD);*/


    /*DummyRadexOut(2);
    double densities[8]={0};
    densities[0]=1e-4;
    densities[1]=1e-8;
    double TB=100; double LC,GH;
        radexinp(0.1,TB,&h2_h_p_rot,densities);
    char com[200]= {0};
    char mm[64]= {0};
        strcpy(mm,h2_h_p_rot.DAT_file_name);
        strtok(mm,".");
        strcat(mm,".inp");
        memset(com, 0, sizeof(com));
        strcat(com, "C:\\Radex\\bin\\radex.exe < ");
        strcat(com,mm);
        system(com);
        double *levels=radexout(&h2_h_p_rot);
        Cooling_heating(&LC,&GH,&h2_h_p_rot,levels,6.3e-7*densities[0],densities,TB);
        printf("\n%le\t%le\t%le",LC,GH,(LC-GH)/(6.3e-7*densities[0])*1e7);
        free(levels);*/

   //LVL(&h2_h,10,3402,10);
   // square_cooling_power(&h2_h,0,3e3,15,"squareH2HCP.txt");

    //square_cooling_power(&h2_h_p_rot,1,100, 1000,30,"squareH2H+CP.txt");
    //square_cooling_power(&h2_h_p_rot,0,100, 1000,10,"squareH2H'CP.txt");
    //square_cooling_power(&h2_h,0,100,1e4,50,"squareH2HCP.txt");

    DATfree(&h2_h);
    DATfree(&hd_h_rot);
    DATfree(&h2_h_p_rot);


    return 0;
}
