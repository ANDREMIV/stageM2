#include <stdio.h>
#include <stdlib.h>
#include "mathutils.h"
#include "simplecosmomodels.h"
#include "parsing.h"
///t is the cosmic time t multiplied by Hubble constant Ho
///a is the scale factor

///DANGERS WITH THE GLOBAL VARS SEE SIMPLECOSMOMODELS.H BE CAREFULL !!!



struct Vector_Data GD;


int main(void)
{
    deLCDM();
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

    GD.a=a;GD.t=t;GD.Tb=Tb;GD.Tr=Tr;


    int imax = expansion_calc(t,a,Tr,Tb);

    XEplot();

    /*FILE* RAD=fopen("radoutComp.txt","w"); //forces new empty file for later cat >> commands
    fclose(RAD);*/

    ROPT("h2-h.dat");
    //ROPT("hd-h-rot.dat");
    //ROPT("h2-h-H+-rot.dat");

    return 0;
}
