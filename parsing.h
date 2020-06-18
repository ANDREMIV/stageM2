#include "simplecosmomodels.h"

#define COL_NB_MAX 10
#define goto_next_line(file) while(fgetc(file)!='\n');
#define skip_n_lines(file,n) {int i=n;for(;i--;)while(fgetc(file)!='\n');}



struct Einstein_col
{
    int nbtrans;
    double* A; //spontaneous emission
    int (*ij)[2];
    double* Bse; // stimulated emission
    double* Babs; // absorption
};

struct col_coefs
{
    int nb_temps;
    int nb_col_trans;
    double* coefs;
    double* temps;
    int (*ij)[2];
    char Col_name[8];
};

struct datfile
{
    char mol_name[8];
    int Col_nb;
    struct col_coefs* C;
    struct Einstein_col A;
    int npop;
    int nlines;
    int (*QNs)[2];
    double *Ener; //in J, if in cm-1, multiply by 100*C*hPl to have in J
    float *Wg;
    char DAT_file_name[64];
};

struct ensemble_data
{
    struct datfile* d;
    double* col_densities;//nb_cols is the number of col_densities
};

void radexinp(double Trad, double Tbar,struct datfile* d,double* densities);
double* radexout( struct datfile* d);
void DATfree(struct datfile* d);
void DATinit(char* DAT_file_name, struct datfile* d);
void DummyRadexOut(int nb_col);
void Cooling_heating(double *cooling, double *heating, struct datfile* d, double* levels, double n, double* densities, double Tb);
void Cooling_heating_power_per_collisionner(double *cooling, double *heating, int which, struct datfile* d, double* levels, double Tb);
void LVL(struct datfile* d);
void cube_cooling_power(struct datfile* d, int which, int Tmin, int Tmax, int STEP, char* outname);
void square_cooling_power(struct datfile* d, int which, int Tmin, int Tmax, int STEP, char* outname);
void deriv_pop_net(double t, double*y, double*dydt,void* params);
