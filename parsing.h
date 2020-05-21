#include "simplecosmomodels.h"

#define COL_NB_MAX 10
#define goto_next_line(file) while(fgetc(file)!='\n');
#define skip_n_lines(file,n) {int i=n;for(;i--;)while(fgetc(file)!='\n');}

struct datfile
{
    char mol_name[8];
    int Col_nb;
    struct col_coefs* C;
    int npop;
    int nlines;
    int (*QNs)[2];
    float *Ener;
    float *Wg;
    char DAT_file_name[64];
};

struct col_coefs
{
    int nb_temps;
    int nb_col_trans;
    float* coefs;
    float* temps;
    int (*ij)[2];
    char Col_name[8];
};

void radexinp(double Trad, double Tbar,struct datfile* d,double* densities);
double* radexout( struct datfile* d);
void DATfree(struct datfile* d);
void DATinit(char* DAT_file_name, struct datfile* d);
void DummyRadexOut(int nb_col);
void Cooling_heating(double *cooling, double *heating, struct datfile* d, double* levels, double n, double* densities, double Tb);
void LVL(struct datfile* d);
