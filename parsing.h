#include "simplecosmomodels.h"

#define COL_NB_MAX 10
#define goto_next_line(file) while(fgetc(file)!='\n');
#define skip_n_lines(file,n) {int i=n;for(;i--;)while(fgetc(file)!='\n');}

void radexinp(int Col_nb, char* outinp, double Trad, double Tbar, char* mol, char** Cols_name,double* densities, int rot);
double* radexout(int Col_nb,int npop,int nlines,int (*QNs)[2]);
void ROPT(char* DAT_file_name);
