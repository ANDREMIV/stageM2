#define MEMORYBLOC 1e5
#define NBMEMBLOC 100
#define D(pointer,index) (*(*(pointer+(index)/n)+((index)%n)))

void rk42(void (*derivs)(double, double*, double*), int n\
          ,double dt, double t, double*yo, double *Oy);



void rk62(void (*derivs)(double, double*, double*), int n\
          ,double dt, double t, double*yo, double *Oy);


void derivs(double t, double*y, double*dydt);

double rk4(double(*f)(double, double), double dt, double t, double a);

double XE(double Xi);

void XEplot();

int expansion_calc(double **t, double **a, double **Trad,double **Tb);
