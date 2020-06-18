#define MEMORYBLOC ((int)1e5)
#define NBMEMBLOC ((int)100)
#define D(pointer,index) (*(*(pointer+(index)/MEMORYBLOC)+((index)%MEMORYBLOC)))

void rk42(void (*derivs)(double, double*, double*, void*), int n\
          ,double dt, double t, double*yo, double *Oy, void* params);



void rk62(void (*derivs)(double, double*, double*, void*), int n\
          ,double dt, double t, double*yo, double *Oy, void* params);



double rk4(double(*f)(double, double), double dt, double t, double a);

double XE(double Xi);

void XEplot();

struct Vector_Data
{
    double** t,**a,**Tr,**Tb;
};

extern struct Vector_Data GD;
