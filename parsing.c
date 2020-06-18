#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parsing.h"
#include "mathutils.h"

void DummyRadexOut(int nb_col)
{
    char com[128]= {0};
    sprintf(com, "cp %dcol.out radex.out",nb_col); //put dummy radex output file in case radex do not produce output
    system(com);
}

///Retrieves useful infos from DAT file
void DATinit(char* DAT_file_name, struct datfile* d)
{
    strcpy(d->DAT_file_name,DAT_file_name);
        char com[200]= {0};
        //copy DAT file here for examination
    memset(com, 0, sizeof(com));
    strcat(com,"cp C:\\Radex\\data\\");
    strcat(com,DAT_file_name);
    strcat(com," ");
    strcat(com,DAT_file_name);
    system(com);
    FILE* FDAT=fopen(DAT_file_name,"r");
    goto_next_line(FDAT);
    fscanf(FDAT,"%s",d->mol_name);
    skip_n_lines(FDAT, 4);

    fscanf(FDAT,"%d",&d->npop);
    int npop=d->npop;
    skip_n_lines(FDAT, 2);
    int (*QNs)[2];//Quantum numbers, supports 2 quantum numbers
    float *Wg;
    double *Ener;
    d->Wg=(float*)calloc(npop,sizeof(float));
    d->Ener=(double*)calloc(npop,sizeof(double));
    d->QNs=(int(*)[2])calloc(npop,sizeof(int[2]));
    Wg=d->Wg;
    Ener=d->Ener;
    QNs=d->QNs;
    {
        int i;
        for(i=0; i<npop; i++)
        {
            fscanf(FDAT,"%*s%lf%f%d%*[_]%d",&Ener[i],&Wg[i],&(QNs[i][0]),&(QNs[i][1]));
            Ener[i]*=100*C*hPl;
        }
    }

    skip_n_lines(FDAT, 2);

    int nlines=0;
    fscanf(FDAT,"%d",&d->nlines);
    nlines=d->nlines;
    skip_n_lines(FDAT, 2);
    d->A.nbtrans=nlines;
    d->A.A=(double*)calloc(nlines,sizeof(double));
    d->A.Babs=(double*)calloc(nlines,sizeof(double));
    d->A.Bse=(double*)calloc(nlines,sizeof(double));
    d->A.ij=(int(*)[2])calloc(nlines,sizeof(int[2]));
    {
        int i;
        for(i=0; i<nlines; i++)
        {
            double nu; //freq in GHz
            fscanf(FDAT,"%*s%d%d%le%lf",&(d->A.ij[i][0]),&(d->A.ij[i][1]), &(d->A.A[i]), &nu);///
            d->A.Bse[i]=d->A.A[i]*C*C/2.0/hPl/nu/nu/nu*1e-27; //freq in GHz
            int l=d->A.ij[i][1]-1;
            int u=d->A.ij[i][0]-1;
            d->A.Babs[i]=d->A.Bse[i]*d->Wg[u]/d->Wg[l];
            goto_next_line(FDAT);

        }
    }
    goto_next_line(FDAT);



    int Col_nb;
    fscanf(FDAT,"%d",&d->Col_nb);
    Col_nb=d->Col_nb;
    struct col_coefs *Col_coefs= (struct col_coefs *) calloc(Col_nb,sizeof(struct col_coefs));

    skip_n_lines(FDAT, 1);

    {
        int i;
        for(i=0; i<Col_nb; i++)
        {
            int ctrans, ctemps;
            skip_n_lines(FDAT, 1);
            fscanf(FDAT,"%*d%*s%*s%s",Col_coefs[i].Col_name);
            skip_n_lines(FDAT, 2);
            fscanf(FDAT,"%d",&ctrans);
            (Col_coefs[i].nb_col_trans)=ctrans;
            skip_n_lines(FDAT, 2);
            fscanf(FDAT,"%d",&ctemps);
            (Col_coefs[i].nb_temps)=ctemps;
            Col_coefs[i].coefs=(double*)calloc(ctrans*ctemps,sizeof(double));
            Col_coefs[i].temps=(double*)calloc(ctemps,sizeof(double));
            Col_coefs[i].ij=(int(*)[2])calloc(ctrans,sizeof(int[2]));
            skip_n_lines(FDAT, 2);
            d->C=Col_coefs;
            {
                int j;
                for(j=0; j<ctemps; j++)
                    fscanf(FDAT,"%lf",&(Col_coefs[i].temps[j]));

            }

            skip_n_lines(FDAT, 2);
            {
                int j;
                for(j=0; j<ctrans; j++)
                {
                    fscanf(FDAT,"%*s%d%d",&(Col_coefs[i].ij[j][0]),&(Col_coefs[i].ij[j][1]));
                    //fscanf(FDAT,"%*c");

                    {
                        int k;
                        for(k=0; k<ctemps; k++)
                            fscanf(FDAT,"%le",&(Col_coefs[i].coefs[j*ctemps+k]));
                        //printf("%e", Col_coefs[i].coefs[j*ctemps+k]);
                    }
                    goto_next_line(FDAT);
                }
            }

        }
    }

//delete copied file

    fclose(FDAT);
    memset(com, 0, sizeof(com));
    strcat(com,"del ");
    strcat(com,DAT_file_name);
    system(com);
}

void DATfree(struct datfile* d)
{
    free(d->Wg);
    free(d->Ener);
    free(d->QNs);
    free(d->A.A);
    free(d->A.Babs);
    free(d->A.Bse);
    free(d->A.ij);
    {
        int i;
        for(i=0; i<d->Col_nb; i++)
        {
            free(d->C[i].temps);
            free(d->C[i].coefs);
            free(d->C[i].ij);
        }
        free(d->C);
    }
}

void square_cooling_power(struct datfile* d, int which, int Tmin, int Tmax, int STEP, char* outname)
{

    FILE* OUT=fopen(outname,"w");
    double TB;
    fprintf(OUT,"TB");
    for(TB=Tmin;TB<=Tmax;TB+=STEP)
    {
    fprintf(OUT,"\n%d|\t",(int)(TB));

        DummyRadexOut(d->Col_nb);
     double LC,GH;double ds[8]={1e-0,1e-4,1,1,1,1,1,1};
        radexinp(0.1,TB,d,ds);
    char com[200]= {0};
    char mm[64]= {0};
        strcpy(mm,d->DAT_file_name);
        strtok(mm,".");
        strcat(mm,".inp");
        memset(com, 0, sizeof(com));
        strcat(com, "C:\\Radex\\bin\\radex.exe < ");
        strcat(com,mm);
        system(com);
        double *levels=radexout(d);
        Cooling_heating_power_per_collisionner(&LC,&GH,which,d,levels,TB);
        fprintf(OUT,"%.1le\t",(LC-GH));
        free(levels);

    }
    fclose(OUT);
}

void cube_cooling_power(struct datfile* d, int which,int Tmin, int Tmax, int STEP, char* outname)
{
    FILE* OUT=fopen(outname,"w");
    double TB,TR;
    fprintf(OUT,"TB\\TR\t");
    for(TR=Tmin;TR<=Tmax;TR+=STEP)
    fprintf(OUT,"%d\t\t",(int)(TR));
    for(TB=Tmin;TB<=Tmax;TB+=STEP)
    {
    fprintf(OUT,"\n%d|\t",(int)(TB));
    for(TR=Tmin;TR<=Tmax;TR+=STEP)
    {

        DummyRadexOut(d->Col_nb);
     double LC,GH;double ds[8]={1,1,1,1,1,1,1,1};
        radexinp(TR,TB,d,ds);
    char com[200]= {0};
    char mm[64]= {0};
        strcpy(mm,d->DAT_file_name);
        strtok(mm,".");
        strcat(mm,".inp");
        memset(com, 0, sizeof(com));
        strcat(com, "C:\\Radex\\bin\\radex.exe < ");
        strcat(com,mm);
        system(com);
        double *levels=radexout(d);
        Cooling_heating_power_per_collisionner(&LC,&GH,which,d,levels,TB);
        fprintf(OUT,"%.1le\t",(LC-GH));
        free(levels);



    }

    }
    fclose(OUT);
}

double Cul_interpol(struct col_coefs* Col_coefs, int k, int i, double Tb)
{//i is transition index, k is collisionner index
    double Cul;
    int j; int ctemps=(Col_coefs[k].nb_temps);
    for(j=0; j<ctemps; j++)
                        if(Tb<=Col_coefs[k].temps[j])
                            break;
                    if(j==ctemps)
                    {
                        /*printf("\nerror Clu(Tb) outside of interpolation range");
                        system("pause");
                        exit(1);*/
                        Cul=Col_coefs[k].coefs[i*ctemps+j-1];//maxout approximation
                    }else{
                    ///linear fit
                    if(j!=0)
                    Cul=Col_coefs[k].coefs[i*ctemps+j]+(Tb-Col_coefs[k].temps[j])*(Col_coefs[k].coefs[i*ctemps+j-1]\
                            -Col_coefs[k].coefs[i*ctemps+j])/(Col_coefs[k].temps[j-1]-Col_coefs[k].temps[j]);
                    else Cul=Tb/Col_coefs[k].temps[j]*Col_coefs[k].coefs[i*ctemps+j];} //minout approximation
    return Cul;

}


void Cooling_heating_power_per_collisionner(double *cooling, double *heating, int which, struct datfile* d, double* levels, double Tb)
//levels are the populations of the collisionned in the context of datfile d
//J.cm^3/s
{
    ///Cooling
        double LC =0;
        double GH=0;
        struct col_coefs* Col_coefs=d->C;
        {
            int i,k;
            k=which;
            {

                int ctrans;
                ctrans=(Col_coefs[k].nb_col_trans);
                for(i=0; i<ctrans; i++)
                {
                    //for(j=1;j<npop;j++)
                    int l=Col_coefs[k].ij[i][1]-1;
                    int u=Col_coefs[k].ij[i][0]-1;
                    double Eul=(d->Ener[u]-d->Ener[l]);
                    double Cul=Cul_interpol(Col_coefs,k,i,Tb);

                    GH+=levels[u]*Cul*Eul;
                    LC+=levels[l]*reciprocal_coef(Tb,Eul,d->Wg[l],d->Wg[u],Cul)*Eul;
                }
            }
        }
    *cooling=LC;
    *heating=GH;
}

void Cooling_heating(double *cooling, double *heating, struct datfile* d, double* levels, double n, double* densities, double Tb)
//n is density of the one being collisionned
//densities are the densities of the collisionners
//levels are the populations of the collisionned in the context of datfile d
//J/cm^3/s
{
    ///Cooling
        double LC =0;
        double GH=0;
        struct col_coefs* Col_coefs=d->C;
        {
            int i,k;
            for(k=0; k<d->Col_nb; k++)
            {

                int ctrans;
                ctrans=(Col_coefs[k].nb_col_trans);
                for(i=0; i<ctrans; i++)
                {
                    //for(j=1;j<npop;j++)
                    int l=Col_coefs[k].ij[i][1]-1;
                    int u=Col_coefs[k].ij[i][0]-1;
                    double Eul=(d->Ener[u]-d->Ener[l]);
                    double Cul=Cul_interpol(Col_coefs,k,i,Tb);

                    GH+=n*densities[k]*levels[u]*Cul*Eul;
                    LC+=n*densities[k]*levels[l]*reciprocal_coef(Tb,Eul,d->Wg[l],d->Wg[u],Cul)*Eul;
                }
            }
        }
    *cooling=LC;
    *heating=GH;
}

///Creates a input file for radex in the context of an .dat file
void radexinp(double Trad, double Tbar,struct datfile* d,double* densities)
{
    char com[64]= {0};
    char mm[64]= {0};
    memset(com, 0, sizeof(com));
    strcpy(mm,d->DAT_file_name);
    strtok(mm,".");
    strcpy(com,mm);
    strcat(com,".inp");
    FILE* radexOIN=fopen(com,"w");//new input file



    fprintf(radexOIN,"%s\nradex.out\n100 400000\n%le\n",d->DAT_file_name,Tbar);
    fprintf(radexOIN,"%d\n",d->Col_nb);
    int j;
    for(j=0; j<d->Col_nb; j++)
        fprintf(radexOIN,"%s\n%le\n",d->C[j].Col_name,densities[j]);
    fprintf(radexOIN,"%le\n1e15\n1.0\n0\n",Trad);

    fclose(radexOIN);
}

int specific_search_int(int (*p)[2],int n,int u, int l)
{
    int j;
    for(j=0; j<n; j++)
        if((p[j][0]==u)&&(p[j][1]==l))
            break;
    if(j==n)return -1;
    return j;
}

///Open radex output file in order to retrieve population levels in the context of a dat file
double* radexout(struct datfile* d)
{
    int Col_nb=d->Col_nb;int npop=d->npop;int nlines=d->nlines; int (*QNs)[2]=d->QNs;
    FILE* Fradexout=fopen("radex.out","r");

    int i;
    skip_n_lines(Fradexout,10+Col_nb);
    double *levels;
    levels=calloc(npop,sizeof(double));

    int UP[2]= {0};
    int DN[2]= {0};
    double levelup;
    double leveldn;
    for(i=0; i<nlines; i++)
    {
        fscanf(Fradexout,"%d%*[_]%d%*[ -]%d%*[_]%d%*s%*s%*s%*s%*s%*s%le%le",&UP[0],&UP[1],&DN[0],&DN[1],&levelup,&leveldn); //discard previous info to obtain only the levels
        int j;
        //printf("%d",UP[1]);system("pause");
        j=specific_search_int(QNs,npop,UP[0],UP[1]);
        levels[j]=levelup;
        j=specific_search_int(QNs,npop,DN[0],DN[1]);
        levels[j]=leveldn;
        goto_next_line(Fradexout);
    }

    fclose(Fradexout);
    return levels;
}

///print levels in the context of a dat file
void LVL(struct datfile* d)
{
    DummyRadexOut(d->Col_nb);
    double densities[8]={0};
    int i,z;
    ///Opens conserrponding level/rop files
    char com[200]= {0};
    char mm[64]= {0};
    memset(com, 0, sizeof(com));
    strcpy(mm,d->DAT_file_name);
    char* unptr=strtok(mm,".");
    sprintf(com, "levels%s.txt",unptr);
    FILE* LEVELS=fopen(com,"w");
    ///write levels info
    for (i = 0; i < d->npop; i++)
        fprintf(LEVELS,"v=%d___j=%d\t",d->QNs[i][0],d->QNs[i][1]);
    fprintf(LEVELS,"\n");
    fprintf(LEVELS,"z\tTrad [K]\tTbar [K]\tHeating-Cooling [J/s/cm^3]\tdot Tbar Heating-Cooling [K/s]\tintegrated Heating cooling over zdiff [K]\t");
    for (i = 0; i < d->Col_nb; i++)
        fprintf(LEVELS,"%s\t",d->C[i].Col_name);
    fprintf(LEVELS,"\n");

    int zdif=5;
    for(z=10; z<200; z+=zdif)
    {
        double LC =0;
        double GH=0;
        for(i=0; 1/D(GD.a,i)-1>z; i++);
        int ic=i;
        for(i=0; 1/D(GD.a,i)-1>z+zdif; i++);
        double dtd= D(GD.t,ic)-D(GD.t,i);

        double H= 0.05*(1+z)*(1+z)*(1+z)*H0*H0*3/8/M_PI/G/MH*1e-6;//*1e-6; in cm^3
        double nH2=6.3e-7*H;

        {
            int i;
            for(i=0; i<d->Col_nb; i++)
                densities[i]= i<1 ? H : H*4e-4;
        }
        double customTbar=D(GD.Tb,ic);
        radexinp(D(GD.Tr,ic),customTbar,d,densities);
        strcpy(mm,d->DAT_file_name);
        unptr=strtok(mm,".");
        strcat(mm,".inp");
        memset(com, 0, sizeof(com));
        strcat(com, "C:\\Radex\\bin\\radex.exe < ");
        strcat(com,mm);
        system(com);

        double* levels;
        levels=radexout(d);

        Cooling_heating(&LC, &GH,d,levels,nH2,densities,D(GD.Tb,ic));

        ///printing
        double KL=(GH-LC)*2/KB/3/H; //Kelvin loss K/s
        double IKL= dtd*KL/H0;

        fprintf(LEVELS,"%d\t%le\t%le\t%le\t%le\t%le\t",z,D(GD.Tr,ic), D(GD.Tb,ic),GH-LC,KL,IKL);
        for (i = 0; i < d->Col_nb; i++)
            fprintf(LEVELS,"%le\t",densities[i]);
        fprintf(LEVELS,"\n");

        for (i = 0; i < d->npop; i++)
            fprintf(LEVELS,"%le\t",levels[i]);
        fprintf(LEVELS,"\n");

        free(levels);
    }
    fclose(LEVELS);


}

double Plank_nu(double nu,double T)
{
    return 2*hPl*nu*nu*nu/C/C/(exp(hPl*nu/KB/T)-1);
}

#define TR 2
#define TB 1

#define IVRS 3
#define NBDENS 2
void deriv_pop_net(double t, double*y, double*dydt,void* params)
{
    //y0=a
    //y1=Tb
    //y2=Tr
    //y3=p0
    //y4=p1
    //...
    //yn=pn
    #ifdef NOSTATEQ
    double DC=8.0/3*TS*AR/Me/C/H0; //Compton coupling constant
    dydt[0]=expansion(0,y[0]);
    dydt[1]=-2*dydt[0]/y[0]*y[1]\
            +DC*y[2]*y[2]*y[2]*y[2]*(y[2]-y[1])*1e-4;
    dydt[2]=-y[2]*dydt[0]/y[0];
    dydt[3]=-3*dydt[0]/y[0]*y[3];
    dydt[4]=1e-4*dydt[3];
    #else
    dydt[0]=0;
    dydt[1]=0;
    dydt[2]=0;
    dydt[3]=0;
    dydt[4]=0;
    #endif
    struct datfile* p = (struct datfile*) params;

    int n=p->npop;
    int i,j,k;
    double r2=0;
    for(i=0;i<n;i++)
    {
        double r=0;
        for(j=0;j<n;j++)
        {
            int m= i > j ? -1 : 1;
            int u= i > j ? i : j;
            int l= i > j ? j : i;



            int i2=specific_search_int(p->A.ij,p->A.nbtrans,u+1,l+1);
        double Eul=(p->Ener[u]-p->Ener[l]);
            if(i2!=-1)
            {
            r+=m* p->A.A[i2]*y[u+IVRS+NBDENS]; //in-out spontaneous emissions
            r-=m*p->A.Babs[i2]*y[l+IVRS+NBDENS]*Plank_nu(Eul/hPl,y[TR]); //in-out absorptions
            r+=m*p->A.Bse[i2]*y[u+IVRS+NBDENS]*Plank_nu(Eul/hPl,y[TR]); //in-out stimulated emissions
            }
            for(k=0;k<p->Col_nb;k++)
            {
            int i5=specific_search_int(p->C[k].ij,p->C[k].nb_col_trans,u+1,l+1);
            if(i5!=-1)
            {
                double Cul=Cul_interpol(p->C,k,i5,y[TB]);
                r-=m*reciprocal_coef(y[TB],Eul,p->Wg[l],p->Wg[u],Cul)*y[l+IVRS+NBDENS]*y[IVRS+k]; //in-out collisional excitations
                r+=m*Cul*y[u+IVRS+NBDENS]*y[IVRS+k]; //in-out collisional desexcitations
            }
            }

        }
    r2+=r;
    dydt[i+IVRS+NBDENS]=r/H0;
    }

}
#undef TR
#undef TB
#undef IVRS
#undef NBDENS
