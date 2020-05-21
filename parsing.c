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
    float *Ener;
    d->Wg=(float*)calloc(npop,sizeof(float));
    d->Ener=(float*)calloc(npop,sizeof(float));
    d->QNs=(int(*)[2])calloc(npop,sizeof(int[2]));
    Wg=d->Wg;
    Ener=d->Ener;
    QNs=d->QNs;
    {
        int i;
        for(i=0; i<npop; i++)
        {
            fscanf(FDAT,"%*s%f%f%d%*[_]%d",&Ener[i],&Wg[i],&(QNs[i][0]),&(QNs[i][1]));
        }
    }

    skip_n_lines(FDAT, 2);

    int nlines=0;
    fscanf(FDAT,"%d",&d->nlines);
    nlines=d->nlines;
    skip_n_lines(FDAT, nlines+3);
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
            Col_coefs[i].coefs=(float*)calloc(ctrans*ctemps,sizeof(float));
            Col_coefs[i].temps=(float*)calloc(ctemps,sizeof(float));
            Col_coefs[i].ij=(int(*)[2])calloc(ctrans,sizeof(int[2]));
            skip_n_lines(FDAT, 2);
            d->C=Col_coefs;
            {
                int j;
                for(j=0; j<ctemps; j++)
                    fscanf(FDAT,"%f",&(Col_coefs[i].temps[j]));

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
                            fscanf(FDAT,"%e",&(Col_coefs[i].coefs[j*ctemps+k]));
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
            int i,j,k;
            for(k=0; k<d->Col_nb; k++)
            {

                int ctrans, ctemps;
                ctrans=(Col_coefs[k].nb_col_trans);
                ctemps=(Col_coefs[k].nb_temps);
                for(i=0; i<ctrans; i++)
                {
                    //for(j=1;j<npop;j++)
                    int l=Col_coefs[k].ij[i][1]-1;
                    int u=Col_coefs[k].ij[i][0]-1;
                    double Eul=100*C*hPl*(d->Ener[u]-d->Ener[l]);
                    double Cul;
                    for(j=0; j<ctemps; j++)
                        if(Tb<=Col_coefs[k].temps[j])
                            break;
                    if(j==ctemps)
                    {
                        printf("\nerror Clu(Tb) outside of interpolation range");
                        system("pause");
                        exit(1);
                    }
                    ///linear fit
                    if(j!=0)
                    Cul=Col_coefs[k].coefs[i*ctemps+j]+(Tb-Col_coefs[k].temps[j])*(Col_coefs[k].coefs[i*ctemps+j-1]\
                            -Col_coefs[k].coefs[i*ctemps+j])/(Col_coefs[k].temps[j-1]-Col_coefs[k].temps[j]);
                    else Cul=Tb/Col_coefs[k].temps[j]*Col_coefs[k].coefs[i*ctemps+j];
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

int specific_search_int(int (*p)[2],int n,int a, int b)
{
    int j;
    for(j=0; j<n; j++)
        if((p[j][0]==a)&&(p[j][1]==b))
            break;
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
    fprintf(LEVELS,"z\tTrad [K]\tTbar [K]\tHeating-Cooling [J/s/cm^3]\tdot Tbar Heating-Cooling [K/s]\t");
    for (i = 0; i < d->Col_nb; i++)
        fprintf(LEVELS,"%s\t",d->C[i].Col_name);
    fprintf(LEVELS,"\n");


    for(z=10; z<200; z+=5)
    {
        double LC =0;
        double GH=0;
        for(i=0; 1/D(GD.a,i)-1>z; i++);
        int ic=i;

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

        fprintf(LEVELS,"%d\t%le\t%le\t%le\t%le\t",z,D(GD.Tr,ic), D(GD.Tb,ic),GH-LC,(GH-LC)*1e6*2/KB/3/H);
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
