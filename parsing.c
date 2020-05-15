#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parsing.h"
#include "mathutils.h"


struct col_coefs
{
    int nb_temps;
    int nb_col_trans;
    float* coefs;
    float* temps;
    int (*ij)[2];
};

void ROPT(char* DAT_file_name)
{
    int Col_nb=0;
    char * mol;
    char *Cols[COL_NB_MAX];
    double densities[COL_NB_MAX];
///find nb of collisionners and chars of collisionners
    char DAT_file_name2[64]= {0}; //carefull do not use this array after, it holds name of COLs
    strcpy(DAT_file_name2,DAT_file_name);
    mol = strtok (DAT_file_name2,"-.");
    for (Col_nb=0; 1; Col_nb++)
        if((Cols[Col_nb] = strtok (NULL, "-.")) ==NULL)
            break;
    Col_nb--;//remove .dat col
    int g=1;
    int z;//g hold condition to say rot was here
    int i;
    for (i=0; i<Col_nb; i++)
        if(!strcmp(Cols[i],"hp"))
            strcpy(Cols[i],"H+");
    if(!strncmp(Cols[Col_nb-1],"rot",3))
        Col_nb--;
    else
        g--; //rot not a collisionner

    ///Opens conserrponding level/rop files
    char com[200]= {0};
    char mm[64]= {0};
    sprintf(com, "cp %dcol.out radex.out",Col_nb); //put dummy radex output file in case radex do not produce output
    system(com);
    memset(com, 0, sizeof(com));
    strcpy(mm,DAT_file_name);
    char* unptr=strtok(mm,".");
    sprintf(com, "levels%s.txt",unptr);
    FILE* LEVELS=fopen(com,"w");
    memset(com, 0, sizeof(com));
    sprintf(com, "ROP%s.txt",unptr);
    FILE* ROP=fopen(com,"w"); //ratio ORTHO PARA


    ///copy DAT file here for examination
    memset(com, 0, sizeof(com));
    strcat(com,"cp C:\\Radex\\data\\");
    strcat(com,DAT_file_name);
    strcat(com," ");
    strcat(com,DAT_file_name);
    system(com);
    FILE* FDAT=fopen(DAT_file_name,"r");
    skip_n_lines(FDAT, 5);
    int npop;
    fscanf(FDAT,"%d",&npop);
    skip_n_lines(FDAT, 2);
    int (*QNs)[2];//Quantum numbers, supports 2 quantum numbers
    float *Wg;
    float *Ener;
    Wg=(float*)calloc(npop,sizeof(float));
    Ener=(float*)calloc(npop,sizeof(float));
    QNs=(int(*)[2])calloc(npop,sizeof(int[2]));
    {
        int i;
        for(i=0; i<npop; i++)
        {
            fscanf(FDAT,"%*s%f%f%d%*[_]%d",&Ener[i],&Wg[i],&(QNs[i][0]),&(QNs[i][1]));
        }
    }

    goto_next_line(FDAT);
    goto_next_line(FDAT);
    int nlines=0;
    fscanf(FDAT,"%d",&nlines);
    skip_n_lines(FDAT, nlines+3);
    fscanf(FDAT,"%d",&Col_nb);
    goto_next_line(FDAT);
    struct col_coefs *Col_coefs= (struct col_coefs *) calloc(Col_nb,sizeof(struct col_coefs));
    {
        int i;
        for(i=0; i<Col_nb; i++)
        {
            int ctrans, ctemps;
            skip_n_lines(FDAT, 3);
            fscanf(FDAT,"%d",&ctrans);
            (Col_coefs[i].nb_col_trans)=ctrans;
            skip_n_lines(FDAT, 2);
            fscanf(FDAT,"%d",&ctemps);
            (Col_coefs[i].nb_temps)=ctemps;
            Col_coefs[i].coefs=(float*)calloc(ctrans*ctemps,sizeof(float));
            Col_coefs[i].temps=(float*)calloc(ctemps,sizeof(float));
            Col_coefs[i].ij=(int(*)[2])calloc(ctrans,sizeof(int[2]));
            skip_n_lines(FDAT, 2);
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


    fclose(FDAT);
    memset(com, 0, sizeof(com));
    strcat(com,"del ");
    strcat(com,DAT_file_name);
    system(com);

    ///write levels info
    for (i = 0; i < npop; i++)
        fprintf(LEVELS,"v=%d___j=%d\t",QNs[i][0],QNs[i][1]);
    fprintf(LEVELS,"\n");
    fprintf(LEVELS,"z\tTrad [K]\tTbar [K]\tHeating-Cooling [J/s/cm^3]\tdot Tbar Heating-Cooling [K/s]\t");
    for (i = 0; i < Col_nb; i++)
        fprintf(LEVELS,"%s\t",Cols[i]);
    fprintf(LEVELS,"\n");


    for(z=10; z<1000; z+=5)
    {
        char outinp[64]= {0};
        for(i=0; 1/D(GD.a,i)-1>z; i++);
        int ic=i;

        double H= 0.05*(1+z)*(1+z)*(1+z)*H0*H0*3/8/M_PI/G/MH*1e-6;//*1e-6; in cm^3
        double nH2=6.3e-7*H;

        {
            int i;
            for(i=0; i<Col_nb; i++)
                densities[i]= i<1 ? H : H*4e-4;
        }
        double customTbar=D(GD.Tb,ic);
        radexinp(Col_nb,outinp,D(GD.Tr,ic),customTbar,mol,Cols,densities,g);
        memset(com, 0, sizeof(com));
        strcat(com, "C:\\Radex\\bin\\radex.exe < ");
        strcat(com,outinp);
        system(com);

        double* levels;
        levels=radexout(Col_nb,npop,nlines,QNs);

        ///Cooling
        double LC =0;
        double GH=0;
        {
            int i,j,k;
            for(k=0; k<Col_nb; k++)
            {

                int ctrans, ctemps;
                ctrans=(Col_coefs[k].nb_col_trans);
                ctemps=(Col_coefs[k].nb_temps);
                for(i=0; i<ctrans; i++)
                {
                    //for(j=1;j<npop;j++)
                    int l=Col_coefs[k].ij[i][1]-1;
                    int u=Col_coefs[k].ij[i][0]-1;
                    double Eul=100*C*hPl*(Ener[u]-Ener[l]);
                    double Cul;
                    for(j=0; j<ctemps; j++)
                        if(D(GD.Tb,ic)<=Col_coefs[k].temps[j])
                            break;
                    if(j==ctemps)
                    {
                        printf("\nerror Clu(Tb) outside of interpolation range");
                        system("pause");
                        exit(1);
                    }
                    ///linear fit
                    if(j!=0)
                    Cul=Col_coefs[k].coefs[i*ctemps+j]+(D(GD.Tb,ic)-Col_coefs[k].temps[j])*(Col_coefs[k].coefs[i*ctemps+j-1]\
                            -Col_coefs[k].coefs[i*ctemps+j])/(Col_coefs[k].temps[j-1]-Col_coefs[k].temps[j]);
                    else Cul=D(GD.Tb,ic)/Col_coefs[k].temps[j]*Col_coefs[k].coefs[i*ctemps+j];
                    GH+=nH2*densities[k]*levels[u]*Cul*Eul;
                    LC+=nH2*densities[k]*levels[l]*reciprocal_coef(customTbar,Eul,Wg[l],Wg[u],Cul)*Eul;
                }
            }
        }

        ///printing

        fprintf(LEVELS,"%d\t%le\t%le\t%le\t%le\t",z,D(GD.Tr,ic), D(GD.Tb,ic),GH-LC,(GH-LC)*1e6*2/KB/3/H);
        for (i = 0; i < Col_nb; i++)
            fprintf(LEVELS,"%le\t",densities[i]);
        fprintf(LEVELS,"\n");

        for (i = 0; i < npop; i++)
            fprintf(LEVELS,"%le\t",levels[i]);
        fprintf(LEVELS,"\n");



        /* double rortho=0,rpara=0;
         for (i = 0; i < npop; i++)
             if(i%2)rortho+=levels[i]; else rpara+=levels[i];
         fprintf(ROP,"%d\t%le\t\n",z,rortho/rpara);*/

        free(levels);
    }
    fclose(LEVELS);
    fclose(ROP);
    free(Wg);
    free(Ener);
    free(QNs);

    {
        int i;
        for(i=0; i<Col_nb; i++)
        {
            free(Col_coefs[i].temps);
            free(Col_coefs[i].coefs);
            free(Col_coefs[i].ij);
        }
        free(Col_coefs);
    }

}

///Creates a input file for radex
void radexinp(int Col_nb, char* outinp, double Trad, double Tbar, char* mol, char** Cols_name,double* densities, int rot)
{
    char com[64]= {0};
    strcat(com,mol);
    int i;
    for(i=0; i<Col_nb; i++)
    {
        strcat(com,"-");
        strcat(com,Cols_name[i]);
    }
    if(rot)
        strcat(com,"-rot");
    strcat(com,".inp");
    FILE* radexOIN=fopen(com,"w");//new input file
    strcpy(outinp,com);
    memset(com, 0, sizeof(com));

    strtok(outinp,".");
    strcat(outinp,".dat");

    fprintf(radexOIN,"%s\nradex.out\n100 400000\n%le\n",outinp,Tbar);
    strtok(outinp,".");
    strcat(outinp,".inp");
    if(!strcmp(Cols_name[Col_nb-1],"rot"))
        Col_nb--;
    fprintf(radexOIN,"%d\n",Col_nb);
    int j;
    for(j=0; j<Col_nb; j++)
        fprintf(radexOIN,"%s\n%le\n",Cols_name[j],densities[j]);
    fprintf(radexOIN,"%le\n1e15\n1.0\n0\n",Trad);

    fclose(radexOIN);
}

specific_search_int(int (*p)[2],int n,int a, int b)
{
    int j;
    for(j=0; j<n; j++)
        if((p[j][0]==a)&&(p[j][1]==b))
            break;
    return j;
}

///Open radex output file in order to retrieve population levels
double* radexout(int Col_nb,int npop,int nlines,int (*QNs)[2])
{

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
