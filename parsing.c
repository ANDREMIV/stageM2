#include <stdio.h>
#include <stdlib.h>

void radexinpHp(double Trad, double Tbar, double densityH, double densityp)
{
    FILE* radexin=fopen("h2.inp","r");
    FILE* radexout=fopen("h2o.inp","w");

    char s[200]= {0};
    int k;
    int i;
    for (k=0,i=0; k<3&&((s[i]=fgetc(radexin))!=EOF); i++)
        if(s[i]=='\n')
            k++;
    fprintf(radexout,"%s%le\n",s,Tbar);
    while(fgetc(radexin)!='\n');

    for(i=0; i<200; i++)
        s[i]=0;
    for (k=3,i=0; k<5&&((s[i]=fgetc(radexin))!=EOF); i++)
        if(s[i]=='\n')
            k++;
    fprintf(radexout,"%s%le\n",s,densityH);
    while(fgetc(radexin)!='\n');

    for(i=0; i<200; i++)
        s[i]=0;
    for (k=5,i=0; k<6&&((s[i]=fgetc(radexin))!=EOF); i++)
        if(s[i]=='\n')
            k++;
    fprintf(radexout,"%s%le\n",s,densityp);
    while(fgetc(radexin)!='\n');

    fprintf(radexout,"%le\n",Trad);
    while(fgetc(radexin)!='\n');
    for(i=0; i<200; i++)
        s[i]=0;
    for (i=0; ((s[i]=fgetc(radexin))!=EOF); i++)
        fputc(s[i],radexout);

    fclose(radexin);
    fclose(radexout);
}

void radexinpH(double Trad, double Tbar, double densityH)
{
    FILE* radexin=fopen("h22.inp","r");
    FILE* radexout=fopen("h22o.inp","w");

    char s[200]= {0};
    int k;
    int i;
    for (k=0,i=0; k<3&&((s[i]=fgetc(radexin))!=EOF); i++)
        if(s[i]=='\n')
            k++;
    fprintf(radexout,"%s%.10lf\n",s,Tbar);
    while(fgetc(radexin)!='\n');

    for(i=0; i<200; i++)
        s[i]=0;
    for (k=3,i=0; k<5&&((s[i]=fgetc(radexin))!=EOF); i++)
        if(s[i]=='\n')
            k++;
    fprintf(radexout,"%s%.10lf\n",s,densityH);
    while(fgetc(radexin)!='\n');

    for(i=0; i<200; i++)
        s[i]=0;
    /*for (k=5,i=0; k<6&&((s[i]=fgetc(radexin))!=EOF); i++)
        if(s[i]=='\n')
            k++;
    fprintf(radexout,"%s%.10lf\n",s,densityp);
    while(fgetc(radexin)!='\n');*/

    fprintf(radexout,"%.10lf\n",Trad);
    while(fgetc(radexin)!='\n');
    for(i=0; i<200; i++)
        s[i]=0;
    for (i=0; ((s[i]=fgetc(radexin))!=EOF); i++)
        fputc(s[i],radexout);

    fclose(radexin);
    fclose(radexout);
}

double* radexoutH(int n)
{

    FILE* radexout=fopen("radex.out","r");

    char s[2000]= {0};
    int k;
    int i;
    for (k=0,i=0; k<11&&((s[i]=fgetc(radexout))!=EOF); i++)
        if(s[i]=='\n')
            k++;
    double *levels;
    levels=malloc(sizeof(double)*2*n);
    for(i=0; i<2*n; i++)
        *(levels+i)=0;

    for(i=0; i<n; i++)
    {
        fscanf(radexout,"%s%s%s%s%s%s%s%s%s%le%le",s,s,s,s,s,s,s,s,s,\
               &(*(levels+2*i+0%2)),&(*(levels+2*i+1%2)));
        while(fgetc(radexout)!='\n');
    }

    fclose(radexout);
    return levels;
}

double* radexoutHp(int n)
{

    FILE* radexout=fopen("radex.out","r");

    char s[2000]= {0};
    int k;
    int i;
    for (k=0,i=0; k<12&&((s[i]=fgetc(radexout))!=EOF); i++)
        if(s[i]=='\n')
            k++;
    double *levels;
    levels=malloc(sizeof(double)*2*n);
    for(i=0; i<2*n; i++)
        *(levels+i)=0;

    for(i=0; i<n; i++)
    {
        fscanf(radexout,"%s%s%s%s%s%s%s%s%s%le%le",s,s,s,s,s,s,s,s,s,\
               &(*(levels+2*i+0%2)),&(*(levels+2*i+1%2)));
        while(fgetc(radexout)!='\n');
    }

    fclose(radexout);
    return levels;
}
