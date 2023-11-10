#include "exportFiles.h"
#include "findPT.h"
#include "fitModels.h"
#include "fitPolytropes.h"


/*
COMPILE with 
 gcc -o taspp main.c EOS_cold.c exportFiles.c findPT.c fitModels.c fitPolytropes.c -lm
*/

int main(){
    //Modify: path to table, name of table, path for files to export
    char pathTables[1000]="E:\\DOCTORADO\\REPO\\95_coldEoSTablesArizona\\fixed\\";
    char pathExport[1000]="C:\\Users\\marin\\Desktop\\TASPP\\";
    char eosName[100]="ps";
    
    //Read table
    char full_path[1000];
    full_path[0]='\0';
    struct Tab table;
    strcat(full_path,pathTables);
    strcat(full_path,eosName);
    strcat(full_path,".txt");
    printf("Working with table at %s\n",full_path);
    prepareTab(full_path, &table);

    //Identify PT
    double *xL, *xR; int n_PT;
    n_PT=obtainPTs(&xL, &xR, table); 
    printf("Number of PTs identified: %d\n",n_PT);
    for (int i=0; i<n_PT; i++) printf("     %e  %e\n", xL[i], xR[i]);

    //Identify the need of an extra polytrope
    double extra_cut=-1;
    if(n_PT==1) extra_cut=needExtraCut(table, xR[0]);
    printf("\nTASPP model does "); if(extra_cut<0) printf("NOT "); printf("need an extra polytrope\n");

    //Construct polytropes beyond crust outside of PTs
    printf("\nFitting polytropes\n");
    double *Kfitted, *Gfitted, rc;
    Kfitted=(double*)malloc((n_PT+1+SIGN(extra_cut))*sizeof(double));
    Gfitted=(double*)malloc((n_PT+1+SIGN(extra_cut))*sizeof(double));
    rc=fitPolytropes(table, n_PT, xL, xR, extra_cut, Kfitted, Gfitted);

    //Obtain polynomial models for PTs
    double *a, *b, *c, *f, *d, *xm;
    a=(double*)malloc(n_PT*sizeof(double));
    b=(double*)malloc(n_PT*sizeof(double));
    c=(double*)malloc(n_PT*sizeof(double));
    f=(double*)malloc(n_PT*sizeof(double));
    d=(double*)malloc(n_PT*sizeof(double));
    xm=(double*)malloc(n_PT*sizeof(double));
    fitPolynomials(table, n_PT, xL, xR, Kfitted, Gfitted, a, b, c, f, d, xm);

    //Export model in par file (valid for TOV solver and hydrodynamics)
    double *rho_PT; int count=0;
    rho_PT=(double*)malloc(n_PT*2*sizeof(double));
    for(int i=0; i<n_PT; i++) {rho_PT[count]=xL[i]; count++; rho_PT[count]=xR[i]; count++;}
    full_path[0]='\0';
    strcat(full_path, pathExport);
    strcat(full_path, eosName);
    strcat(full_path, "_TASPP.txt");
    exportEOS_TASPP(full_path, n_PT, d, xm, rho_PT, a, b, c, f, extra_cut, rc, Kfitted, Gfitted);

    struct TASPP eos;
    double **params;
    params=(double**)malloc(n_PT*sizeof(double*));
    for (int k=0; k<n_PT; k++){
        params[k]=(double*)malloc(4*sizeof(double));
        params[k][0]=a[k]; params[k][1]=b[k]; params[k][2]=c[k]; params[k][3]=f[k];
    }
    prepareTASPP(rho_PT, params, d, xm, n_PT, &eos, extra_cut, rc, Kfitted, Gfitted);
    full_path[0]='\0';
    strcat(full_path, pathExport);
    strcat(full_path, eosName);
    strcat(full_path, "_quant.txt");
    exportQuantities(full_path, rc, table.rhotab[table.lenEos-1], eos);


    free(a); free(b); free(c); free(d); free(f); free(xm);
    free(Kfitted); free(Gfitted); free(rho_PT);
    for (int k=0; k<n_PT; k++) free(params[k]);
    cleanTab(&table);





}