#include "exportFiles.h"

void exportEOS_TASPP(char *filename, int n_PT, double *d, double *xm, double *rho_PT, double *a, double *b, double *c, double *f, double extra_cut, double rc, double *Kfitted, double *Gfitted){
    FILE *fout; 
    fout=fopen(filename,"w");

    fprintf(fout, "#EoS type (PP, TASPP, Tabulated)\nTASPP\n\n");
    fprintf(fout,"#(PP) EoS name, (Tabulated) Path to table, (TASPP) number of PTs\n");
    fprintf(fout,"%d\n\n",n_PT);
    fprintf(fout, "#(Only TASPP from here on)\n#d\n");
    for(int i=0; i<n_PT-1; i++) {fprintf(fout,"%.15f\t", d[i]);} fprintf(fout,"%.15f\n\n", d[n_PT-1]);
    fprintf(fout,"#xm\n");
    for(int i=0; i<n_PT-1; i++) {fprintf(fout,"%.15e\t", xm[i]);} fprintf(fout,"%.15f\n\n", xm[n_PT-1]);
    fprintf(fout,"#rhoPT (left, right)\n");
    for(int i=0; i<2*n_PT-1; i++) {fprintf(fout,"%.15e\t", rho_PT[i]);} fprintf(fout,"%.15f\n\n", rho_PT[2*n_PT-1]);
    fprintf(fout,"#polynomial parameters\n");
    for(int i=0; i<n_PT; i++) {fprintf(fout,"%.15f\t%.15f\t%.15f\t%.15f\n", a[i], b[i], c[i], f[i]);} fprintf(fout, "\n\n");
    fprintf(fout, "#rho_cut and extra_cut (if there's none, specify -1)\n");
    fprintf(fout, "%.15e\t%.15e\n\n", rc, extra_cut);
    fprintf(fout, "#K polytropes\n");
    for(int i=0; i<2+SIGN(extra_cut)-1; i++) {fprintf(fout, "%.15e\t", Kfitted[i]);}fprintf(fout, "%.15e\n\n", Kfitted[2+SIGN(extra_cut)-1]);
    fprintf(fout, "#Gamma polytropes\n");
    for(int i=0; i<2+SIGN(extra_cut)-1; i++) {fprintf(fout, "%.15e\t", Gfitted[i]);}fprintf(fout, "%.15e\n\n", Gfitted[2+SIGN(extra_cut)-1]);

}

void exportQuantities(char *filename, double rc, double fin, struct TASPP eos){
    FILE *fout; 
    fout=fopen(filename,"w");

    double step=(log10(fin)-log10(rc))/5001;
    double rho;
    fprintf(fout,"#1-rho  2-G  3-h  4-P  5-csc\n");
    for (int j=0; j<5001; j++){
        rho=pow(10, log10(rc)+j*step);
        fprintf(fout, "%e %f %e %e %f\n", rho, GTASPP(rho, eos), hTASPP(rho,eos), pressureTASPP(rho, eos), cscTASPP(rho, eos));
    }
    fclose(fout);
}