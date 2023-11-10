#ifndef EOS_COLD
#define EOS_COLD

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SIGN(a) ( ( (a)<(0) ) ? (0) :(1) )

#define SLY 1
#define FPS 2
#define BBB2 3
#define GS1 4
#define BGN1H1 5
#define H6 6
#define ALF1 7 
#define PAL6 8
#define BPAL12 9
#define MS2 10
#define GS2 11
#define GNH3 12
#define H4 13
#define H5 14
#define H7 15
#define ALF2 16
#define ALF3 17
#define ALF4 18
#define APR1 19
#define APR2 20
#define APR3 21
#define APR4 22
#define WFF1 23
#define WFF2 24
#define WFF3 25
#define ENG  26
#define MPA1 27
#define MS1 28
#define MS1B 29
#define PS 30
#define H1 31
#define H2 32
#define H3 33
#define PCL2 34

struct TASPP{
    double *a, *b, *c, *f, *d, *xm; //params of polynomial pieces (length n_PT)
    int n_PT; //number of PTs in the EoS
    double *G; //array of adiabatic index polytrope ordered by use 
    double *K; //array of constant K polytrope ordered by use 
    double *rmd; //array of dividing densities, increasing order
    int Nrmd; //length of densities vector
    int Nparam; //length of G and K vectors 
    int Npieces; //length of pieces vector
    int *pieces; //pieces[i]>0, contains index of G, K arrays to use for rho<rmd[i]. If negative, -(pieces+1) contains index of a,b,d,xm to use
    double *cteps; //constant for internal energy of polytrope
    double *Pcut; //Pressure at piece change
    int *i_ph; //contains index in rmd of xL of the PTs (length n_PT). -(pieces+1) also points here
};

struct PP{
    double *G; //array of adiabatic index polytrope ordered by use
    double *K; //array of constant K polytrope ordered by use
    double *rmd; //array of dividing densities, increasing order
    double *cteps; //constant for internal energy of polytrope
    double *Pcut; //Pressure at piece change
    int *pieces;// contains index of G, K arrays to use for rho<rmd[i].
    int Nrmd; //length of densities vector
    int Nparam; //length of G and K vectors
};

struct Tab{
    double *rhotab, *fundtab, *aditab, *Ptab, *epstab, *csctab, *csrtab, *grutab;
    //quantities to read from the table
    double *logrhotab, *logfundtab, *logaditab, *logPtab, *logepstab, *logcsctab, *logcsrtab, *loggrutab;
    //quantities calculated to optimize computations
    int lenEos; //number of density points in the table
    double minG; //contains the minimum value of the fundamental derivative along the table
};

int stringToEoS(char *name);
void preparePP(struct PP *eos,int which);
void prepareTASPP(double *rho_PT,double **params, double *denter,
                      double *xmenter,int n_PT, struct TASPP *eos, 
                      double extra_cut, double rc, double *Kfitted, double *Gfitted);
void prepareTab(char*namefile,struct Tab *eos);


void cleanTASPP(struct TASPP *eos);
void cleanPP(struct PP *eos);
void cleanTab(struct Tab *eos);

//for tabulated
void prepareTab(char *name, struct Tab *eos);
int locateQuantity(double *quant, double looked, int lenEos);
double extrapolateTab(int lenEos, double locate, double *quantGiven, double *quantGet);
double extrapolateGTab(int lenEos, double locate, double *quantGiven, double *quantGet, double minG);

//for polytropes
double epsilonpoly(double rho, double k, double G, double cte, double rp);
double Ppoly(double rho, double k, double G);
double cspoly(double rho, double k, double G, double cte, double rp);
double Gpoly(double G);

//for polynomial g234
int wherePT_g234(double x, struct TASPP eos);
double poli_g234(double x,struct TASPP eos, int pt);
double intpoli_g234(double x,struct TASPP eos, int pt);
double dpoli_g234(double x,struct TASPP eos, int pt);
double Pfromcs_g234(double x, struct TASPP eos, int pt);
double fundfromcs_g234(double x, struct TASPP eos, int pt);
double gammafromcs_g234(double x, struct TASPP eos, int pt);
double epsilonfromcs_g234(double x, struct TASPP eos, int pt);
double csrelfromcs_g234(double x, struct TASPP eos, int pt);

//TASPP
double pressureTASPP(double rho, struct TASPP eos);
double densityTASPP(double P, struct TASPP eos);
double epsilonTASPP(double rho, struct TASPP eos);
double hTASPP(double rho, struct TASPP eos);
double csTASPP(double rho, struct TASPP eos);
double cscTASPP(double rho, struct TASPP eos);
double GTASPP(double rho, struct TASPP eos);
double nonlfactorTASPP(double rho, struct TASPP eos);
double kfancyTASPP(double rho,  struct TASPP eos);
double adiabaticexponentTASPP(double rho, struct TASPP eos);

//PP
double pressurePP(double rho, struct PP eos);
double densityPP(double P, struct PP eos);
double epsilonPP(double rho, struct PP eos);
double hPP(double rho, struct PP eos);
double csPP (double rho, struct PP eos);
double cscPP (double rho, struct PP eos);
double GPP(double rho, struct PP eos);
double nonlfactorPP(double rho, struct PP eos);
double kfancyPP(double rho, double h,  struct PP eos);
double adiabaticExponentPP(double rho, struct PP eos);
double paramGPP(double rho, struct PP eos);
double paramKPP(double rho, struct PP eos);

//Tab
double pressureTab(double rho, struct Tab eos);
double densityTab(double P, struct Tab eos);
double epsilonTab(double rho, struct Tab eos);
double hTab(double rho, struct Tab eos);
double csTab(double rho, struct Tab eos);
double GTab(double rho, struct Tab eos);
double nonlfactorTab(double rho, struct Tab eos);
double kfancyTab(double rho, struct Tab eos);


double nonlfactorsignPP(double rho, struct PP eos);
double nonlfactorsignTASPP(double rho, struct TASPP eos);
double nonlfactorsignTab(double rho, struct Tab eos);

#endif