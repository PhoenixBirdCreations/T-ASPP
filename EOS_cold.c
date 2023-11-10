#include "EOS_cold.h"


/*prepareTab
Description: Load tabulated EoS file and construct the EoS object
Input: path to the file, pointer to Tab object
Output: -
*/
void prepareTab(char *name, struct Tab *eos){ 
    FILE *fp = fopen(name, "r");
    
    const unsigned MAX_LENGTH = 256;
    char buffer[MAX_LENGTH], line2[MAX_LENGTH];
    int count=0;
    int element=0;
    char *token;

    int numberLines=-1;
    while (fgets(buffer, MAX_LENGTH, fp)) numberLines++;
    eos->lenEos=numberLines;
    printf("#File read: %s. #lines %d\n",name,numberLines);
    rewind(fp);
    eos->rhotab=(double *)malloc(numberLines*sizeof(double));
    eos->fundtab=(double *)malloc(numberLines*sizeof(double));
    eos->aditab=(double *)malloc(numberLines*sizeof(double));
    eos->Ptab=(double *)malloc(numberLines*sizeof(double));
    eos->epstab=(double *)malloc(numberLines*sizeof(double));
    eos->csctab=(double *)malloc(numberLines*sizeof(double));
    eos->csrtab=(double *)malloc(numberLines*sizeof(double));
    eos->grutab=(double *)malloc(numberLines*sizeof(double));
    eos->logrhotab=(double *)malloc(numberLines*sizeof(double));
    eos->logfundtab=(double *)malloc(numberLines*sizeof(double));
    eos->logaditab=(double *)malloc(numberLines*sizeof(double));
    eos->logPtab=(double *)malloc(numberLines*sizeof(double));
    eos->logepstab=(double *)malloc(numberLines*sizeof(double));
    eos->logcsctab=(double *)malloc(numberLines*sizeof(double));
    eos->logcsrtab=(double *)malloc(numberLines*sizeof(double));
    eos->loggrutab=(double *)malloc(numberLines*sizeof(double));

    fgets(buffer, MAX_LENGTH, fp); //skip header
    while (fgets(buffer, MAX_LENGTH, fp)){
        memcpy(line2, buffer, MAX_LENGTH);
        token = strtok(line2, "\t");
        element=0;
        while( token != NULL ) {
            if(element==0) eos->rhotab[count]=strtod(token,NULL);
            if(element==1) eos->fundtab[count]=strtod(token,NULL);
            if(element==2) eos->aditab[count]=strtod(token,NULL);
            if(element==3) eos->Ptab[count]=strtod(token,NULL);
            if(element==4) eos->epstab[count]=strtod(token,NULL);
            if(element==5) eos->csrtab[count]=strtod(token,NULL);
            if(element==6) eos->csctab[count]=strtod(token,NULL);
            if(element==7) eos->grutab[count]=strtod(token,NULL);
            token = strtok(NULL, "\t");
            element++;
        }
        count++;
    }
    eos->minG=50;
    for (register int i=0; i<numberLines; i++){
      if(eos->fundtab[i]<eos->minG){
        eos->minG=eos->fundtab[i];
      }
    } //if EoS is nonconvex, G needs a shift before taking logarithms
    if(eos->minG<0) eos->minG=fabs(eos->minG);
    else eos->minG=-1; //removes transformation

    for (register int i=0; i<numberLines; i++){
      eos->logrhotab[i]=log10(eos->rhotab[i]);
      eos->logfundtab[i]=log10(eos->fundtab[i]+eos->minG+1);
      eos->logaditab[i]=log10(eos->aditab[i]);
      eos->logPtab[i]=log10(eos->Ptab[i]);
      eos->logepstab[i]=log10(eos->epstab[i]);
      eos->logcsrtab[i]=log10(eos->csrtab[i]);
      eos->logcsctab[i]=log10(eos->csctab[i]);
      eos->loggrutab[i]=log10(eos->grutab[i]);
    }
fclose(fp);

}

/*locateQuantity
Description: locate a quantity in the tabulated data
Input: array of the tabulated quantity, number to look in it, length of the array
Output: index of the array after which the quantity is located
*/
int locateQuantity(double *quant, double looked, int lenEos){ 
    if(looked<quant[0]) return 0;
    register int i;
    for(i=0; i<lenEos-1; i++){
        if(quant[i]<=looked && looked<=quant[i+1])
            return i;
    }
    return i;
}

/*extrapolateTab
Description: logarithm interpolation of any tabulated quantity
Input: length of the table, point where to extrapolate, array of the quantity x for extrapolation, array of the quantity y for extrapolation
Output: extrapolated quantity
*/
double extrapolateTab(int lenEos, double locate, double *quantGiven, double *quantGet){
    int i=locateQuantity(quantGiven, log10(locate), lenEos);
    if(i==(lenEos-1)){   
        return pow(10,quantGet[i]+(quantGet[i]-quantGet[i-1])/(quantGiven[i]-quantGiven[i-1])*(log10(locate)-quantGiven[i]));
    }
    return pow(10,quantGet[i]+(quantGet[i+1]-quantGet[i])/(quantGiven[i+1]-quantGiven[i])*(log10(locate)-quantGiven[i]));
}

/*extrapolateGTab
Description: logarithm interpolation of fundamentald derivative (may need transformation)
Input: length of the table, point where to extrapolate, array of the quantity x for extrapolation, array of fundamental derivative, minimum value of G in the table
Output: extrapolated quantity
*/
double extrapolateGTab(int lenEos, double locate, double *quantGiven, double *quantGet, double minG){
    int i=locateQuantity(quantGiven, log10(locate), lenEos);
    if(i==(lenEos-1)){  
        return pow(10,quantGet[i]+(quantGet[i]-quantGet[i-1])/(quantGiven[i]-quantGiven[i-1])*(log10(locate)-quantGiven[i]))-1-minG;
    }
    return pow(10,quantGet[i]+(quantGet[i+1]-quantGet[i])/(quantGiven[i+1]-quantGiven[i])*(log10(locate)-quantGiven[i]))-1-minG;
}

/*parametersEoS
Description: implements the PP parameters from Read et al.
Input: name of EoS, pointer to the address of logP1, pointer to first adiabatic index, second, and third
Output:-
*/
void parametersEoS(int which, double *logp1, double *G1, double *G2, double *G3){
  if (which==SLY){ 
    *logp1=34.384;
    *G1=3.005;
    *G2=2.988;
    *G3=2.851;
    printf("PrepaprePP is using SLy\n");
  }
  else if(which==FPS){
    *logp1=34.283;
    *G1=2.985;
    *G2=2.863;
    *G3=2.600;
    printf("PrepaprePP is using FPS\n");
  }
  else if(which==BBB2){ 
    *logp1=34.331;
    *G1=3.418;
    *G2=2.835;
    *G3=2.832;
    printf("PrepaprePP is using BBB2\n");
    }
  else if(which==GS1){
    *logp1=34.504;
    *G1=2.350;
    *G2=1.267;
    *G3=2.421;
    printf("PrepaprePP is using GS1\n");
  }
  else if(which==BGN1H1){
    *logp1=34.623;
    *G1=3.258;
    *G2=1.472;
    *G3=2.464;
    printf("PrepaprePP is using BGN1H1\n");
  }
  else if (which==H6){
    *logp1=34.593;
    *G1=2.637;
    *G2=2.121;
    *G3=2.067;
    printf("PrepaprePP is using H6\n");
  }
  else if(which==ALF1){
    *logp1=34.055;
    *G1=2.013;
    *G2=3.389;
    *G3=2.033;
    printf("PrepaprePP is using ALF1\n");
  }
  else if(which==PAL6){
    *logp1=34.380;
    *G1=2.227;
    *G2=2.189;
    *G3=2.159;
    printf("PrepaprePP is using PAL6\n");
  }
  else if(which==BPAL12){
    *logp1=34.358;
    *G1=2.209;
    *G2=2.201;
    *G3=2.176;
    printf("PrepaprePP is using BPAL12\n");
  }
  else if(which==MS2){ 
    *logp1=34.605;
    *G1=2.447;
    *G2=2.184;
    *G3=1.855;
    printf("PrepaprePP is using MS2\n");
  }
  else if(which==GS2){
    *logp1=34.642;
    *G1=2.519;
    *G2=1.571;
    *G3=2.314;
    printf("PrepaprePP is using GS2\n");
  }
  else if(which==GNH3){ 
    *logp1=34.648;
    *G1=2.664;
    *G2=2.194;
    *G3=2.304;
    printf("PrepaprePP is using GNH3\n");
  }
  else if(which==H4){
    *logp1=34.669;
    *G1=2.909;
    *G2=2.246;
    *G3=2.144;
    printf("PrepaprePP is using H4\n");
  }
  else if(which==H5){ 
    *logp1=34.609;
    *G1=2.793;
    *G2=1.974;
    *G3=1.915;
    printf("PrepaprePP is using H5\n");
  }
  else if(which==H7){
    *logp1=34.559;
    *G1=2.621;
    *G2=2.048;
    *G3=2.006;
    printf("PrepaprePP is using H7\n");
  }
  else if(which==ALF2){
    *logp1=34.616;
    *G1=4.070;
    *G2=2.411;
    *G3=1.890;
    printf("PrepaprePP is using ALF2\n");
  }
  else if(which==ALF3){
    *logp1=34.283;
    *G1=2.883;
    *G2=2.653;
    *G3=1.952;
    printf("PrepaprePP is using ALF3\n");
  }
  else if(which==ALF4){
    *logp1=34.314;
    *G1=3.009;
    *G2=3.438;
    *G3=1.803;
    printf("PrepaprePP is using ALF4\n");
  }
  else if(which==APR1){
    *logp1=33.943;
    *G1=2.442;
    *G2=3.256;
    *G3=2.908;
    printf("PrepaprePP is using APR1\n");
  }
  else if(which==APR1){
    *logp1=33.943;
    *G1=2.442;
    *G2=3.256;
    *G3=2.908;
    printf("PrepaprePP is using APR1\n");
  }
  else if(which==APR2){
    *logp1=34.126;
    *G1=2.643;
    *G2=3.014;
    *G3=2.945;
    printf("PrepaprePP is using APR2\n");
  }
  else if(which==APR3){
    *logp1=34.392;
    *G1=3.166;
    *G2=3.573;
    *G3=3.281;
    printf("PrepaprePP is using APR3\n");
  }
  else if(which==APR4){
    *logp1=34.269;
    *G1=2.830;
    *G2=3.445;
    *G3=3.348;
    printf("PrepaprePP is using APR4\n");
  }
  else if(which==WFF1){
    *logp1=34.031;
    *G1=2.519;
    *G2=3.791;
    *G3=3.660;
    printf("PrepaprePP is using WFF1\n");
  }
  else if(which==WFF2){
    *logp1=34.233;
    *G1=2.888;
    *G2=3.475;
    *G3=3.517;
    printf("PrepaprePP is using WFF2\n");
  }
  else if(which==WFF3){
    *logp1=34.283;
    *G1=3.329;
    *G2=2.952;
    *G3=2.589;
    printf("PrepaprePP is using WFF3\n");
  }
  else if(which==ENG){
    *logp1=34.437;
    *G1=3.514;
    *G2=3.130;
    *G3=3.168;
    printf("PrepaprePP is using ENG\n");
  }
  else if(which==MPA1){
    *logp1=34.495;
    *G1=3.446;
    *G2=3.572;
    *G3=2.887;
    printf("PrepaprePP is using MPA1\n");
  }
  else if(which==MS1){
    *logp1=34.858;
    *G1=3.224;
    *G2=3.033;
    *G3=1.325;
    printf("PrepaprePP is using MS1\n");
  }
  else if(which==MS1B){
    *logp1=34.855;
    *G1=3.456;
    *G2=3.011;
    *G3=1.425;
    printf("PrepaprePP is using MS1B\n");
  }
  else if(which==PS){
    *logp1=34.671;
    *G1=2.216;
    *G2=1.640;
    *G3=2.365;
    printf("PrepaprePP is using PS\n");
  }
  else if(which==H1){
    *logp1=34.564;
    *G1=2.595;
    *G2=1.845;
    *G3=1.897;
    printf("PrepaprePP is using H1\n");
  }
  else if(which==H2){
    *logp1=34.617;
    *G1=2.775;
    *G2=1.855;
    *G3=1.858;
    printf("PrepaprePP is using H2\n");
  }
  else if(which==H3){
    *logp1=34.646;
    *G1=2.787;
    *G2=1.951;
    *G3=1.901;
    printf("PrepaprePP is using H3\n");
  }
  else if(which==PCL2){
    *logp1=34.507;
    *G1=2.554;
    *G2=1.880;
    *G3=1.977;
    printf("PrepaprePP is using PCL2\n");
  }
  else {
    printf("Please provide an EoS name that has been implemented\n");
    exit(2);
  }

}

/*stringToEoS
Description: for reading the type of EoS from the par file
Input: the string of the name of the EoS, all in caps
Output: the quantity defined in the program that refers to that EoS
*/
int stringToEoS(char *name){
  if(strstr(name, "SLY") != NULL) return SLY;
  if(strstr(name, "FPS") != NULL) return FPS;
  if(strstr(name, "BBB2") != NULL) return BBB2;
  if(strstr(name, "GS1") != NULL) return GS1;
  if(strstr(name, "BGN1H1") != NULL) return BGN1H1;
  if(strstr(name, "H6") != NULL) return H6;
  if(strstr(name, "ALF1") != NULL) return ALF1;
  if(strstr(name, "PAL6") != NULL) return PAL6;
  if(strstr(name, "BPAL12") != NULL) return BPAL12;
  if(strstr(name, "MS2") != NULL) return MS2;
  if(strstr(name, "GS2") != NULL) return GS2;
  if(strstr(name, "GNH3") != NULL) return GNH3;
  if(strstr(name, "H4") != NULL) return H4;
  if(strstr(name, "H5") != NULL) return H5;
  if(strstr(name, "H7") != NULL) return H7;
  if(strstr(name, "ALF2") != NULL) return ALF2;
  if(strstr(name, "ALF3") != NULL) return ALF3;
  if(strstr(name, "ALF4") != NULL) return ALF4;
  if(strstr(name, "APR1") != NULL) return APR1;
  if(strstr(name, "APR2") != NULL) return APR2;
  if(strstr(name, "APR3") != NULL) return APR3;
  if(strstr(name, "APR4") != NULL) return APR4;
  if(strstr(name, "WFF1") != NULL) return WFF1;
  if(strstr(name, "WFF2") != NULL) return WFF2;
  if(strstr(name, "WFF3") != NULL) return WFF3;
  if(strstr(name, "ENG") != NULL) return ENG;
  if(strstr(name, "MPA1") != NULL) return MPA1;
  if(strstr(name, "MS1") != NULL) return MS1;
  if(strstr(name, "MS1B") != NULL) return MS1B;
  if(strstr(name, "PS") != NULL) return PS;
  if(strstr(name, "H1") != NULL) return H1;
  if(strstr(name, "H2") != NULL) return H2;
  if(strstr(name, "H3") != NULL) return H3;
  if(strstr(name, "PCL2") != NULL) return PCL2;
  else return -1;
}
/*preparePP
Description: build the object of a PP EoS
Input: pointer to the PP object, name of EoS to use
Output: - 
*/
void preparePP(struct PP *eos, int which){
  double logp1, G1, G2, G3;
  parametersEoS(which, &logp1, &G1, &G2, &G3);
  double c2=898755178736817640000.0;
 
 //low density crust SLy4 
  double Kl3=6.8011e-9;
  double Gl3=1.58425;
  double rl3=2.44034e7;

  double Kl2=1.06186e-6;
  double Gl2=1.28733;
  double rl2=3.78358e11;

  double Kl1=53.2697;
  double Gl1=0.62223;
  double rl1=2.62780e12;

  double Gl=1.35692;
  double Kl=3.99874e-8;
 
 //High density selected
  double P=pow(10,logp1);
  double d1=pow(10,14.7);
  double d2=1.0e15;
  double K1=P/pow(d1,G1)/c2;
  double K2=P/pow(d1,G2)/c2;
  double K3=K2*pow(d2,G2)/pow(d2,G3); 
  double rl=pow(Kl/K1,1.0/(G1-Gl));

  eos->G=(double*)malloc(7*sizeof(double));
  eos->K=(double*)malloc(7*sizeof(double));
  eos->cteps=(double*)malloc(7*sizeof(double));
  eos->Pcut=(double*)malloc(6*sizeof(double));
  eos->rmd=(double*)malloc(6*sizeof(double));
  eos->pieces=(int*)malloc(7*sizeof(int));
  eos->Nrmd=6;
  eos->Nparam=7;

  eos->rmd[0]=rl3; eos->rmd[1]=rl2; eos->rmd[2]=rl1; eos->rmd[3]=rl; eos->rmd[4]=d1; eos->rmd[5]=d2;
  eos->G[0]=Gl3; eos->G[1]=Gl2; eos->G[2]=Gl1; eos->G[3]=Gl; eos->G[4]=G1; eos->G[5]=G2; eos->G[6]=G3;
  eos->K[0]=Kl3; eos->K[1]=Kl2; eos->K[2]=Kl1; eos->K[3]=Kl; eos->K[4]=K1; eos->K[5]=K2; eos->K[6]=K3;
  for(int i=0;i<7;i++) eos->pieces[i]=i;

  //Dividing pressures
  eos->Pcut[0]=Kl3*pow(rl3,Gl3);
  eos->Pcut[1]=Kl2*pow(rl2,Gl2);
  eos->Pcut[2]=Kl1*pow(rl1,Gl1);
  eos->Pcut[3]=K1*pow(rl,G1);
  eos->Pcut[4]=K1*pow(d1,G1);
  eos->Pcut[5]=K2*pow(d2,G2);

 //Construction of constant for internal energy
  eos->cteps[0]=0.0;
  eos->cteps[1]= eos->cteps[0]+Kl3/(Gl3-1)*pow(rl3,Gl3-1);
  eos->cteps[2]=eos->cteps[1]+Kl2/(Gl2-1)*(pow(rl2,Gl2-1)-pow(rl3,Gl2-1));
  eos->cteps[3]=eos->cteps[2]+Kl1/(Gl1-1)*(pow(rl1,Gl1-1)-pow(rl2,Gl1-1));
  eos->cteps[4]=eos->cteps[3]+Kl/(Gl-1)*(pow(rl,Gl-1)-pow(rl1,Gl-1));
  eos->cteps[5]=eos->cteps[4]+K1/(G1-1)*(pow(d1,G1-1)-pow(rl,G1-1));
  eos->cteps[6]=eos->cteps[5]+K2/(G2-1)*(pow(d2,G2-1)-pow(d1,G2-1));
}

/*prepareTASPP
Description: build the object of a TASPP EoS
Input: array of PT extremes, matrix of coefficient of polynomials, array d term of polynomials, array xm term of polynomials,
       number of PTs, pointer to the TASPP object, matching density for extra high density polytrope (-1 if none),
       matching density with crust model, array of K of polytropes fitted, array of Gamma of polytropes fitted
Output: - 
*/
void prepareTASPP(double *rho_PT,double **paramsenter, double *denter,
                      double *xmenter,int n_PT, struct TASPP *eos, 
                      double extra_cut, double rc, double *Kfitted, double *Gfitted){
  
  if(n_PT<1){
    printf("Please provide a valid n_PT. A TASPP model has at least 1 PT.\n");
    exit(3);
  }

 //low density crust SLy4 
  double Kl3=6.8011e-9;
  double Gl3=1.58425;
  double rl3=2.44034e7;

  double Kl2=1.06186e-6;
  double Gl2=1.28733;
  double rl2=3.78358e11;

  double Kl1=53.2697;
  double Gl1=0.62223;
  double rl1=2.62780e12;

  double Gl=1.35692;
  double Kl=3.99874e-8;

  //Entered params for the polynomials
  eos->a=(double *)malloc(n_PT*sizeof(double));
  eos->b=(double *)malloc(n_PT*sizeof(double));
  eos->c=(double *)malloc(n_PT*sizeof(double));
  eos->f=(double *)malloc(n_PT*sizeof(double));
  eos->d=(double *)malloc(n_PT*sizeof(double));
  eos->xm=(double *)malloc(n_PT*sizeof(double));
  for(int i=0; i<n_PT; i++){
    eos->a[i]=paramsenter[i][0];
    eos->b[i]=paramsenter[i][1];
    eos->c[i]=paramsenter[i][2];
    eos->f[i]=paramsenter[i][3];
  }
  memcpy(eos->d, denter, n_PT*sizeof(double));
  memcpy(eos->xm, xmenter, n_PT*sizeof(double));

  //Construction of densities in order
    //reserve memory of auxiliar vectors
  double *rho_aux, *rho_PP;
  rho_aux=(double*)malloc((4+SIGN(extra_cut)+2*n_PT)*sizeof(double));
  rho_PP=(double*)malloc((4+SIGN(extra_cut))*sizeof(double));
  rho_PP[0]=rl3; rho_PP[1]=rl2; rho_PP[2]=rl1; rho_PP[3]=rc; 
  if(extra_cut>0) rho_PP[4]=extra_cut;
  int *iL, *iR;
  iL=(int*)malloc(n_PT*sizeof(int));
  iR=(int*)malloc(n_PT*sizeof(int));

    //create array of all densities in order, save location of PTs
  int i_PT=0, i_PP=0, i_aux=0, i_loc=0;
  while(i_aux<(4+SIGN(extra_cut)+2*n_PT)){ 
      if(i_PT<n_PT*2){ //if there is a PT limit still to consider
      if(i_PP<4+SIGN(extra_cut)){ //if there is a polytrope limit still to consider. Aquí tiene que ser 4
        if(rho_PT[i_PT]<rho_PP[i_PP]){ //in order there is now a extreme of a PT
          rho_aux[i_aux]=rho_PT[i_PT];
          if((i_PT%2)==0){ //we save the index of the extreme as left...
            iL[i_loc]=i_aux;
          }
          else{ //... or right
            iR[i_loc]=i_aux;
            i_loc++; //only increases after exiting a PT
          }  
          i_aux++; //and increase counters
          i_PT++;
        }
        else{  //in order there is now a dividing density of the PP
          rho_aux[i_aux]=rho_PP[i_PP];
          i_aux++;
          i_PP++;
        }
      }
      else{ //there is a PT limit to consider but no a polytrope limit
        rho_aux[i_aux]=rho_PT[i_PT];
        if((i_PT%2)==0){ //we save the index of the extreme as left...
          iL[i_loc]=i_aux;
        }
        else{ //... or right
          iR[i_loc]=i_aux;
          i_loc++; //only increases after exiting a PT
        }  
        i_aux++; //and increase counters
        i_PT++;
      }
    }
    else{  //there is no more PT limit to consider, just polytropes
      rho_aux[i_aux]=rho_PP[i_PP];
      i_aux++;
      i_PP++;
    }
  }

//  /*DEBUG*/  printf("rho_aux: ");for(int k=0; k<4+2*n_PT; k++)printf("%e ",rho_aux[k]);printf("\n");
  
  eos->Nrmd=4+2*n_PT+SIGN(extra_cut);

    //save the final density vector. Fill location of PT in rmd
  eos->rmd=(double *)malloc(eos->Nrmd*sizeof(double));
  eos->i_ph=(int *)malloc((n_PT)*sizeof(int));
  eos->n_PT=n_PT;
  i_aux=0; i_loc=0; int i_final, i_rho=0; 
  for (i_final=0; i_final<eos->Nrmd; i_final++){
    if(i_loc<n_PT){
      if (i_aux>iL[i_loc]){
        eos->i_ph[i_rho]=i_final-1;
        i_rho++;
        i_aux=iR[i_loc];
        i_loc++;
      }
    }
    eos->rmd[i_final]=rho_aux[i_aux];
    i_aux++;
  }
  
//  /*DEBUG*/printf("rmd: ");for(int k=0; k<eos->Nrmd; k++)printf("%e ",eos->rmd[k]);printf("\n");
//  /*DEBUG*/printf("i_ph: ");for(int k=0; k<eos->n_PT; k++)printf("%d ",eos->i_ph[k]);printf("\n");
  
      //Fill the used params: G, K
  double *Gaux, *Kaux; int normal=4+1+n_PT+SIGN(extra_cut);
  Gaux=(double *)malloc(normal*sizeof(double)); 
  Kaux=(double *)malloc(normal*sizeof(double));
  Gaux[0]=Gl3; Gaux[1]=Gl2; Gaux[2]=Gl1; Gaux[3]=Gl; 
  Kaux[0]=Kl3; Kaux[1]=Kl2; Kaux[2]=Kl1; Kaux[3]=Kl; 
  for (int i=4; i<normal; i++){
    Gaux[i]=Gfitted[i-4];
    Kaux[i]=Kfitted[i-4];
  }

//  /*DEBUG*/ printf("Gaux: ");for (int i=0; i<normal; i++) printf("%f ",Gaux[i]); printf("\n");
  
  eos->Nparam=normal;
  eos->G=(double *)malloc(eos->Nparam*sizeof(double));
  eos->K=(double *)malloc(eos->Nparam*sizeof(double)); 
  memcpy(eos->G, Gaux, normal*sizeof(double));
  memcpy(eos->K, Kaux, normal*sizeof(double));
  
// /*DEBUG*/printf("params: ");for(int k=0; k<eos->Nparam; k++) printf("%f ",eos->G[k]);printf("\n");
    
    //Fill in the index of the pieces
  eos->Npieces=eos->Nparam+n_PT;
  eos->pieces=(int *)malloc((eos->Npieces)*sizeof(int));
  i_aux=0; i_loc=0; int p_PP=0, p_pol=-1;
  for (i_final=0; i_final<=eos->Nrmd; i_final++){
    if(i_loc<n_PT){
      if (i_aux>iL[i_loc]){
        i_aux=iR[i_loc];
        i_loc++;
        eos->pieces[i_final]=p_pol;
        p_pol--;
      }
      else{
        eos->pieces[i_final]=p_PP;
        p_PP++;
      }
    }
    else{
      eos->pieces[i_final]=p_PP;
      p_PP++;
    }
    i_aux++;
  }

//  /*DEBUG*/printf("pieces: ");for(int k=0; k<=eos->Nrmd; k++) printf("%d ",eos->pieces[k]);printf("\n");
 
    //Dividing pressures
  eos->Pcut=(double *)malloc(eos->Nrmd*sizeof(double));
  int p;
  i_aux=0;
  for(i_final=0; i_final<eos->Nrmd; i_final++){
    p=eos->pieces[i_aux];
    if(p<0) {
      i_aux++;
      p=eos->pieces[i_aux];
      i_aux--;
    }
    eos->Pcut[i_final]=eos->K[p]*pow(eos->rmd[i_final],eos->G[p]);
    i_aux++;
  }

// /*DEBUG*/printf("Pcut: ");for(int k=0; k<eos->Nrmd; k++)printf("%e ",eos->Pcut[k]);printf("\n");

 //Construction of constant for internal energy
  eos->cteps=(double *)malloc(eos->Nparam*sizeof(double));
  eos->cteps[0]=0.0;
  i_aux=0; i_loc=0; i_rho=0;
  double aux_rho;
  for (i_final=1; i_final<eos->Nparam; i_final++, i_aux++){
    p=eos->pieces[i_aux];
    if(eos->pieces[i_aux+1]>=0){
      if(i_rho>0) aux_rho=eos->rmd[i_rho-1]; else aux_rho=0.0;
      eos->cteps[i_final]=eos->cteps[i_final-1]+eos->K[p]/(eos->G[p]-1)*(pow(eos->rmd[i_rho],eos->G[p]-1)-pow(aux_rho,eos->G[p]-1));
      i_rho++;
    }
    else{
      i_rho++; i_aux++;
      eos->cteps[i_final]=epsilonfromcs_g234(eos->rmd[i_rho],*eos,-(eos->pieces[i_aux]+1));
      i_rho++;
    }
  }
 
// /*DEBUG*/printf("ctes: ");for(int k=0; k<eos->Nparam; k++) printf("%f ",eos->cteps[k]);printf("\n");

  free(rho_aux); free(Gaux); free(Kaux); free(rho_PP); free(iL); free(iR);
}

///Functions to free memory of the EoS object
void cleanTASPP(struct TASPP *eos){
  free(eos->a); free(eos->b); free(eos->c); free(eos->f);
  free(eos->d); free(eos->xm);
  free(eos->G); free(eos->K); free(eos->rmd);
  free(eos->pieces); free(eos->cteps); free(eos->Pcut); free(eos->i_ph);
}
void cleanPP(struct PP *eos){
  free(eos->G); free(eos->K); free(eos->rmd);
  free(eos->pieces); free(eos->cteps); free(eos->Pcut); 
}
void cleanTab(struct Tab *eos){
  free(eos->rhotab); free(eos->fundtab); free(eos->aditab);
  free(eos->Ptab); free(eos->epstab); free(eos->csctab); 
  free(eos->csrtab); free(eos->grutab); 
  free(eos->logrhotab); free(eos->logfundtab); free(eos->logaditab);
  free(eos->logPtab); free(eos->logepstab); free(eos->logcsctab); 
  free(eos->logcsrtab); free(eos->loggrutab); 
}

//functions for polytropes
double epsilonpoly(double rho, double k, double G, double cte, double rp){
    return cte+k/(G-1)*(pow(rho,G-1)-pow(rp,G-1));
}
double Ppoly(double rho, double k, double G){
  return k*pow(rho,G);
}
double cspoly(double rho, double k, double G, double cte, double rp){
  double P=Ppoly(rho,k,G);
  return G*P/(rho+P+rho*epsilonpoly(rho,k,G,cte,rp));
}
double Gpoly(double G){
  return (G+1)*0.5;
}



//functions for polynomial g234

/*wherePT_g234
Description: locate the index of the piece of the PT where a density value is
Input: density value, EoS object
Output: index of pieces array 
*/
int wherePT_g234(double x, struct TASPP eos){
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(x>=eos.rmd[i]){
          j=eos.pieces[i+1];
          break;
      }
  }
  return -(j+1);
}

/*poli_g234
Description: evaluate a polynomial piece at a density value
Input: density value, EoS object, index of piece
Output: value of the polynomial (cs^2) 
*/
double poli_g234(double x, struct TASPP eos, int pt){
  if(pt<0) pt=wherePT_g234(x,eos);
  double xL=eos.rmd[eos.i_ph[pt]], xR=eos.rmd[eos.i_ph[pt]+1], dif=xR-xL;
  double x2=x*x, xl2=xL*xL, xl4=xl2*xl2, xr2=xR*xR, xr4=xr2*xr2, xm2=eos.xm[pt]*eos.xm[pt];
  return eos.d[pt]+
          eos.f[pt]*(x-(xr2-xl2)/(dif*2))/eos.xm[pt]+
          eos.c[pt]*(x2-(xr2*xR-xl2*xL)/(dif*3))/xm2+
          eos.b[pt]*(x2*x-(xr4-xl4)/(dif*4))/(xm2*eos.xm[pt])+
          eos.a[pt]*(x2*x2-(xr4*xR-xl4*xL)/(dif*5))/(xm2*xm2);
}

/*dpoli_g234
Description: evaluate the derivative of a polynomial piece at a density value
Input: density value, EoS object, index of piece
Output: value of the derivative of the polynomial
*/
double dpoli_g234(double x, struct TASPP eos, int pt){
  if(pt<0) pt=wherePT_g234(x,eos);
  double xm=eos.xm[pt], xm2=xm*xm;
  return eos.a[pt]*4*(x*x*x)/(xm2*xm2)+
          eos.b[pt]*3*x*x/(xm2*xm)+
          eos.c[pt]*2*x/xm2+
          eos.f[pt]/xm;
}

/*Pfromcs_g234
Description: evaluate the pressure from a polynomial piece at a density value
Input: density value, EoS object, index of piece
Output: value of the pressure 
*/
double Pfromcs_g234(double x, struct TASPP eos, int pt){
  if(pt<0) pt=wherePT_g234(x,eos);
  double t1, t2, t3, t4, t5; int indexPT=eos.i_ph[pt];
  double xL=eos.rmd[indexPT], xR=eos.rmd[indexPT+1];
  double xm2=eos.xm[pt]*eos.xm[pt], xl2=xL*xL, xl4=xl2*xl2, xr2=xR*xR, xr4=xr2*xr2, x2=x*x, x4=x2*x2, xmxl=x-xL, factor=xmxl/(xR-xL);
  t1=eos.d[pt]*xmxl;
  t2=eos.f[pt]/eos.xm[pt]*((x2-xl2)*0.5-(xr2-xl2)*factor/(2));
  t3=eos.c[pt]/xm2*((x2*x-xl2*xL)/3-(xr2*xR-xl2*xL)*factor/(3));
  t4=eos.b[pt]/(xm2*eos.xm[pt])*((x4-xl4)*0.25-(xr4-xl4)*factor/(4));
  t5=eos.a[pt]/(xm2*xm2)*((x4*x-xl4*xL)*0.2-(xr4*xR-xl4*xL)*factor/(5));
  
  return eos.Pcut[indexPT]+t1+t2+t3+t4+t5;
  
}

/*fundfromcs_g234
Description: evaluate the fundamental derivative from a polynomial piece at a density value
Input: density value, EoS object, index of piece
Output: value of the fundamental derivative 
*/
double fundfromcs_g234(double x, struct TASPP eos, int pt){
    if(pt<0) pt=wherePT_g234(x,eos);
    return 1.0+0.5*x*dpoli_g234(x,eos,pt)/poli_g234(x,eos,pt);
}

/*gammafromcs_g234
Description: evaluate the adiabatic index from a polynomial piece at a density value
Input: density value, EoS object, index of piece
Output: value of the adiabatic index
*/
double gammafromcs_g234(double x, struct TASPP eos, int pt){
  if(pt<0) pt=wherePT_g234(x,eos);
  return x*poli_g234(x,eos,pt)/Pfromcs_g234(x,eos,pt);
}

/*epsilonfromcs_g234
Description: evaluate the internal energy from a polynomial piece at a density value
Input: density value, EoS object, index of piece
Output: value of the internal energy
*/
double epsilonfromcs_g234(double x, struct TASPP eos, int pt){
  if(pt<0) pt=wherePT_g234(x,eos);
  double xL=eos.rmd[eos.i_ph[pt]], xR=eos.rmd[eos.i_ph[pt]+1];
  int p; p=eos.pieces[eos.i_ph[pt]];
  double alfa, beta, xi, delta;

  double xr2=xR*xR, xr4=xr2*xr2, xl2=xL*xL, xl3=xL*xl2, xl4=xl2*xl2, xl5=xl4*xL;
  double dif=xR-xL;
  double xm=eos.xm[pt], xm2=xm*xm;
  double x2=x*x, x4=x2*x2;

  alfa=(xr4*xR-xl5)/(5*dif);
  beta=(xr4-xl4)/(4*dif);
  xi=(xr2*xR-xl3)/(3*dif);
  delta=(xr2-xl2)/(2*dif);

  double factor=x*log(xL/x)+(x-xL);

  return epsilonpoly(xL,eos.K[p],eos.G[p],eos.cteps[p],eos.rmd[eos.i_ph[pt]-1])
          -eos.d[pt]*factor/x
          +eos.f[pt]/(2*xm*x)*(xl2-2*xL*x+x2+2*delta*factor)
          +eos.c[pt]/(6*x*xm2)*(2*xl3-3*xl2*x+x2*x+6*xi*factor)
          +eos.b[pt]/(12*x*xm2*xm)*(3*xl4-4*xl3*x+x4+12*beta*factor)
          +eos.a[pt]/(20*x*xm2*xm2)*(4*xl5-5*xl4*x+x4*x+20*alfa*factor)
          -Ppoly(xL,eos.K[p],eos.G[p])*(xL-x)/(xL*x);
}

/*epsilonfromcs_g234
Description: evaluate the relativistic sound speed from a polynomial piece at a density value
Input: density value, EoS object, index of piece
Output: value of the relativistic sound speed
*/
double csrelfromcs_g234(double x, struct TASPP eos, int pt){
  if(pt<0) pt=wherePT_g234(x,eos);
  double eps=epsilonfromcs_g234(x,eos,pt);
  double P=Pfromcs_g234(x,eos,pt);
  double h=1+eps+P/x;
  return poli_g234(x,eos,pt)/h;
}


//functions for TASPP. The functions calculate the quantity of their names, taking as input density and the EoS object
double pressureTASPP(double rho, struct TASPP eos){
  int i,j;        
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          if(j<0){
            return Pfromcs_g234(rho,eos,-(j+1));
          }
          else{
            return Ppoly(rho,eos.K[j],eos.G[j]);
          } 
      }
  }
  j=eos.pieces[0];
  return Ppoly(rho,eos.K[j],eos.G[j]);
}
double epsilonTASPP(double rho, struct TASPP eos){
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          if(j<0)
            return epsilonfromcs_g234(rho,eos,-(j+1));
          else
            return epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],eos.rmd[i]);
      }
  }
  j=eos.pieces[0];
  return epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],0.0);
}
double hTASPP(double rho, struct TASPP eos){
  int i,j; double eps,P=-3.0;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          if(j<0){
            eps=epsilonfromcs_g234(rho,eos,-(j+1));
            P=Pfromcs_g234(rho,eos,-(j+1));    
            break;
          }
          else{
              eps=epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],eos.rmd[i]);
              P=Ppoly(rho,eos.K[j],eos.G[j]);
              break;
          }
      }
  }
  j=eos.pieces[0];
  if(P<0){
      eps=epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],0.0);
      P=Ppoly(rho,eos.K[j],eos.G[j]);
  }
  return 1.0+eps+P/rho;

}
double csTASPP(double rho, struct TASPP eos){ 
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          if(j<0)
            return csrelfromcs_g234(rho,eos,-(j+1));
          else
            return cspoly(rho,eos.K[j],eos.G[j],eos.cteps[j],eos.rmd[i]);
      }
  }
  j=eos.pieces[0];
  return cspoly(rho,eos.K[j],eos.G[j],eos.cteps[j],0.0);
}
double cscTASPP(double rho, struct TASPP eos){
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          if(j<0)
            return poli_g234(rho,eos,-(j+1));
          else
            return Ppoly(rho,eos.K[j],eos.G[j])/rho*eos.G[j];  
      }
  }
  j=eos.pieces[0];
  return Ppoly(rho,eos.K[j],eos.G[j])/rho*eos.G[j];
}
double GTASPP(double rho, struct TASPP eos){
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          if(j<0)
            return fundfromcs_g234(rho,eos,-(j+1));
          else
            return Gpoly(eos.G[j]);
      }
  }
  j=eos.pieces[0];
  return Gpoly(eos.G[j]);
}
double nonlfactorTASPP(double rho, struct TASPP eos){
  double css=csTASPP(rho,eos);
  return  sqrt(css)/rho*(GTASPP(rho,eos)-1.5*css);
}
double kfancyTASPP(double rho,struct TASPP eos){ 
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          if(j<0)
            return 1+Pfromcs_g234(rho,eos,-(j+1))/rho/(1+epsilonfromcs_g234(rho,eos,-(j+1)));
          else
            return 1+(Ppoly(rho,eos.K[j],eos.G[j])/rho)/(1.0+epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],eos.rmd[i]));
      }
  }
  j=eos.pieces[0];
  return 1+(Ppoly(rho,eos.K[j],eos.G[j])/rho)/(1.0+epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],0.0));
}
double adiabaticexponentTASPP(double rho, struct TASPP eos){
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          if(j<0)
            return gammafromcs_g234(rho,eos,-(j+1));
          else
            return eos.G[j];
      }
  }
  j=eos.pieces[0];
  return eos.G[j];
}
double nonlfactorsignTASPP(double rho, struct TASPP eos){
  double css=csTASPP(rho,eos);
  return  (GTASPP(rho,eos)-1.5*css);
}

//functions for PP. The functions calculate the quantity of their names, taking as input density and the EoS object
double pressurePP(double rho, struct PP eos){
  int i,j;
  for(i=0; i<eos.Nrmd;i++){
      if(rho<eos.rmd[i]){
          j=eos.pieces[i];
          return Ppoly(rho,eos.K[j],eos.G[j]);
      }
  }
  j=eos.pieces[i];
  return Ppoly(rho,eos.K[j],eos.G[j]);
 
}
double densityPP(double P, struct PP eos){
    int i,j;
  for(i=0; i<eos.Nrmd;i++){
      if(P<eos.Pcut[i]){
          j=eos.pieces[i];
          return pow(P/eos.K[j],1.0/eos.G[j]);
      }
  }
  j=eos.pieces[i];
  return pow(P/eos.K[j],1.0/eos.G[j]);
}
double epsilonPP(double rho, struct PP eos){
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          return epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],eos.rmd[i]);
      }
  }
  j=eos.pieces[0];
  return epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],0.0);
}
double hPP(double rho, struct PP eos){
  double eps,P=-3.0;
    int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          eps=epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],eos.rmd[i]);
          P=Ppoly(rho,eos.K[j],eos.G[j]);
          break;
      }
  }
  if(P<0){
      j=eos.pieces[0];
      eps=epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],0.0); 
     P=Ppoly(rho,eos.K[j],eos.G[j]);
  }
  
  return 1.0+eps+P/rho;
}
double csPP (double rho, struct PP eos){ 
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          return cspoly(rho,eos.K[j],eos.G[j],eos.cteps[j],eos.rmd[i]);
      }
  }
  j=eos.pieces[0];
  return cspoly(rho,eos.K[j],eos.G[j],eos.cteps[j],0.0);
}
double cscPP (double rho, struct PP eos){
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          return Ppoly(rho,eos.K[j],eos.G[j])/rho*eos.G[j];
      }
  }
  j=eos.pieces[0];
  return Ppoly(rho,eos.K[j],eos.G[j])/rho*eos.G[j];
}
double GPP(double rho, struct PP eos){
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          return Gpoly(eos.G[j]);
      }
  }
  j=eos.pieces[0];
  return Gpoly(eos.G[j]);
}
double nonlfactorPP(double rho, struct PP eos){
  double css=csPP(rho,eos);
  return  sqrt(css)/rho*(GPP(rho,eos)-1.5*css);
}
double kfancyPP(double rho, double h, struct PP eos){
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          return h/(1.0+epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],eos.rmd[i]));
      }
  }
  j=eos.pieces[0];
  return h/(1.0+epsilonpoly(rho,eos.K[j],eos.G[j],eos.cteps[j],0.0));
}
double paramGPP(double rho, struct PP eos){
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          return eos.G[j];
      }
  }
  j=eos.pieces[0];
  return eos.G[j];
}
double paramKPP(double rho, struct PP eos){
  int i,j;
  for(i=eos.Nrmd-1; i>=0;i--){
      if(rho>=eos.rmd[i]){
          j=eos.pieces[i+1];
          return eos.K[j];
      }
  }
  j=eos.pieces[0];
  return eos.K[j];
}
double nonlfactorsignPP(double rho, struct PP eos){
  double css=csPP(rho,eos);
  return  (GPP(rho,eos)-1.5*css);
}


//functions for tabulated. The functions calculate the quantity of their names, taking as input density and the EoS object
double pressureTab(double rho, struct Tab eos){
  return extrapolateTab(eos.lenEos, rho, eos.logrhotab, eos.logPtab);
}
double densityTab(double P, struct Tab eos){
  return extrapolateTab(eos.lenEos, P, eos.logPtab, eos.logrhotab);
}
double epsilonTab(double rho, struct Tab eos){
  double eps=extrapolateTab(eos.lenEos, rho, eos.logrhotab, eos.logepstab);
  if (eps<0) return 0; else return eps;
}
double hTab(double rho, struct Tab eos){
  double eps=epsilonTab(rho,eos);
  double P=extrapolateTab(eos.lenEos, rho, eos.logrhotab, eos.logPtab);
  return 1+eps+P/rho;
}
double csTab(double rho, struct Tab eos){
  return extrapolateTab(eos.lenEos, rho, eos.logrhotab, eos.logcsrtab);
}
double GTab(double rho, struct Tab eos){
  return extrapolateGTab(eos.lenEos, rho, eos.logrhotab, eos.logfundtab, eos.minG);
}
double nonlfactorTab(double rho, struct Tab eos){
  double cs=csTab(rho,eos);
  if(cs<0) {printf("cs<0 in tabular data!!\n"); exit(1);}
  return sqrt(cs)/rho*(GTab(rho,eos)-1.5*cs);
}
double kfancyTab(double rho, struct Tab eos){
  double k=extrapolateTab(eos.lenEos, rho, eos.logrhotab, eos.loggrutab);
  return k/(k-csTab(rho,eos));
}
double nonlfactorsignTab(double rho, struct Tab eos){
  double cs=csTab(rho,eos);
  if(cs<0) {printf("cs<0 in tabular data!!\n"); exit(1);}
  return (GTab(rho,eos)-1.5*cs);
}



