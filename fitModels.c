#include "fitModels.h"


int minvector(int N, double *v){
    int k; double min=1e90;
    for(int i=0; i<N; i++){
        if (v[i]<min){
            min=v[i];
            k=i;
        }
    }
    return k;
}
int restructure(double left, double right, struct Tab table, double **rhoaux, double **cscaux, double **Paux){
    int imin,imax;

    for(int i=0; i<table.lenEos; i++) if (table.rhotab[i]>=left) {imin=i; break;}
    for(int i=imin; i<table.lenEos; i++) if (table.rhotab[i]>=right) {imax=i; break;}
    
    int N=imax-imin+1;
    
    if(((*rhoaux)=(double *)malloc(N*sizeof(double)))==NULL) {printf("error rho again\n"); exit(1);}
    if(((*cscaux)=(double *)malloc(N*sizeof(double)))==NULL) {printf("error csc again\n"); exit(1);}
    if(((*Paux)=(double *)malloc(N*sizeof(double)))==NULL) {printf("error P again\n"); exit(1);}
    for(int i=0; i<N; i++){
        (*rhoaux)[i]=table.rhotab[imin+i];
        (*cscaux)[i]=table.csctab[imin+i];
        (*Paux)[i]=table.Ptab[imin+i];
    }

    return N;

}
int find_min(double *csc, int N){
    double min=3, index=-1;
    for(int i=0; i<N;i++){
        if(csc[i]<min) {min=csc[i]; index=i;}
    }
    return index;
}
double Peasy(double x, double *params, double d, double xm, double xL, double xR, double PL){
    double t1=params[0]/pow(xm,4)*((pow(x,5)-pow(xL,5))/5-(pow(xR,5)-pow(xL,5))*(x-xL)/(5*(xR-xL)));
    double t2=params[1]/pow(xm,3)*((pow(x,4)-pow(xL,4))/4-(pow(xR,4)-pow(xL,4))*(x-xL)/(4*(xR-xL)));
    double t3=params[2]/pow(xm,2)*((pow(x,3)-pow(xL,3))/3-(pow(xR,3)-pow(xL,3))*(x-xL)/(3*(xR-xL)));
    double t4=params[3]/xm*((pow(x,2)-pow(xL,2))/2-(pow(xR,2)-pow(xL,2))*(x-xL)/(2*(xR-xL)));
    double t5=d*(x-xL);
    return PL+t1+t2+t3+t4+t5;
}

double c_fit(double x, double xm, double PR, double PL, double xL, double xR, double *params){
    double dif=xR-xL;
    return (PR-PL)/dif+
            params[0]/pow(xm,4)*(pow(x,4)-(pow(xR,5)-pow(xL,5))/(5*dif))+
            params[1]/pow(xm,3)*(pow(x,3)-(pow(xR,4)-pow(xL,4))/(4*dif))+
            params[2]/pow(xm,2)*(pow(x,2)-(pow(xR,3)-pow(xL,3))/(3*dif))+
            params[3]/xm*(x-(pow(xR,2)-pow(xL,2))/(2*dif));
}

double cost(double *rho, double *csc, int N, double xm, double PR, double PL, double xL, double xR, double *params){
    double sum=0.0, cvalue;
    int i;
    for(i=0; i<N;i++) {
        cvalue=c_fit(rho[i],xm,PR,PL,xL,xR,params);
        sum+=(csc[i]-cvalue)*(csc[i]-cvalue)/(csc[i]*csc[i]);
    }
    sum*=0.5;
    return sum;
}
double cost_P(double *rho, double *P, int N, double *params, double d, double xm, double xL, double xR, double PL){
    double sum=0.0, cvalue;
    int i;
    for(i=0; i<N;i++) {
        cvalue=Peasy(rho[i],params,d,xm,xL,xR,PL);
        sum+=(P[i]-cvalue)*(P[i]-cvalue)/(P[i]*P[i]);
    }
    sum*=0.5;
    return sum;
}
double cost_notrelative(double *rho, double *csc, int N, double xm, double PR, double PL, double xL, double xR, double *params){
    double sum=0.0, cvalue;
    int i;
    for(i=0; i<N;i++) {
        cvalue=c_fit(rho[i],xm,PR,PL,xL,xR,params);
        sum+=(csc[i]-cvalue)*(csc[i]-cvalue);
    }
    sum*=0.5;
    return sum;
}

double getabc_parabola(double *params, double xL,double xR,double PL,double PR,double xm,double CL,double CR,double CM,
             int N, double *rho, double *csc, double *P){
    
    double D, xa, xb, xc; double coste[6],costeP[6];
    D=(PR-PL)/(xR-xL); xa=(pow(xR,3)-pow(xL,3))/(3*(xR-xL));
    xb=(pow(xR,2)-pow(xL,2))/(2*(xR-xL));

    double **localparams;
    localparams=(double **)malloc(6*sizeof(double*));
    for(int j=0; j<6;j++){
        localparams[j]=(double *)malloc(4*sizeof(double));
        localparams[j][0]=0.0;
        localparams[j][1]=0.0;
    }

    double Cmu,xmu, Cnu,xnu;double ybr,yba, dem;
    //OP1: igualamos en CL, CM
    Cmu=CL; xmu=xL; Cnu=CM; xnu=xm;
    dem=(xmu-xb)/xm;
    ybr=(Cmu-D)/dem;
    yba=-(xmu*xmu-xa)/(xm*xm)/dem;
    dem=(xnu*xnu-xa)/(xm*xm)+yba*(xnu-xb)/xm;
    localparams[0][2]=(Cnu-D-ybr*(xnu-xb)/xm)/dem;
    localparams[0][3]=ybr+localparams[0][2]*yba;
    coste[0]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[0]);
    costeP[0]=cost_P(rho,P,N,localparams[0],D,xm,xL,xR,PL);

    //OP2: igualamos en CL, CR
    Cmu=CL; xmu=xL; Cnu=CR; xnu=xR;
    dem=(xmu-xb)/xm;
    ybr=(Cmu-D)/dem;
    yba=-(xmu*xmu-xa)/(xm*xm)/dem;
    dem=(xnu*xnu-xa)/(xm*xm)+yba*(xnu-xb)/xm;
    localparams[1][2]=(Cnu-D-ybr*(xnu-xb)/xm)/dem;
    localparams[1][3]=ybr+localparams[1][2]*yba;
    coste[1]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[1]);
    costeP[1]=cost_P(rho,P,N,localparams[1],D,xm,xL,xR,PL);

    //OP3: igualamos en CM, CR
    Cmu=CR; xmu=xR; Cnu=CM; xnu=xm;
    dem=(xmu-xb)/xm;
    ybr=(Cmu-D)/dem;
    yba=-(xmu*xmu-xa)/(xm*xm)/dem;
    dem=(xnu*xnu-xa)/(xm*xm)+yba*(xnu-xb)/xm;
    localparams[2][2]=(Cnu-D-ybr*(xnu-xb)/xm)/dem;
    localparams[2][3]=ybr+localparams[2][2]*yba;
    coste[2]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[2]);
    costeP[2]=cost_P(rho,P,N,localparams[2],D,xm,xL,xR,PL);

    //OP4: igualamos en CL, vértice en xm
    Cmu=CL; xmu=xL; 
    dem=(xmu-xb)/xm;
    ybr=(Cmu-D)/dem;
    yba=-(xmu*xmu-xa)/(xm*xm)/dem;
    localparams[3][2]=(-ybr)/(2+yba);
    localparams[3][3]=ybr+localparams[3][2]*yba;
    coste[3]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[3]);
    costeP[3]=cost_P(rho,P,N,localparams[3],D,xm,xL,xR,PL);

    //OP5: igualamos en CM, vértice en xm
    Cmu=CM; xmu=xm; 
    dem=(xmu-xb)/xm;
    ybr=(Cmu-D)/dem;
    yba=-(xmu*xmu-xa)/(xm*xm)/dem;
    localparams[4][2]=(-ybr)/(2+yba);
    localparams[4][3]=ybr+localparams[4][2]*yba;
    coste[4]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[4]);
    costeP[4]=cost_P(rho,P,N,localparams[4],D,xm,xL,xR,PL);

    //OP6: igualamos en CR, vértice en xm
    Cmu=CR; xmu=xR; 
    dem=(xmu-xb)/xm;
    ybr=(Cmu-D)/dem;
    yba=-(xmu*xmu-xa)/(xm*xm)/dem;
    localparams[5][2]=(-ybr)/(2+yba);
    localparams[5][3]=ybr+localparams[5][2]*yba;
    coste[5]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[5]);
    costeP[5]=cost_P(rho,P,N,localparams[5],D,xm,xL,xR,PL);

    //Checking thermodynamic consistency
    double xv,yv,testhere;
    for (int model=0;model<6;model++){
        xv=-(localparams[model][3]*xm)/localparams[model][2]*0.5;
        if (xv<xL) testhere=rho[0];
        else if (xv>xR) testhere=rho[N-1];
        else testhere=xv;
        yv=c_fit(testhere,xm,PR,PL,xL,xR,localparams[model]);
        if(yv<0){
            printf("Model %d thermo inconsistent: at %e, value %e\n",model+1,xv,yv);
            coste[model]=1e30;
            costeP[model]=1e50;
        }    
    }

    //causality: acausal waves will be killed by viscosity of the Riemann solver

    //leading term in parabola positive to ensure shape
    for (int model=0;model<6;model++){
        if(localparams[model][2]<0){
            coste[model]=1e29;
            costeP[model]=1e49;
        }
    }       

    double costeMedio[6]; for(int k=0;k<6;k++) costeMedio[k]=(coste[k]+costeP[k])*0.5;
    int this=minvector(6,costeMedio);
    
    params[0]=0.0;
    params[1]=0.0;
    params[2]=localparams[this][2];
    params[3]=localparams[this][3];

    return costeMedio[this];

}   

double getabc_cubic(double *params, double xL,double xR,double PL,double PR,double xm,double CL,double CR,double CM,
             int N, double *rho, double *csc, double *P){
    double D, xa, xb, xc; double coste[10], costeP[10];
    
    D=(PR-PL)/(xR-xL); xa=(pow(xR,4)-pow(xL,4))/(4*(xR-xL));
    xb=(pow(xR,3)-pow(xL,3))/(3*(xR-xL)); xc=(pow(xR,2)-pow(xL,2))/(2*(xR-xL));
    double yar,yab,yac,ybr, yba, ybc, dem;
    double Cmu,xmu, Cnu,xnu, xalfa;

    double **localparams;
    localparams=(double **)malloc(10*sizeof(double*));
    for(int j=0; j<10;j++){
        localparams[j]=(double *)malloc(4*sizeof(double));
        localparams[j][0]=0.0;
    }

    //OP1: CL, CM, CR
    yar=(CL-D)*pow(xm,3)/(pow(xL,3)-xa);
    yab=(xL*xL-xb)*xm/(pow(xL,3)-xa);
    yac=(xL-xc)*xm*xm/(pow(xL,3)-xa);
    ybr=((CR-D)*xm*xm-yar*(xR*xR*xR-xa)/xm)/(xR*xR-xb-yab*(xR*xR*xR-xa)/xm);
    ybc=(yac*(xR*xR*xR-xa)/xm-xm*(xR-xc))/(xR*xR-xb-yab*(xR*xR*xR-xa)/xm);
    localparams[0][3]=(CM-D-(1-xa/pow(xm,3))*(yar-yab*ybr)-ybr*(1-xb/xm/xm))/((1-xa/pow(xm,3))*(-yab*ybc-yac)+(1-xb/xm/xm)*ybc+1-xc/xm);
    localparams[0][2]=ybr+ybc*localparams[0][3];
    localparams[0][1]=yar-localparams[0][2]*yab-localparams[0][3]*yac;
    coste[0]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[0]);
    costeP[0]=cost_P(rho,P,N,localparams[0],D,xm,xL,xR,PL);

    //OP2: extremos, CL
    Cmu=CL; xmu=xL;
    dem=(pow(xmu,3)-xa)/pow(xm,3)+1.5*(xL*xL-xm*xm)*(xmu*xmu-xb)/(pow(xm,3)*(xm-xL))-3*(xmu-xc)/xm*(1+(xL*xL-xm*xm)/(xm*(xm-xL)));
    localparams[1][1]=(Cmu-D)/dem;
    localparams[1][2]=1.5*localparams[1][1]*(xL*xL-xm*xm)/(xm*(xm-xL));
    localparams[1][3]=-3*localparams[1][1]-2*localparams[1][2];
    coste[1]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[1]);
    costeP[1]=cost_P(rho,P,N,localparams[1],D,xm,xL,xR,PL);

    //OP3: extremos, CM
    Cmu=CM; xmu=xm;
    dem=(pow(xmu,3)-xa)/pow(xm,3)+1.5*(xL*xL-xm*xm)*(xmu*xmu-xb)/(pow(xm,3)*(xm-xL))-3*(xmu-xc)/xm*(1+(xL*xL-xm*xm)/(xm*(xm-xL)));
    localparams[2][1]=(Cmu-D)/dem;
    localparams[2][2]=1.5*localparams[2][1]*(xL*xL-xm*xm)/(xm*(xm-xL));
    localparams[2][3]=-3*localparams[2][1]-2*localparams[2][2];
    coste[2]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[2]);
    costeP[2]=cost_P(rho,P,N,localparams[2],D,xm,xL,xR,PL);

    //OP4: extremos, CR
    Cmu=CR; xmu=xR;
    dem=(pow(xmu,3)-xa)/pow(xm,3)+1.5*(xL*xL-xm*xm)*(xmu*xmu-xb)/(pow(xm,3)*(xm-xL))-3*(xmu-xc)/xm*(1+(xL*xL-xm*xm)/(xm*(xm-xL)));
    localparams[3][1]=(Cmu-D)/dem;
    localparams[3][2]=1.5*localparams[3][1]*(xL*xL-xm*xm)/(xm*(xm-xL));
    localparams[3][3]=-3*localparams[3][1]-2*localparams[3][2];
    coste[3]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[3]);
    costeP[3]=cost_P(rho,P,N,localparams[3],D,xm,xL,xR,PL);

    //OP5:extremo xL, CL,CM
    xalfa=xL; Cmu=CL; xmu=xL; Cnu=CM; xnu=xm;
    dem=(xmu*xmu-xb)/(xm*xm)-2*(xmu-xc)*xalfa/(xm*xm);
    ybr=(Cmu-D)/dem;
    yba=-(xmu*xmu*xmu-xa-3*xalfa*xalfa*(xmu-xc))/(pow(xm,3)*dem);
    dem=(xnu*xnu*xnu-xa)/pow(xm,3)+yba*(xnu*xnu-xb)/(xm*xm)-(xnu-xc)*(3*xalfa*xalfa/(xm*xm)+2*yba*xalfa/xm)/xm;
    localparams[4][1]=(Cnu-D-ybr*(xnu*xnu-xb)/(xm*xm)+2*ybr*xalfa*(xnu-xc)/(xm*xm))/dem;
    localparams[4][2]=ybr+yba*localparams[4][1];
    localparams[4][3]=-3*localparams[4][1]*xalfa*xalfa/(xm*xm)-2*localparams[4][2]*xalfa/xm;
    coste[4]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[4]);
    costeP[4]=cost_P(rho,P,N,localparams[4],D,xm,xL,xR,PL);

    //OP6:extremo xL, CL,CR
    xalfa=xL; Cmu=CL; xmu=xL; Cnu=CR; xnu=xR;
    dem=(xmu*xmu-xb)/(xm*xm)-2*(xmu-xc)*xalfa/(xm*xm);
    ybr=(Cmu-D)/dem;
    yba=-(xmu*xmu*xmu-xa-3*xalfa*xalfa*(xmu-xc))/(pow(xm,3)*dem);
    dem=(xnu*xnu*xnu-xa)/pow(xm,3)+yba*(xnu*xnu-xb)/(xm*xm)-(xnu-xc)*(3*xalfa*xalfa/(xm*xm)+2*yba*xalfa/xm)/xm;
    localparams[5][1]=(Cnu-D-ybr*(xnu*xnu-xb)/(xm*xm)+2*ybr*xalfa*(xnu-xc)/(xm*xm))/dem;
    localparams[5][2]=ybr+yba*localparams[5][1];
    localparams[5][3]=-3*localparams[5][1]*xalfa*xalfa/(xm*xm)-2*localparams[5][2]*xalfa/xm;
    coste[5]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[5]);
    costeP[5]=cost_P(rho,P,N,localparams[5],D,xm,xL,xR,PL);

    //OP7:extremo xL, CM, CR
    xalfa=xL; Cmu=CM; xmu=xm; Cnu=CR; xnu=xR;
    dem=(xmu*xmu-xb)/(xm*xm)-2*(xmu-xc)*xalfa/(xm*xm);
    ybr=(Cmu-D)/dem;
    yba=-(xmu*xmu*xmu-xa-3*xalfa*xalfa*(xmu-xc))/(pow(xm,3)*dem);
    dem=(xnu*xnu*xnu-xa)/pow(xm,3)+yba*(xnu*xnu-xb)/(xm*xm)-(xnu-xc)*(3*xalfa*xalfa/(xm*xm)+2*yba*xalfa/xm)/xm;
    localparams[6][1]=(Cnu-D-ybr*(xnu*xnu-xb)/(xm*xm)+2*ybr*xalfa*(xnu-xc)/(xm*xm))/dem;
    localparams[6][2]=ybr+yba*localparams[6][1];
    localparams[6][3]=-3*localparams[6][1]*xalfa*xalfa/(xm*xm)-2*localparams[6][2]*xalfa/xm;
    coste[6]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[6]);
    costeP[6]=cost_P(rho,P,N,localparams[6],D,xm,xL,xR,PL);

    //OP8:extremo xm, CL,CM
    xalfa=xm; Cmu=CL; xmu=xL; Cnu=CM; xnu=xm;
    dem=(xmu*xmu-xb)/(xm*xm)-2*(xmu-xc)*xalfa/(xm*xm);
    ybr=(Cmu-D)/dem;
    yba=-(xmu*xmu*xmu-xa-3*xalfa*xalfa*(xmu-xc))/(pow(xm,3)*dem);
    dem=(xnu*xnu*xnu-xa)/pow(xm,3)+yba*(xnu*xnu-xb)/(xm*xm)-(xnu-xc)*(3*xalfa*xalfa/(xm*xm)+2*yba*xalfa/xm)/xm;
    localparams[7][1]=(Cnu-D-ybr*(xnu*xnu-xb)/(xm*xm)+2*ybr*xalfa*(xnu-xc)/(xm*xm))/dem;
    localparams[7][2]=ybr+yba*localparams[7][1];
    localparams[7][3]=-3*localparams[7][1]*xalfa*xalfa/(xm*xm)-2*localparams[7][2]*xalfa/xm;
    coste[7]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[7]);
    costeP[7]=cost_P(rho,P,N,localparams[7],D,xm,xL,xR,PL);

    //OP9:extremo xm, CL,CR
    xalfa=xm; Cmu=CL; xmu=xL; Cnu=CR; xnu=xR;
    dem=(xmu*xmu-xb)/(xm*xm)-2*(xmu-xc)*xalfa/(xm*xm);
    ybr=(Cmu-D)/dem;
    yba=-(xmu*xmu*xmu-xa-3*xalfa*xalfa*(xmu-xc))/(pow(xm,3)*dem);
    dem=(xnu*xnu*xnu-xa)/pow(xm,3)+yba*(xnu*xnu-xb)/(xm*xm)-(xnu-xc)*(3*xalfa*xalfa/(xm*xm)+2*yba*xalfa/xm)/xm;
    localparams[8][1]=(Cnu-D-ybr*(xnu*xnu-xb)/(xm*xm)+2*ybr*xalfa*(xnu-xc)/(xm*xm))/dem;
    localparams[8][2]=ybr+yba*localparams[8][1];
    localparams[8][3]=-3*localparams[8][1]*xalfa*xalfa/(xm*xm)-2*localparams[8][2]*xalfa/xm;
    coste[8]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[8]);
    costeP[8]=cost_P(rho,P,N,localparams[8],D,xm,xL,xR,PL);

    //OP10:extremo xm, CM, CR
    xalfa=xm; Cmu=CM; xmu=xm; Cnu=CR; xnu=xR;
    dem=(xmu*xmu-xb)/(xm*xm)-2*(xmu-xc)*xalfa/(xm*xm);
    ybr=(Cmu-D)/dem;
    yba=-(xmu*xmu*xmu-xa-3*xalfa*xalfa*(xmu-xc))/(pow(xm,3)*dem);
    dem=(xnu*xnu*xnu-xa)/pow(xm,3)+yba*(xnu*xnu-xb)/(xm*xm)-(xnu-xc)*(3*xalfa*xalfa/(xm*xm)+2*yba*xalfa/xm)/xm;
    localparams[9][1]=(Cnu-D-ybr*(xnu*xnu-xb)/(xm*xm)+2*ybr*xalfa*(xnu-xc)/(xm*xm))/dem;
    localparams[9][2]=ybr+yba*localparams[9][1];
    localparams[9][3]=-3*localparams[9][1]*xalfa*xalfa/(xm*xm)-2*localparams[9][2]*xalfa/xm;
    coste[9]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[9]);
    costeP[9]=cost_P(rho,P,N,localparams[9],D,xm,xL,xR,PL);

    //removing models with shape that is not suitable for the PT, and physically invalid
    double det,r1,r2;
    for(int model=0; model<10; model++){
        det=localparams[model][2]*localparams[model][2]-3*localparams[model][1]*localparams[model][3];
        if(det<0){ //no extremes in the model, we need a minimum
            coste[model]=1e25;
            costeP[model]=1e25;
        }
        else{
            r1=xm*(-localparams[model][2]+sqrt(det))/(3*localparams[model][1]);
            r2=xm*(-localparams[model][2]-sqrt(det))/(3*localparams[model][1]);
            if(  ((r1<xL || r1>xR) && (r2<xL || r2>xR))  ||  ((r1>xL && r1<xR) && (r2>xL && r2<xR)) ){
                //both roots inside or both roots outside
                coste[model]=1e25;
                costeP[model]=1e25;
            }
            else{ //one root inside
                if(r1>xL && r1<xR){ //r1 is the one inside
                    if ((6*r1*localparams[model][1]/pow(xm,3)+2*localparams[model][2]/pow(xm,2))<0){ //the root inside is a maximum
                        coste[model]=1e25;
                        costeP[model]=1e25;
                    }
                    else{ //this is the minimum. Check thermo consistency and causality
                        if(c_fit(r1,xm,PR,PL,xL,xR,localparams[model])<0){
                            coste[model]=1e30;
                            costeP[model]=1e30;
                        }
                    }
                }
                else{ //r2 is the one inside
                    if ((6*r2*localparams[model][1]/pow(xm,3)+2*localparams[model][2]/pow(xm,2))<0){ //the root inside is a maximum
                        coste[model]=1e25;
                        costeP[model]=1e25;
                    }
                    else{ //this is the minimum. Check thermo consistency and causality
                        if(c_fit(r2,xm,PR,PL,xL,xR,localparams[model])<0){
                            coste[model]=1e30;
                            costeP[model]=1e30;
                        }
                    }
                }
            }

        }
    }
    
    double costeMedio[10]; for(int k=0;k<10;k++) costeMedio[k]=(coste[k]+costeP[k])*0.5;
    int this=minvector(10,costeMedio);    

    params[0]=0.0;
    params[1]=localparams[this][1];
    params[2]=localparams[this][2];
    params[3]=localparams[this][3];

    return costeMedio[this];

}   

double getabc_quartic(double *params, double xL,double xR,double PL,double PR,double xm,double CL,double CR,double CM,
             int N, double *rho, double *csc, double *P){
    double D, alfa, beta, xi, delta; double coste[5],costeP[5];
    D=(PR-PL)/(xR-xL); 
    alfa=(pow(xR,5)-pow(xL,5))/(5*(xR-xL));
    beta=(pow(xR,4)-pow(xL,4))/(4*(xR-xL));
    xi=(pow(xR,3)-pow(xL,3))/(3*(xR-xL));
    delta=(pow(xR,2)-pow(xL,2))/(2*(xR-xL));

    int model;

    double **localparams;
    localparams=(double **)malloc(10*sizeof(double*));
    for(int j=0; j<5;j++){
        localparams[j]=(double *)malloc(4*sizeof(double));
    }


    double yfr,yfc,yfb,yfa, ycr,ycb,yca, ybr,yba, yar, dem;
    double x1, x2, x3, x4, c1, c2, c4;

    //OP1: CL, CM, CR, c'(xL)=0
    x1=xL; c1=CL; x2=xm; c2=CM; x3=xL; x4=xR; c4=CR; model=0;

    dem=xm/(x1-delta);
    yfr=(c1-D)*dem;
    yfa=(pow(x1,4)-alfa)/pow(xm,4)*dem;
    yfb=(pow(x1,3)-beta)/pow(xm,3)*dem;
    yfc=(pow(x1,2)-xi)/pow(xm,2)*dem;

    dem=(pow(x2,2)-xi)/pow(xm,2)-yfc*(x2-delta)/xm;
    ycr=(c2-D-yfr*(x2-delta)/xm)/dem;
    yca=((pow(x2,4)-alfa)/pow(xm,4)-yfa*(x2-delta)/xm)/dem;
    ycb=((pow(x2,3)-beta)/pow(xm,3)-yfb*(x2-delta)/xm)/dem;

    dem=3*xm*x3*x3-2*x3*xm*xm*ycb-xm*xm*xm*yfb+xm*xm*xm*yfc*ycb;
    ybr=(2*x3*xm*xm*ycr+xm*xm*xm*yfr-yfc*ycr)/dem;
    yba=(4*x3*x3*x3-2*x3*xm*xm*yca-xm*xm*xm*yfb+xm*xm*xm*yfc*ycb)/dem;
    //-----all the same until here----
    dem=(pow(x4,4)-alfa)/pow(xm,4)-yba*(pow(x4,3)-beta)/pow(xm,3)+(pow(x4,2)-xi)/pow(xm,2)*(-yca+ycb*yba)+(x4-delta)/xm*(-yfa+yfb*yba+yfc*yca-yfc*ycb*yba);
    yar=c4-D+(pow(x4,3)-beta)/pow(xm,3)*ybr+(pow(x4,2)-xi)/pow(xm,2)*(-ycr-ycb*ybr)+(x4-delta)/xm*(-yfr-yfb*ybr+yfc*ycr+yfc*ycb*ybr);
    
    localparams[model][0]=yar/dem;
    localparams[model][1]=-ybr-yba*localparams[model][0];
    localparams[model][2]=ycr-yca*localparams[model][0]-ycb*localparams[model][1];
    localparams[model][3]=yfr-yfa*localparams[model][0]-yfb*localparams[model][1]-yfc*localparams[model][2];
   coste[model]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[model]);
    costeP[model]=cost_P(rho,P,N,localparams[model],D,xm,xL,xR,PL);

    //OP2: CL, CM, CR, c'(xL)=0
    x1=xL; c1=CL; x2=xm; c2=CM; x3=xm; x4=xR; c4=CR; model=1;

    dem=xm/(x1-delta);
    yfr=(c1-D)*dem;
    yfa=(pow(x1,4)-alfa)/pow(xm,4)*dem;
    yfb=(pow(x1,3)-beta)/pow(xm,3)*dem;
    yfc=(pow(x1,2)-xi)/pow(xm,2)*dem;

    dem=(pow(x2,2)-xi)/pow(xm,2)-yfc*(x2-delta)/xm;
    ycr=(c2-D-yfr*(x2-delta)/xm)/dem;
    yca=((pow(x2,4)-alfa)/pow(xm,4)-yfa*(x2-delta)/xm)/dem;
    ycb=((pow(x2,3)-beta)/pow(xm,3)-yfb*(x2-delta)/xm)/dem;

    dem=3*xm*x3*x3-2*x3*xm*xm*ycb-xm*xm*xm*yfb+xm*xm*xm*yfc*ycb;
    ybr=(2*x3*xm*xm*ycr+xm*xm*xm*yfr-yfc*ycr)/dem;
    yba=(4*x3*x3*x3-2*x3*xm*xm*yca-xm*xm*xm*yfb+xm*xm*xm*yfc*ycb)/dem;
    //-----all the same until here----
    dem=(pow(x4,4)-alfa)/pow(xm,4)-yba*(pow(x4,3)-beta)/pow(xm,3)+(pow(x4,2)-xi)/pow(xm,2)*(-yca+ycb*yba)+(x4-delta)/xm*(-yfa+yfb*yba+yfc*yca-yfc*ycb*yba);
    yar=c4-D+(pow(x4,3)-beta)/pow(xm,3)*ybr+(pow(x4,2)-xi)/pow(xm,2)*(-ycr-ycb*ybr)+(x4-delta)/xm*(-yfr-yfb*ybr+yfc*ycr+yfc*ycb*ybr);
    
    localparams[model][0]=yar/dem;
    localparams[model][1]=-ybr-yba*localparams[model][0];
    localparams[model][2]=ycr-yca*localparams[model][0]-ycb*localparams[model][1];
    localparams[model][3]=yfr-yfa*localparams[model][0]-yfb*localparams[model][1]-yfc*localparams[model][2];
   coste[model]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[model]);
    costeP[model]=cost_P(rho,P,N,localparams[model],D,xm,xL,xR,PL);

    //OP3:  CL, CM, extremos
    x1=xL; c1=CL; x2=xm; c2=CM; x3=xL; x4=xm; model=2;

    dem=xm/(x1-delta);
    yfr=(c1-D)*dem;
    yfa=(pow(x1,4)-alfa)/pow(xm,4)*dem;
    yfb=(pow(x1,3)-beta)/pow(xm,3)*dem;
    yfc=(pow(x1,2)-xi)/pow(xm,2)*dem;

    dem=(pow(x2,2)-xi)/pow(xm,2)-yfc*(x2-delta)/xm;
    ycr=(c2-D-yfr*(x2-delta)/xm)/dem;
    yca=((pow(x2,4)-alfa)/pow(xm,4)-yfa*(x2-delta)/xm)/dem;
    ycb=((pow(x2,3)-beta)/pow(xm,3)-yfb*(x2-delta)/xm)/dem;

    dem=3*xm*x3*x3-2*x3*xm*xm*ycb-xm*xm*xm*yfb+xm*xm*xm*yfc*ycb;
    ybr=(2*x3*xm*xm*ycr+xm*xm*xm*yfr-yfc*ycr)/dem;
    yba=(4*x3*x3*x3-2*x3*xm*xm*yca-xm*xm*xm*yfb+xm*xm*xm*yfc*ycb)/dem;
    //-----all the same until here----
    dem=4*x4*x4*x4-3*xm*x4*x4*yba+2*x4*xm*xm*(-yca+ycb*yba)+xm*xm*xm*(-yfa+yfb*yba+yfc*yca-yfc*ycb*yba);
    yar=3*xm*x4*x4*ybr+2*x4*xm*xm*(-ycr-ycb*ybr)+xm*xm*xm*(-yfr-yfb*ybr+yfc*ycr+yfc*ycb*ybr);
    
    localparams[model][0]=yar/dem;
    localparams[model][1]=-ybr-yba*localparams[model][0];
    localparams[model][2]=ycr-yca*localparams[model][0]-ycb*localparams[model][1];
    localparams[model][3]=yfr-yfa*localparams[model][0]-yfb*localparams[model][1]-yfc*localparams[model][2];
    coste[model]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[model]);
    costeP[model]=cost_P(rho,P,N,localparams[model],D,xm,xL,xR,PL);

    //OP4: CL, CR, extremos
    x1=xL; c1=CL; x2=xR; c2=CR; x3=xL; x4=xm; model=3;

    dem=xm/(x1-delta);
    yfr=(c1-D)*dem;
    yfa=(pow(x1,4)-alfa)/pow(xm,4)*dem;
    yfb=(pow(x1,3)-beta)/pow(xm,3)*dem;
    yfc=(pow(x1,2)-xi)/pow(xm,2)*dem;

    dem=(pow(x2,2)-xi)/pow(xm,2)-yfc*(x2-delta)/xm;
    ycr=(c2-D-yfr*(x2-delta)/xm)/dem;
    yca=((pow(x2,4)-alfa)/pow(xm,4)-yfa*(x2-delta)/xm)/dem;
    ycb=((pow(x2,3)-beta)/pow(xm,3)-yfb*(x2-delta)/xm)/dem;

    dem=3*xm*x3*x3-2*x3*xm*xm*ycb-xm*xm*xm*yfb+xm*xm*xm*yfc*ycb;
    ybr=(2*x3*xm*xm*ycr+xm*xm*xm*yfr-yfc*ycr)/dem;
    yba=(4*x3*x3*x3-2*x3*xm*xm*yca-xm*xm*xm*yfb+xm*xm*xm*yfc*ycb)/dem;
    //-----all the same until here----
    dem=4*x4*x4*x4-3*xm*x4*x4*yba+2*x4*xm*xm*(-yca+ycb*yba)+xm*xm*xm*(-yfa+yfb*yba+yfc*yca-yfc*ycb*yba);
    yar=3*xm*x4*x4*ybr+2*x4*xm*xm*(-ycr-ycb*ybr)+xm*xm*xm*(-yfr-yfb*ybr+yfc*ycr+yfc*ycb*ybr);
    localparams[model][0]=yar/dem;
    localparams[model][1]=-ybr-yba*localparams[model][0];
    localparams[model][2]=ycr-yca*localparams[model][0]-ycb*localparams[model][1];
    localparams[model][3]=yfr-yfa*localparams[model][0]-yfb*localparams[model][1]-yfc*localparams[model][2];
   coste[model]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[model]);
    costeP[model]=cost_P(rho,P,N,localparams[model],D,xm,xL,xR,PL);

    //OP5: CM, CR, extremos
    x1=xm; c1=CM; x2=xR; c2=CR; x3=xL; x4=xm; model=4;

    dem=xm/(x1-delta);
    yfr=(c1-D)*dem;
    yfa=(pow(x1,4)-alfa)/pow(xm,4)*dem;
    yfb=(pow(x1,3)-beta)/pow(xm,3)*dem;
    yfc=(pow(x1,2)-xi)/pow(xm,2)*dem;

    dem=(pow(x2,2)-xi)/pow(xm,2)-yfc*(x2-delta)/xm;
    ycr=(c2-D-yfr*(x2-delta)/xm)/dem;
    yca=((pow(x2,4)-alfa)/pow(xm,4)-yfa*(x2-delta)/xm)/dem;
    ycb=((pow(x2,3)-beta)/pow(xm,3)-yfb*(x2-delta)/xm)/dem;

    dem=3*xm*x3*x3-2*x3*xm*xm*ycb-xm*xm*xm*yfb+xm*xm*xm*yfc*ycb;
    ybr=(2*x3*xm*xm*ycr+xm*xm*xm*yfr-yfc*ycr)/dem;
    yba=(4*x3*x3*x3-2*x3*xm*xm*yca-xm*xm*xm*yfb+xm*xm*xm*yfc*ycb)/dem;
    //-----all the same until here----
    dem=4*x4*x4*x4-3*xm*x4*x4*yba+2*x4*xm*xm*(-yca+ycb*yba)+xm*xm*xm*(-yfa+yfb*yba+yfc*yca-yfc*ycb*yba);
    yar=3*xm*x4*x4*ybr+2*x4*xm*xm*(-ycr-ycb*ybr)+xm*xm*xm*(-yfr-yfb*ybr+yfc*ycr+yfc*ycb*ybr);
    localparams[model][0]=yar/dem;
    localparams[model][1]=-ybr-yba*localparams[model][0];
    localparams[model][2]=ycr-yca*localparams[model][0]-ycb*localparams[model][1];
    localparams[model][3]=yfr-yfa*localparams[model][0]-yfb*localparams[model][1]-yfc*localparams[model][2];
    coste[model]=cost(rho,csc,N,xm,PR,PL,xL,xR,localparams[model]);
    costeP[model]=cost_P(rho,P,N,localparams[model],D,xm,xL,xR,PL);

    //check appropiate shape for the model, thermo consistency and causality
    double det, theta, Q,R,S1,S2, a1,a2,a3, r1,r2,r3;
    for(int model=0; model<5; model++){
        a1=3*localparams[model][1]*xm/(4*localparams[model][0]);
        a2=localparams[model][2]*pow(xm,2)/(2*localparams[model][0]);
        a3=localparams[model][3]*pow(xm,3)/(4*localparams[model][0]);

        Q=(3*a2-a1*a1)/9; R=(9*a1*a2-27*a3-2*a1*a1*a1)/54;
        det=Q*Q*Q+R*R;


        if(det>0){ //two complex roots. Only interested in the real one
            S1=cbrt(R+sqrt(det));
            S2=cbrt(R-sqrt(det));
            r1=S1+S2-a1/3;
            if(r1>xL && r1<xR){ //r1 is inside
                if ((12*r1*r1*localparams[model][0]/pow(xm,4)+6*r1*localparams[model][1]/pow(xm,3)+2*localparams[model][2]/pow(xm,2))<0){ //the root inside is a maximum
                    coste[model]=1e25;
                    costeP[model]=1e25;
                }
                else{ //this is the minimum. Check thermo consistency and causality
                    if(c_fit(r1,xm,PR,PL,xL,xR,localparams[model])<0){
                        coste[model]=1e30;
                        costeP[model]=1e30;
                    }
                }
            }
        }
        else{ //three real roots
            theta=acos(-R/sqrt(-Q*Q*Q));
            r1=-2*sqrt(-Q)*cos(theta/3)-a1/3;
            r2=-2*sqrt(-Q)*cos((theta+2*M_PI)/3)-a1/3;
            r3=-2*sqrt(-Q)*cos((theta+4*M_PI)/3)-a1/3;

           // printf("    3 roots: %e %e and %e\n",r1,r2,r3);

            if((r1>xL && r1<xR) && (r2<xL || r2>xR) && (r3<xL || r3>xR)){ //r1 is inside, the other two outside
                if ((12*r1*r1*localparams[model][0]/pow(xm,4)+6*r1*localparams[model][1]/pow(xm,3)+2*localparams[model][2]/pow(xm,2))<0){ //the root inside is a maximum
                    coste[model]=1e25;
                    costeP[model]=1e25;
                }
                else{ //this is the minimum. Check thermo consistency and causality
                    if(c_fit(r1,xm,PR,PL,xL,xR,localparams[model])<0){
                        coste[model]=1e30;
                        costeP[model]=1e30;
                    }
                }
            }
            else if((r2>xL && r2<xR) && (r1<xL || r1>xR) && (r3<xL || r3>xR)){ //r2 is inside, the other two outside
                if ((12*r2*r2*localparams[model][0]/pow(xm,4)+6*r2*localparams[model][1]/pow(xm,3)+2*localparams[model][2]/pow(xm,2))<0){ //the root inside is a maximum
                    coste[model]=1e25;
                    costeP[model]=1e25;
                }
                else{ //this is the minimum. Check thermo consistency and causality
                    if(c_fit(r2,xm,PR,PL,xL,xR,localparams[model])<0){
                        coste[model]=1e30;
                        costeP[model]=1e30;
                    }
                }
            }
            else if((r3>xL && r3<xR) && (r2<xL || r2>xR) && (r1<xL || r1>xR)){ //r3 is inside, the other two outside
                if ((12*r3*r3*localparams[model][0]/pow(xm,4)+6*r3*localparams[model][1]/pow(xm,3)+2*localparams[model][2]/pow(xm,2))<0){ //the root inside is a maximum
                    coste[model]=1e25;
                    costeP[model]=1e25;
                }
                else{ //this is the minimum. Check thermo consistency and causality
                    if(c_fit(r3,xm,PR,PL,xL,xR,localparams[model])<0){
                        coste[model]=1e30;
                        costeP[model]=1e30;
                    }
                }
            }
            else{ //more than 1 root inside
                coste[model]=1e26;
                costeP[model]=1e26;
            }
        }
    }

    double costeMedio[5]; for(int k=0;k<5;k++) costeMedio[k]=(coste[k]+costeP[k])*0.5;
    int this=minvector(5,costeMedio);
    
    params[0]=localparams[this][0];
    params[1]=localparams[this][1];
    params[2]=localparams[this][2];
    params[3]=localparams[this][3];

    return costeMedio[this];

}   



void fitPolynomials(struct Tab table, int n_PT, double *xL, double *xR, double *Kfitted, double *Gfitted, double *a, double *b, double *c, double *f, double *d, double *xm){
    double left, right, PL, PR, CM, CL, CR, xmpt, dpt;
    double *rhoaux, *cscaux, *Paux, model_cost[3], parabola[4], cubic[4], quartic[4];
    int ptlen, index_min, selected;
    int i;
    for(i=0; i<n_PT; i++){
        left=xL[i]; right=xR[i];
        PL=Kfitted[i]*pow(left, Gfitted[i]);
        PR=Kfitted[i+1]*pow(right, Gfitted[i+1]);
        ptlen=restructure(left, right, table, &rhoaux, &cscaux, &Paux);
        index_min=find_min(cscaux,ptlen);
        xmpt=rhoaux[index_min];
        CL=cscaux[0]; CR=cscaux[ptlen-1]; CM=cscaux[index_min];
        dpt=(PR-PL)/(right-left);

        xm[i]=xmpt;
        d[i]=dpt;

        model_cost[0]=getabc_parabola(parabola, left, right, PL, PR, xmpt, CL, CR, CM, ptlen, rhoaux, cscaux, Paux);
        model_cost[1]=getabc_cubic(cubic, left, right, PL, PR, xmpt, CL, CR, CM, ptlen, rhoaux, cscaux, Paux);
        model_cost[2]=getabc_quartic(quartic, left, right, PL, PR, xmpt, CL, CR, CM, ptlen, rhoaux, cscaux, Paux);
        
        selected=minvector(3,model_cost);
        printf("---Modeling PT number %d\n", i+1);
        printf("    Cost of models: %e %e %e:\n",model_cost[0],model_cost[1],model_cost[2]);
        if(selected==0) {
            printf("    Chosen parabola\n");
            a[i]=parabola[0];
            b[i]=parabola[1];
            c[i]=parabola[2];
            f[i]=parabola[3];
        }
        if(selected==1) {
            printf("    Chosen cubic\n");
            a[i]=cubic[0];
            b[i]=cubic[1];
            c[i]=cubic[2];
            f[i]=cubic[3];
        }
        if(selected==2) {
            printf("    Chosen quartic\n");
            a[i]=quartic[0];
            b[i]=quartic[1];
            c[i]=quartic[2];
            f[i]=quartic[3];
        }

    }

}