#include "findPT.h"

//Find Phase transitions 
int findCSCchanging(int imin,struct Tab table){
    int i;
    double q1,q2,q3,extra;
    if(imin>=table.lenEos-3) return -1;
    q1=table.csctab[imin];
    q2=table.csctab[imin+1];
    q3=table.csctab[imin+2];
    for (i=imin+1; i<table.lenEos-3; i++){
        if((q3-q2)<=0 && (q2-q1)>=0){
            extra=table.csctab[i+2];
            if((extra-q3)*(q3-q2)>=0)
                break;
            else{
                q1=q2;
                q2=q3;
                q3=table.csctab[i+2];
            }
        }
        else{
            q1=q2;
            q2=q3;
            q3=table.csctab[i+2];
        }
    }
    if(i==(table.lenEos-3)) i=-1;
    return i;
}
double TV(int i, double *v, int lenEos){
    double variation=0.0;
    if(i<lenEos-4){
        for (int j=i; j<i+3; j++){
            variation+=fabs(v[j+1]-v[j]);
        }
        return variation;
    }
    else{
        for (int j=i; j<lenEos-1; j++){
            variation+=fabs(v[j+1]-v[j]);
        }
        return variation;
    }
}
int findKinkAdiabaticIndex(int index,struct Tab table){
    double vprevprev, vprev, vhere, vnext, vnextnext;
    if(index+2>table.lenEos-4) return 0;
    vhere=TV(index,table.aditab,table.lenEos);
    vnext=TV(index+1,table.aditab,table.lenEos);
    vprev=TV(index-1,table.aditab,table.lenEos);
    if(vhere<0 || vprev<0 || vnext<0) return 0;

    if (vhere>vnext && vhere>vprev){
        return 2;
    }
    else{//due to oscillations, we allow the peak to be right before or after
        vprevprev=TV(index-2,table.aditab,table.lenEos);
        vnextnext=TV(index+2,table.aditab,table.lenEos);
        if(vprevprev<0 || vnextnext<0) return 0;
        if(vprev>vprevprev && vprev>vhere)
            return 1;
        else if (vnext>vnextnext && vnext>vhere)
            return 1;
        else
            return 0;
    }    
}
int findCSCsmooth(int imin,struct Tab table){
    double tvhere, tvnext;
    tvhere=TV(imin,table.csctab,table.lenEos);
    int j;

    for(j=imin+1; j<table.lenEos-3; j++){
        tvnext=TV(j,table.csctab,table.lenEos);
        if(tvnext<=tvhere && table.csctab[j]>=table.csctab[imin]){
            return j;
        }
    }
   return -1;
}
int checkOscillations(int index, int s, struct Tab table){
    int changes=0; double mon, monprev;
    if (index<10 || index>table.lenEos-10) return 0;  //do not put a PT here, it is at the border
    
    else if(s<0){
        monprev=table.aditab[index-1]-table.aditab[index];
        for (int i=index-1; i>index-10; i--){
            mon=table.aditab[i-1]-table.aditab[i];
            if(mon*monprev<0) changes++;
            monprev=mon;
        }
    }
    else{
        monprev=table.aditab[index+1]-table.aditab[index];
        for (int i=index+1; i>index+10; i++){
            mon=table.aditab[i+1]-table.aditab[i];
            if(mon*monprev<0) changes++;
            monprev=mon;
        }
    }
    if (changes>=3) return 0;  //do not put a PT here: the table oscillates
    else return 1;
}
int findPTs(double *xL, double *xR, struct Tab table){
    int i,ok,index, n_PT=0,k=0;
    for(i=0; i<table.lenEos-4; i++){
        ok=0;
        if ((index=findCSCchanging(i,table))>=0){
            if(findKinkAdiabaticIndex(index,table)){
                ok=1;
            }
            i=index;
        }
        if(ok){
            if((index=findCSCsmooth(index, table))>=0){
                if(checkOscillations(i, -1, table) && checkOscillations(index, 1, table)){
                    if(xL) xL[k]=table.rhotab[i];
                    if(xR) xR[k]=table.rhotab[index];
                    k++;
                    n_PT++;
                    i=index;
                }
            }
        }
    }
    return  n_PT; 
}
int obtainPTs(double **xL, double **xR, struct Tab table){
    int N;
    N=findPTs(NULL, NULL, table);
    (*xL)=(double *)malloc(N*sizeof(double));
    (*xR)=(double *)malloc(N*sizeof(double));
    findPTs(*xL, *xR, table);
    return N;
}

//Find the need of extra cut for polytropes

double sigma(double *v, double m, double N){
    double deviation, sumsqr=0;
    for (int i = 1 ; i<= N; i++){
        deviation = v[i] - m;
        sumsqr += deviation * deviation;
    }
    return sqrt(sumsqr/N) ;
}

double needExtraCut(struct Tab table, double xR){
    int min_index;
    for(int i=0; i<table.lenEos; i++){
        if(table.rhotab[i]>xR){
            min_index=i;
            break;
        }
    }

    double *ltv;
    ltv=(double*)malloc((table.lenEos-min_index+1)*sizeof(double));
    double decreasing_rho=-1;

    for (int i=min_index; i<table.lenEos-1; i++){
        if(table.csctab[i-1]<table.csctab[i] && table.csctab[i]>table.csctab[i+1]){
            if(table.csctab[i+1]>table.csctab[i+2])
               decreasing_rho=table.rhotab[i];
        }
    }
    int count=0; double mean=0.0;
    for (int i=min_index; i<table.lenEos; i++){
        ltv[count]=TV(i, table.csctab, table.lenEos); 
        mean+=ltv[count];
        count++;
    }
    mean=mean/(table.lenEos-min_index+1);
    double s; s=sigma(ltv, mean, table.lenEos-min_index+1);

    int ok1=1, ok2=1; int index; double max=-1;
    for (int i=1; i<table.lenEos-min_index; i++){
        if(ltv[i-1]<ltv[i] && ltv[i]>ltv[i+1]){
            ok1=1; ok2=1;
            if(i>3){
                if((ltv[i-3]-ltv[i-2])*(ltv[i-2]-ltv[i-1])<0){
                    ok1=0;
                }
            }
            else{ok1=0;}
            if(i<table.lenEos-3){
                if((ltv[i+1]-ltv[i+2])*(ltv[i+2]-ltv[i+3])<0){
                    ok2=0;
                } 
            }
            else{ok2=0;}
            if(ok1 && ok2){
                if(ltv[i]>max){
                    max=ltv[i];
                    index=i;
                }
            }
        }  
    }
    free(ltv);
    if(max>(mean+2*s)){
        return table.rhotab[min_index+index];
    }
    else if (decreasing_rho>0){
        return decreasing_rho;
    }
    else return -1;

}














