#include "fitPolytropes.h"

#define MAX_ITER 1000
#define TOL 1e-6

double log_polytrope(double x, double logK, double Gamma){  //x is log(rho)
    return logK+Gamma*x;
}
double cost_log_polytrope(double *x, double *y, double *beta, int N){
    double logk=beta[0], Gamma=beta[1];
    double term, sum=0.0; int points=0;
    for (int i=0; i<N; i++){
        term=y[i]-log_polytrope(x[i], logk, Gamma);
        points++;
        sum+=term*term;
    }
    return sum;
}
double cost_nolog(double *x, double *y, int N, double KL, double GL, double *beta){
    double sum=0.0, term; int points=0;
    for(int i=0; i<N; i++){
        if(x[i]>beta[0]){
            term=log10(y[i])-log10(KL*pow(beta[0],GL-beta[1])*pow(x[i],beta[1]));
            points++;
        }
        else{
            term=0;
        }
        sum+=term*term;    
    }
    return sum;
}


void LM_general(double *x, double *y, int N, double *Kc, double *ga){ //perfecto!
    double beta[2]; double lr=0.8;
    beta[0]=-25; beta[1]=2.5;
    double *A, *B;

    if((B=(double*)malloc(N*sizeof(double)))==NULL){printf("error allocating\n"); exit(1);}

    //prepare inverse matrix
    for (int i=0; i<N; i++){
        //A[i]=1;
        B[i]=x[i];
    }

        //all components
    double a,b,c,d;
    double sum=0.0;
    a=N;
    sum=0.0; for(int i=0; i<N; i++) sum+=1*B[i];
    b=sum;
    c=b;
    sum=0.0; for(int i=0; i<N; i++) sum+=B[i]*B[i];
    d=sum;

    double det=a*d-b*c; det=1.0/det;
    double dbeta[2];
    //Delta y
    double *Delta;
    if((Delta=(double*)malloc(N*sizeof(double)))==NULL){printf("error allocating\n"); exit(1);}
    double TA,TB;
    
    double new_cost=1, last_cost=1;
    
    new_cost=cost_log_polytrope(x, y, beta, N);
    int counter=0; double old_beta0=1;
    while(fabs(old_beta0-beta[1])>1e-13){//fabs(new_cost-last_cost)>1e-10){
        counter++;
        last_cost=new_cost;
        old_beta0=beta[1];
        for(int i=0; i<N; i++) Delta[i]=y[i]-log_polytrope(x[i],beta[0],beta[1]);

        sum=0.0; for(int i=0; i<N; i++) sum+=Delta[i];
        TA=sum;
        sum=0.0; for(int i=0; i<N; i++) sum+=B[i]*Delta[i];
        TB=sum;
        dbeta[0]=(TA*(d)+TB*(-b))*det;
        dbeta[1]=(TA*(-c)+TB*(a))*det;

        beta[0]+=lr*dbeta[0]; beta[1]+=lr*dbeta[1]; 
        
        new_cost=cost_log_polytrope(x, y, beta, N);
        if (!(counter%10000))
            printf("%e %e  -> %.15e\n",beta[0],beta[1],new_cost);
    }
    printf("got it: log(kappa)=%f Gamma=%f  -> error=%e\n",beta[0],beta[1],new_cost, counter);
    free(B); free(Delta);
    
    (*Kc)=pow(10, beta[0]); (*ga)=beta[1];

}
void LM_gamma(double *x, double *y, int N, double GL, double KL, double rc, double *K, double *G){
    double beta[1]; double lr=1.0;
    beta[0]=2.5;
    double aux[2];
    double *B;

    if((B=(double*)malloc(N*sizeof(double)))==NULL){printf("error allocating\n"); exit(1);}

    //prepare inverse matrix
    for (int i=0; i<N; i++){
        B[i]=x[i]-log10(rc);
    }

            //all components
        double a,b,c,d;
        double sum=0.0;
        sum=0.0; for(int i=0; i<N; i++) sum+=B[i]*B[i];
        d=sum;
        double dbeta[1];
        //Delta y
        double *Delta;
        if((Delta=(double*)malloc(N*sizeof(double)))==NULL){printf("error allocating\n"); exit(1);}
        double TA,TB;
        
        double new_cost=1, last_cost=1;
        aux[0]=log10(KL)+(GL-beta[0])*log10(rc); aux[1]=beta[0];
        new_cost=cost_log_polytrope(x, y, aux, N);
        int counter=0; double old_beta0=1;
        while(fabs(old_beta0-beta[0])>1e-10){
            counter++;
            last_cost=new_cost;
            old_beta0=beta[0];
            for(int i=0; i<N; i++) Delta[i]=y[i]-log_polytrope(x[i],log10(KL)+(GL-beta[0])*log10(rc),beta[0]);
        
            sum=0.0; for(int i=0; i<N; i++) sum+=B[i]*Delta[i];
            TA=sum;
            dbeta[0]=TA/(d); 
            beta[0]+=lr*dbeta[0]; //beta[1]+=lr*dbeta[1]; 
            aux[0]=log10(KL)+(GL-beta[0])*log10(rc); aux[1]=beta[0];
            new_cost=cost_log_polytrope(x, y, aux, N);
        }
        aux[0]=log10(KL)+(GL-beta[0])*log10(rc); aux[1]=beta[0];
        printf("got it: Gamma=%f  -> error=%e\n",beta[0],cost_log_polytrope(x, y, aux, N));
        (*K)=pow(10,aux[0]); (*G)=aux[1];
        free(B); free(Delta);
    
}
void LM_notlog(double *x, double *y, int N, double GL, double KL, double *rc, double *Kc, double *ga ){
    double beta[2]; double lr=0.8;
    beta[0]=1e13; beta[1]=2.5;

    double *A, *B;

    if((A=(double*)malloc(N*sizeof(double)))==NULL){printf("error allocating\n"); exit(1);}
    if((B=(double*)malloc(N*sizeof(double)))==NULL){printf("error allocating\n"); exit(1);}

    //prepare inverse matrix
    for (int i=0; i<N; i++){
        A[i]=KL*(GL-beta[1])*pow(x[i], beta[1])*pow(beta[0], GL-beta[1]-1);
        B[i]=KL*pow(x[i], beta[1])*pow(beta[0], GL-beta[1])*(log(x[i]/beta[0]));
    }

        //all components
    double a,b,c,d;
    double sum=0.0;
    for(int i=0; i<N; i++) sum+=A[i]*A[i];
    a=sum;
    sum=0.0; for(int i=0; i<N; i++) sum+=A[i]*B[i];
    b=sum;
    c=b;
    sum=0.0; for(int i=0; i<N; i++) sum+=B[i]*B[i];
    d=sum;
    
    double det=a*d-b*c; det=1.0/det;

    double dbeta[2];
    //Delta y
    double *Delta;
    if((Delta=(double*)malloc(N*sizeof(double)))==NULL){printf("error allocating\n"); exit(1);}
    double Deltaxl, Deltaxr, Deltaxm;
    double TA,TB;
    
    double new_cost=1, last_cost=1;
    
    new_cost=cost_nolog(x, y, N, KL, GL,beta);
    int counter=0; double old_beta0=1;
    while(fabs(old_beta0-beta[1])>1e-10){
        counter++;
        last_cost=new_cost;
        old_beta0=beta[1];
        for(int i=0; i<N; i++) Delta[i]=y[i]-KL*pow(beta[0],GL-beta[1])*pow(x[i], beta[1]);
    
        sum=0.0; for(int i=0; i<N; i++) sum+=A[i]*Delta[i];
        TA=sum;
        sum=0.0; for(int i=0; i<N; i++) sum+=B[i]*Delta[i];
        TB=sum;
        dbeta[0]=(TA*(d)+TB*(-b))*det; 
        dbeta[1]=(TA*(-c)+TB*(a))*det;
        beta[0]+=lr*dbeta[0]; beta[1]+=lr*dbeta[1]; 
        
        new_cost=cost_nolog(x, y, N, KL, GL,beta);;

        if (!(counter%10000))
            printf("%e %e  -> %.15e\n",beta[0],beta[1],new_cost);
    }
    printf("got it: rc=%e Gamma=%f  -> error=%e\n",beta[0],beta[1],new_cost);

    (*rc)= beta[0]; (*ga)=beta[1]; (*Kc)=KL*pow(beta[0], GL-beta[1]);

    free(A); free(B); free(Delta);
}


int fillAuxiliar(double **rho, double **P, double minrho, double maxrho, struct Tab table){
    int min_index, max_index;
    for(int i=0; i<table.lenEos; i++){
    if(table.rhotab[i]>minrho){
        min_index=i;
        break;
    }
    } max_index=table.lenEos;
    for(int i=min_index; i<table.lenEos; i++){
        if(table.rhotab[i]>=maxrho){
            max_index=i-1;
            break;
        }
    }
    (*P)=(double*)malloc((max_index-min_index+1)*sizeof(double));
    (*rho)=(double*)malloc((max_index-min_index+1)*sizeof(double));
    int count=0;
    for(int i=min_index; i<max_index+1; i++){
        (*P)[count]=table.Ptab[i];
        (*rho)[count]=table.rhotab[i];
        count++;
    }
    return max_index-min_index+1;
}
int fillAuxiliarLog(double **rho, double **P, double minrho, double maxrho, struct Tab table){
    int min_index, max_index;
    for(int i=0; i<table.lenEos; i++){
    if(table.rhotab[i]>minrho){
        min_index=i;
        break;
    }
    }
    max_index=table.lenEos;
    for(int i=min_index; i<table.lenEos; i++){
        if(table.rhotab[i]>=maxrho){
            max_index=i-1;
            break;
        }
    }
    (*P)=(double*)malloc((max_index-min_index+1)*sizeof(double));
    (*rho)=(double*)malloc((max_index-min_index+1)*sizeof(double));
    int count=0;
    for(int i=min_index; i<max_index+1; i++){
        (*P)[count]=log10(table.Ptab[i]);
        (*rho)[count]=log10(table.rhotab[i]);
        count++;
    }
    return max_index-min_index+1;
}


double fitPolytropes(struct Tab table, int n_PT, double *xL, double *xR, double extra_cut, double *Kfitted, double *Gfitted){
    double *aux_rho, *aux_P; int len_aux;
    double rc;
    double crust_start=2.62780e12;
    //first polytrope between crust and first PT
    //last piece of crust starts at : 2.62780e12, and the polytrope is kappa=3.99874e-8, Gamma=1.35692
    len_aux=fillAuxiliar(&aux_rho, &aux_P, crust_start, xL[0], table);
    LM_notlog(aux_rho, aux_P, len_aux, 1.35692, 3.99874e-8, &rc, &Kfitted[0], &Gfitted[0]);
    free(aux_P); free(aux_rho);
    int last_index_poly=0;

    for(int j=0; j<n_PT-1; j++){
        len_aux=fillAuxiliarLog(&aux_rho, &aux_P, xR[j], xL[j+1], table);
        LM_general(aux_rho, aux_P, len_aux, &Kfitted[last_index_poly+1], &Gfitted[last_index_poly+1]); 
        last_index_poly++;
        free(aux_P); free(aux_rho);
    }

    //extra cut
    if(extra_cut>0){
        len_aux=fillAuxiliarLog(&aux_rho, &aux_P, xR[n_PT-1], extra_cut, table);
        LM_general(aux_rho, aux_P, len_aux, &Kfitted[last_index_poly+1], &Gfitted[last_index_poly+1]); 
        last_index_poly++;
        free(aux_P); free(aux_rho); 
        len_aux=fillAuxiliarLog(&aux_rho, &aux_P, extra_cut, 6e15 , table);
        LM_gamma(aux_rho, aux_P, len_aux, Gfitted[last_index_poly], Kfitted[last_index_poly], extra_cut, &Kfitted[last_index_poly+1], &Gfitted[last_index_poly+1]); 
        free(aux_P); free(aux_rho);
    }
    else{
        len_aux=fillAuxiliarLog(&aux_rho, &aux_P, xR[n_PT-1], table.rhotab[table.lenEos-1], table);
        LM_general(aux_rho, aux_P, len_aux, &Kfitted[last_index_poly+1], &Gfitted[last_index_poly+1]); 
        last_index_poly++;
        free(aux_P); free(aux_rho);
    }

    return rc;
}

 