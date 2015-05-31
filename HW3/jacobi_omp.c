#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"


int main (int argc, char *argv[])
{
    int n=atol(argv[1]); //Mesh size n+2;
    int iMax=atol(argv[2]);//max number of loops;
    
   int i,iter=0,nthread;
   double h=1.0/(n+1), h2=h*h, tol=1e-5,r,res,temp;
   double *u, *u1, *f;
    
    u=calloc(n+2, sizeof(double));
    u1=(double *)malloc(sizeof(double)*(n+2));
    f=(double *)malloc(sizeof(double)*(n+2));
    
    //Initialize:
    for (i=0;i<n+2;i++)
        f[i]=1;
    res=sqrt(n);
   // printf("res0=%e\n", res);
 
    timestamp_type time1, time2;
    get_timestamp(&time1);
    
#pragma omp parallel
{
    #pragma omp master
    {  nthread=omp_get_num_threads();
        printf("Number of threads = %d\n", nthread);
    }
}
  
    while (res>tol && iter<iMax){
#pragma omp parallel for shared(u1,u,f,n,h2) private(i)
        for(i=1;i<n+1;i++){
            u1[i]=f[i]*h2/2+(u[i+1]+u[i-1])/2;
        }
        
#pragma omp parallel for shared(u1,u,n) private(i)
        for(i=1;i<n+1;i++){
            u[i]=u1[i];
        }
        
    //get residue:
    res=0;
    for(i=1;i<n+1;i++){
    #pragma omp critical
        res+=((2*u[i]-u[i+1]-u[i-1])/h2-f[i])*((2*u[i]-u[i+1]-u[i-1])/h2-f[i]);
    }
    res=sqrt(res);
        
    iter++;
    }
    
    get_timestamp(&time2);
    double T=timestamp_diff_in_seconds(time1,time2);


    printf("Residue after %d loops is %e\n",iter,res);
    printf("Time elapsed is %e sec.\n",T);
  


return 0;



}

