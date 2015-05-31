#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"


int main (int argc, char *argv[])
{  int n=atol(argv[1]); //Mesh size n+2;
   int iMax=atol(argv[2]);//max number of loops;
    
   int i,iter=0,nthread;
   double h=1.0/(n+1), h2=h*h, tol=1e-5,r,res,temp;
   double *u, *f;
    
    u=calloc(n+2, sizeof(double));
    f=(double *)malloc(sizeof(double)*(n+2));
    
    //Initialize:
    for (i=0;i<n+2;i++)
        f[i]=1;
    res=sqrt(n);
    
    timestamp_type time1, time2;
    get_timestamp(&time1);
    
#pragma omp parallel
{
    #pragma omp master
    {  nthread=omp_get_num_threads();
        //printf("Number of threads = %d\n", nthread);
    }
}
    
    while (res>tol && iter<iMax){
        
        //odd swipe:
#pragma omp parallel for shared(u,f,n,h2) private(i)
        for(i=1;i<n+1;i+=2){
            u[i]=f[i]*h2/2+(u[i+1]+u[i-1])/2;
        }
        
        //even swipe:
#pragma omp parallel for shared(u,f,n,h2) private(i)
        for(i=2;i<n+1;i+=2){
            u[i]=f[i]*h2/2+(u[i+1]+u[i-1])/2;
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


    //printf("Residue after %d loops is %e\n",iter,res);
    //printf("Time elapsed is %e sec.\n",T);
    
    //print to file:
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename,256,"gs_output%02d.txt",nthread);
    fd=fopen(filename,"w+");
    
    if(NULL == fd){
        printf("Error opening file!\n");
        return 1;
    }
    
    fprintf(fd, "Mesh size n=%d\n", n);
    fprintf(fd,"Number of threads=%d\n",nthread);
    fprintf(fd,"Residue after %d loops is %e\n",iter,res);
    fprintf(fd,"Time elapsed is %e sec.\n",T);
    fclose(fd);

return 0;



}

