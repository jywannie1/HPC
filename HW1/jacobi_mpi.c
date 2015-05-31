#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Get residue:
double gres(double *u, int n, double h2){
    double r=0,temp=0,res;
    int i;
    //r= ((2*u[0]-u[1])/h2-1)*((2*u[0]-u[1])/h2-1)+((2*u[n-1]-u[n-2])/h2-1)*((2*u[n-1]-u[n-2])/h2-1);
    for(i=1;i<n+1;i++){
        temp=(2*u[i]-u[i-1]-u[i+1])/h2-1;
        r+=temp*temp;
	}
    MPI_Allreduce(&r,&res,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return sqrt(res);
         
}

int main (int argc, char *argv[])
{  int N=atol(argv[1]); //Mesh size;
   int M=atol(argv[2]);//max number of loops;
 
   int rank,i,k=0,tag=111,p,n;
   double T1,T2, h=1.0/N, h2=h*h, r0=1e-5,r, res;//res< r0 is the stopping criterion;
    MPI_Request req_sl,req_sr,req_rl,req_rr;
    MPI_Status stat_s,stat_r;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (N%p!=0) {
        printf("Mesh size N must be divisible by the number of processors p!\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    n=N/p;
    //printf("n=%d\n",n);
    
    T1=MPI_Wtime();
    double * u=(double *)calloc(sizeof(double), n+2);
    double * u1=(double *)calloc(sizeof(double), n+2);
    //printf("vectors created\n");
    //First iteration, with zero initial guess;
    for (i=0;i<n+2;i++){
        u[i]=h2/2.0;
        //printf("%f ",u[i]);
    }
    //boundary conditions:
    if (rank==0)
	u[0]=0;
    else if (rank==p)
	u[n+1]=0;    

    while (k<M){
        /*
	if (rank==0)
            u1[1]=(u[2]+h2)/2;
        else
            u1[1]=(u[0]+u[2]+h2)/2;
        
        if (rank==p-1)
            u1[n]=(u[n-1]+h2)/2;
        else
	*/
	 u1[1]=(u[0]+u[2]+h2)/2;
	 u1[n]=(u[n-1]+u[n+1]+h2)/2;
        
        //passing to the right process:
        if (rank<p-1){
            MPI_Isend(&u1[n],1,MPI_DOUBLE,rank+1,tag,MPI_COMM_WORLD,&req_sr);
            MPI_Irecv(&u1[n+1],1,MPI_DOUBLE,rank+1,tag,MPI_COMM_WORLD,&req_rr);
        }
        
        //passing to the left process:
        if (rank>0){
            MPI_Isend(&u1[1],1,MPI_DOUBLE,rank-1,tag,MPI_COMM_WORLD,&req_sl);
            MPI_Irecv(&u1[0],1,MPI_DOUBLE,rank-1,tag,MPI_COMM_WORLD,&req_rl);
        }
        
        for(i=2;i<n;i++)
            u1[i]=(u[i-1]+u[i+1]+h2)/2;
        
        if (rank<p-1){
            MPI_Wait(&req_sr,&stat_s);
            MPI_Wait(&req_rr,&stat_r);
        }
        if (rank>0){
            MPI_Wait(&req_sl,&stat_s);
            MPI_Wait(&req_rl,&stat_r);
        }
       
        for (i=0;i<n+2;i++)
            u[i]=u1[i];
        //printf("%f\n",u[0]);
        
        res=gres(u,n,h2);
        //MPI_Allreduce(&r,&res,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        
        if (rank==0 && k%50==0) 
            printf("Loop %d, residue=%f\n",k,res);

	if ( res<r0){
	    if (rank==0) 
                printf("Stopping at loop %d, final residue=%f\n",k,res);
	    break;
	}
        k++;
        //printf("k=%d\n",k);
    }

/* Clean up */
  free(u);
  free(u1);

/* timing */
  MPI_Barrier(MPI_COMM_WORLD);
  T2=MPI_Wtime();
  if (rank==0){
     printf("Residue after %d loops is %f\n",k,res);
     printf("Time elapsed is %f sec.\n",T2-T1);
  }
  MPI_Finalize();
//*/
    return 0;
    
}

