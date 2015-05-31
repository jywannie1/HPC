#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>


int main (int argc, char **argv)
{  int N=atol(argv[1]); //number of loops;
   int size,rank,i, sum=0, tag=111,dest,src;
   double T1,T2;
   MPI_Status status;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  T1=MPI_Wtime();
  dest=(rank+1)%size;
  src=(rank+size-1)%size;

  if (rank==0){
    for (i=0;i<N;i++){
      MPI_Send(&sum,1,MPI_INT,dest,tag,MPI_COMM_WORLD);
      MPI_Recv(&sum,1,MPI_INT,src,tag,MPI_COMM_WORLD,&status);
    }
  }
  else {
    for (i=0;i<N;i++){
      MPI_Recv(&sum,1,MPI_INT,src,tag,MPI_COMM_WORLD,&status);
      sum+=rank;
      MPI_Send(&sum,1,MPI_INT,dest,tag,MPI_COMM_WORLD);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  T2=MPI_Wtime();
  if (rank==0){
      if (sum==size*(size-1)/2*N)
        printf("Good!\n");
      else
        printf("Bad!\n");
      //printf("Sum=%d, theoretical value=%d\n", sum, size*(size-1)/2*N);
      printf("Time for N loops= %e sec.\n", T2-T1);
      printf("Time per communication=%e sec.\n", (T2-T1)/(2*N*size));
  }

MPI_Finalize();
return 0;
}
