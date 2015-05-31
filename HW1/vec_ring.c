#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define MSGSIZE 500000

int main (int argc, char **argv)
{  int N=atol(argv[1]); //number of loops;
   int size,rank,i,tag=111,dest,src;
   int data[MSGSIZE];
   double T1,T2;
   MPI_Status status;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  T1=MPI_Wtime();
  dest=(rank+1)%size;
  src=(rank+size-1)%size;
    //printf("rank=%d, dest=%d, src=%d\n", rank, dest, src);

  if (rank==0){
   // Initialize data:
    for(i=0; i<MSGSIZE; i++)
        data[i] = i;
    //printf("data initialized.\n");
    for (i=0;i<N;i++){
        MPI_Send(data,MSGSIZE,MPI_INT,dest,tag,MPI_COMM_WORLD);
        //printf("process %d sent data to process %d\n", rank,dest);
        MPI_Recv(data,MSGSIZE,MPI_INT,src,tag,MPI_COMM_WORLD,&status);
        //printf("process %d received data from process %d\n", src,rank);
      }
  }
  else {
    for (i=0;i<N;i++){
       MPI_Recv(data,MSGSIZE,MPI_INT,src,tag,MPI_COMM_WORLD,&status);
       //printf("process %d received data from process %d\n", src,rank);
       MPI_Send(data,MSGSIZE,MPI_INT,dest,tag,MPI_COMM_WORLD);
       //printf("process %d sent data to process %d\n", rank,dest);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  T2=MPI_Wtime();
  
  if (rank==0){
    printf("Time for N loops= %f sec.\n", T2-T1);
      //printf("Sum=%d, theoretical value=%d\n", sum, size*(size-1)/2*N);
    printf("Data communicated per sec=%f GB/s.\n", N*size*MSGSIZE*sizeof(int)/(T2-T1)/1e9);
  }
MPI_Finalize();
return 0;
}
