/******************************************************************************
* FILE: mpi_bug1.c
* DESCRIPTION: 
*   This program has a bug that causes it to hang.
* AUTHOR: Blaise Barney 
* LAST REVISED: 04/13/05
******************************************************************************/

/* Comment:
The program hangs because processor 0 is waiting for receiving from processsor 1, while processor 1 is waiting for receiving from processor 0 before sending to processor 0, thus both processors are waiting. 
Changing these into isend/irecv can solve the problem.
But I didn't figure out how to add the print "count" command...
*/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[])
{
int numtasks, rank, dest, tag, source, rc, count;
char inmsg, outmsg='x';
MPI_Status stat[2];//stat_in[2], stat_out[2];
MPI_Request req[2];//req_in[2], req_out[2];

MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
printf("Task %d starting...\n",rank);

if (rank == 0) 
  if (numtasks > 2) 
    printf("Numtasks=%d. Only 2 needed. Ignoring extra...\n",numtasks);

if (rank<2){
  dest = (rank + 1)%2;
  source = dest;
  tag = rank;
  //rc = MPI_Isend(&outmsg, 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD,&req_out[rank]);
  rc = MPI_Isend(&outmsg, 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD,&req[0]);
  MPI_Wait(&req[0],&stat[0]);
  printf("Task %d sent to task %d...\n",rank,dest);
  //rc = MPI_Irecv(&inmsg, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, &req_in[rank]);
  rc = MPI_Irecv(&inmsg, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, &req[1]);
  printf("Task %d received from task %d...\n",rank,source);
}
  //printf("Task %d: Received %d char(s) from task %d with tag %d \n",rank, count, stat.MPI_SOURCE, stat.MPI_TAG);
  MPI_Finalize();
  return 0;
}
