#include <mpi.h>
using buf_dtype = float; /*type of transfered message*/
#include <bits/stdc++.h>
void RMA_Bcast_Linear(buf_dtype *origin_addr, int my_rank,
                      int i,
                      int nproc,
                      int j,
                      int mid,
                      int length, std::fstream &file,
                      MPI_Win win, MPI_Comm comm);
