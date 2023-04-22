#include <mpi.h>
using buf_dtype = float; /*type of transfered message*/
#include <bits/stdc++.h>
int RMA_Bcast_Linear(buf_dtype *origin_addr, buf_dtype *rcv_buf, int my_rank,
                     int i,
                     int nproc,
                     int j,
                     int mid,
                     int length, std::fstream &file,
                     MPI_Win win, MPI_Comm comm);
static int Compare_data(buf_dtype *buf, buf_dtype *origin_addr, int put_rank, int myrank, std::fstream &file, MPI_Win win);