#include <mpi.h>
using buf_dtype = float; /*type of transfered message*/
#include <bits/stdc++.h>
int BinaryTreeBcast(buf_dtype *origin_addr, buf_dtype *rcv_buf, int my_rank,
                    int root, int nproc,
                    MPI_Win win, MPI_Comm comm);
int send_data_binary(buf_dtype *origin_addr, buf_dtype *rcv_buf, int my_rank,
                     int root, int nproc,
                     MPI_Win win, MPI_Comm comm);