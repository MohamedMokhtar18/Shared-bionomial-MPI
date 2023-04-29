#pragma once

#include <mpi.h>
using buf_dtype = float; /*type of transfered message*/
#include "utils.h"
#include <bits/stdc++.h>

int RMA_Bcast_binomial(buf_dtype *origin_addr, buf_dtype *rcv_buf, int my_rank,
                       descr_t &descr, int nproc,
                       MPI_Win win, MPI_Comm comm);

int comp_srank(int myrank, int root, int nproc);
int comp_rank(int srank, int root, int nproc);
int send_loop(buf_dtype *origin_addr, buf_dtype *rcv_buf, int my_rank,
              const descr_t &descr, int nproc,
              MPI_Win win, MPI_Comm comm);