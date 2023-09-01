#pragma once

#include <mpi.h>
#include "utils.h"

using buf_dtype = float; /*type of transfered message*/
#include <bits/stdc++.h>
int RMA_Bcast_binary_OneSide(buf_dtype *origin_addr, MPI_Datatype origin_datatype, MPI_Aint target_disp, int my_rank,
                     descr_t descr, int nproc, MPI_Win win, MPI_Comm comm);;
int send_data_binary_one_side(buf_dtype *origin_addr, MPI_Datatype origin_datatype, int my_rank,
              const descr_t &descr, int nproc,
              MPI_Win win, MPI_Comm comm);