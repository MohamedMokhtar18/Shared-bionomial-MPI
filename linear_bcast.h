#pragma once

#include <mpi.h>
#include "utils.h"
using buf_dtype = float; /*type of transfered message*/
#include <bits/stdc++.h>
#include <thread>
#include <boost/thread/scoped_thread.hpp>

int RMA_Bcast_Linear(buf_dtype *origin_addr, MPI_Datatype origin_datatype, MPI_Aint target_disp,
                     descr_t descr, int nproc, MPI_Win win, MPI_Comm comm);
// static int Compare_data(buf_dtype *buf, buf_dtype *origin_addr, int put_rank, int myrank, MPI_Win win);
bool bcast_isdone(int comp_data, int myrank, MPI_Win &done_data_win);