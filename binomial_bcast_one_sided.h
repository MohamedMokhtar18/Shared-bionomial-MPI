#pragma once

#include <mpi.h>
using buf_dtype = float; /*type of transfered message*/
#include "utils.h"
#include <bits/stdc++.h>
void RMA_Bcast_binomial_OneSide(buf_dtype *origin_addr, int my_rank,
                                int i,
                                const descr_t &descr, int nproc,
                                int j,
                                int mid,
                                int length, std::fstream &file,
                                MPI_Win win, MPI_Comm comm);
int comp_srank(int myrank, int root, int nproc);
int comp_rank(int srank, int root, int nproc);