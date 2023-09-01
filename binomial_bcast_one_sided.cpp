#include "binomial_bcast_one_sided.h"
#define max_length 8388608 /* ==> 2 x 32 MB per process */

void RMA_Bcast_binomial_OneSide(buf_dtype *origin_addr, MPI_Datatype origin_datatype, MPI_Aint target_disp, int my_rank,
                     descr_t descr, int nproc, MPI_Win win, MPI_Comm comm)
{
     int result;
    for (int i = 0; i < nproc; i++)
    {
        descr.root = i;
        my_rank = i;
        result = send_loop(origin_addr,origin_datatype ,my_rank, descr, nproc, win, comm);
    }
}

int send_loop(buf_dtype *origin_addr, MPI_Datatype origin_datatype,int my_rank,
              const descr_t &descr, int nproc,
              MPI_Win win, MPI_Comm comm)
{
    //? declare arguments
    int result;
    int srank = comp_srank(my_rank, descr.root, nproc); // Compute rank relative to root
    auto mask = 1;
    while (mask < nproc)
    {
        if ((srank & mask) == 0)
        { // send data to the next process if bit is not set
            auto rank = srank | mask;
            if (rank < nproc)
            {
                rank = comp_rank(rank, descr.root, nproc); // Compute rank from srank
                 MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, 0, win);

                // ! assign values to rcv_buf pointer
                // MPI_Put(rcv_buf, descr.message_length, MPI_FLOAT, rank, (MPI_Aint)0, descr.message_length, MPI_FLOAT, win);
                MPI_Put(origin_addr, descr.message_length, origin_datatype, rank, 0, descr.message_length, origin_datatype, win);

                MPI_Win_unlock(rank, win);
            }
            else
            {

                // If bit is set, break
                // (in original non-RMA algorithm it's the receive phase)
                break;
            }

            mask = mask << 1;
        }
    }
    return MPI_SUCCESS;
}
