#include "binomial_bcast.h"
#define max_length 8388608 /* ==> 2 x 32 MB per process */
#include "mcs_lock.h"

MCS_Mutex hdl_binomial; /* Mutex handle */

int RMA_Bcast_binomial(buf_dtype *origin_addr, buf_dtype *rcv_buf, int my_rank,
                       descr_t &descr, int nproc,
                       MPI_Win win, MPI_Comm comm)
{
    int result;
    for (int i = 0; i < nproc; i++)
    {
        descr.root = i;
        my_rank = i;
        result = send_loop(origin_addr, rcv_buf, my_rank, descr, nproc, win, comm);
    }
    return result;
}
// comp_srank: Compute rank relative to root
int comp_srank(int myrank, int root, int nproc)
{
    return (myrank - root + nproc) % nproc;
}

// comp_rank: Compute rank from srank
int comp_rank(int srank, int root, int nproc)
{
    return (srank + root) % nproc;
}

int send_loop(buf_dtype *origin_addr, buf_dtype *rcv_buf, int my_rank,
              const descr_t &descr, int nproc,
              MPI_Win win, MPI_Comm comm)
{
    //? declare arguments
    int master_root = 0;
    int result;
    int srank = comp_srank(my_rank, descr.root, nproc); // Compute rank relative to root
    auto mask = 1;
    int offset;
    while (mask < nproc)
    {
        if ((srank & mask) == 0)
        { // send data to the next process if bit is not set
            auto rank = srank | mask;
            if (rank < nproc)
            {
                rank = comp_rank(rank, descr.root, nproc); // Compute rank from srank
                /*//! offset  is defined so that rcv_buf(xxx+offset) in process 'my_rank' is the same location as
/*                              //!  rcv_buf(xxx) in process 'rank':                                       */
                offset = +(rank - master_root) * max_length;
                if (rcv_buf[offset] != 0)
                    break;
                result = move_data(origin_addr, rcv_buf, offset, rank, descr.message_length, win);
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
