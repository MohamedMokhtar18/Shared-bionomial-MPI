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

    // MCS_Mutex_unlock(hdl_binomial, my_rank);
    // MCS_Mutex_free(&hdl_binomial);
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
    // MCS_Mutex_create(my_rank, comm, &hdl_binomial);
    // MCS_Mutex_lock(hdl_binomial, my_rank);
    while (mask < nproc)
    {
        if ((srank & mask) == 0)
        { // send data to the next process if bit is not set
            auto rank = srank | mask;
            if (rank < nproc)
            {
                rank = comp_rank(rank, descr.root, nproc); // Compute rank from srank
                /*//! offset_left  is defined so that rcv_buf(xxx+offset) in process 'my_rank' is the same location as
/*                              //!  rcv_buf(xxx) in process 'rank':                                       */
                offset = +(rank - master_root) * max_length;
                result = MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, win);
                // ! assign values to rcv_buf pointer
                result = MPI_Win_sync(win);
                for (int k = 0; k < descr.message_length; k++)
                {
                    rcv_buf[k + offset] = origin_addr[k];
                }

                // *(rcv_buf + (rank - my_rank)) = *(origin_addr);
                result = MPI_Win_sync(win);
                if (result != MPI_SUCCESS)
                {
                    MPI_Abort(comm, MPI_ERR_OTHER);
                }
                result = MPI_Win_unlock(rank, win);
                if (result != MPI_SUCCESS)
                {
                    MPI_Abort(comm, MPI_ERR_OTHER);
                }
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
    // MCS_Mutex_unlock(hdl_binomial, my_rank);
    // MCS_Mutex_free(&hdl_binomial);
    return MPI_SUCCESS;
}