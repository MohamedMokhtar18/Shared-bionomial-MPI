#include "binary_bcast_one_side.h"
// #include "mcs_lock.h"

#define max_length 8388608 /* ==> 2 x 32 MB per process */
// MCS_Mutex hdl_binary;      /* Mutex handle */
int RMA_Bcast_binary_OneSide(buf_dtype *origin_addr, MPI_Datatype origin_datatype, MPI_Aint target_disp, int my_rank,
                     descr_t descr, int nproc, MPI_Win win, MPI_Comm comm)
{
    int result;
    for (int i = 0; i < nproc; i++)
    {
        descr.root = i;
        my_rank = i;
        result = send_data_binary_one_side(origin_addr, origin_datatype, my_rank, descr, nproc, win, comm);
    }

    return result;
}
int send_data_binary_one_side(buf_dtype *origin_addr, MPI_Datatype origin_datatype, int my_rank,
              const descr_t &descr, int nproc,
              MPI_Win win, MPI_Comm comm)
{
    int result;
    int rank = (my_rank - descr.root + nproc) % nproc; // reorder, so fake root is zero
    int child1 = 2 * rank + 1;
    int child2 = 2 * rank + 2;

    if (child1 < nproc)
    {
        child1 = (child1 + descr.root) % nproc; // mapping to real
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, 0, win);

                // ! assign values to rcv_buf pointer
                MPI_Put(origin_addr, descr.message_length, origin_datatype, child1, 0, descr.message_length, origin_datatype, win);

                MPI_Win_unlock(rank, win);
    }
    if (child2 < nproc)
    {

        child2 = (child2 + descr.root) % nproc; // mapping to real
                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, 0, win);

                // ! assign values to rcv_buf pointer
                MPI_Put(origin_addr, descr.message_length, origin_datatype, child2, 0, descr.message_length, origin_datatype, win);

                MPI_Win_unlock(rank, win);
    }
    return result;
}