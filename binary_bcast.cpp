#include "binary_bcast.h"
#include "mcs_lock.h"

#define max_length 8388608 /* ==> 2 x 32 MB per process */
MCS_Mutex hdl_binary;      /* Mutex handle */
int BinaryTreeBcast(buf_dtype *origin_addr, buf_dtype *rcv_buf, int my_rank,
                    descr_t descr, int nproc,
                    MPI_Win win, MPI_Comm comm)
{
    int result;
    for (int i = 0; i < nproc; i++)
    {
        descr.root = i;
        my_rank = i;
        result = send_data_binary(origin_addr, rcv_buf, my_rank, descr, nproc, win, comm);
    }

    // MCS_Mutex_unlock(hdl_binary, my_rank);
    // MCS_Mutex_free(&hdl_binary);
    return result;
}
int send_data_binary(buf_dtype *origin_addr, buf_dtype *rcv_buf, int my_rank,
                     descr_t descr, int nproc,
                     MPI_Win win, MPI_Comm comm)
{
    int result;
    int rank = (my_rank - descr.root + nproc) % nproc; // reorder, so fake root is zero
    int offset;
    int master_root = 0;
    int child1 = 2 * rank + 1;
    int child2 = 2 * rank + 2;
    // MCS_Mutex_create(my_rank, comm, &hdl_binary);
    // MCS_Mutex_lock(hdl_binary, my_rank);

    if (child1 < nproc)
    {
        child1 = (child1 + descr.root) % nproc; // mapping to real
        offset = +(child1 - master_root) * max_length;
        if (rcv_buf[offset] != 0)
        {
            return MPI_SUCCESS;
        }
        result = MPI_Win_lock(MPI_LOCK_SHARED, child1, 0, win);
        for (int i = 0; i < descr.message_length; i++)
        {
            rcv_buf[i + offset] = origin_addr[i];
        }
        // *(rcv_buf + (child1 - my_rank)) = *(origin_addr);
        result = MPI_Win_sync(win);
        if (result != MPI_SUCCESS)
        {
            MPI_Abort(comm, result);
        }
        result = MPI_Win_unlock(child1, win);
        if (result != MPI_SUCCESS)
        {
            MPI_Abort(comm, result);
        }
    }
    if (child2 < nproc)
    {

        child2 = (child2 + descr.root) % nproc; // mapping to real
        offset = +(child2 - master_root) * max_length;
        if (rcv_buf[offset] != 0)
        {
            return MPI_SUCCESS;
        }

        result = MPI_Win_lock(MPI_LOCK_SHARED, child2, 0, win);
        for (int i = 0; i < descr.message_length; i++)
        {
            rcv_buf[i + offset] = origin_addr[i];
        }
        // *(rcv_buf + (child2 - my_rank)) = *(origin_addr);
        result = MPI_Win_sync(win);
        if (result != MPI_SUCCESS)
        {
            MPI_Abort(comm, result);
        }
        result = MPI_Win_unlock(child2, win);
        if (result != MPI_SUCCESS)
        {
            MPI_Abort(comm, result);
        }
    }
    return MPI_SUCCESS;
}