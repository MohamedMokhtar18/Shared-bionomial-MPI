#include "binary_bcast.h"
#include "mcs_lock.h"
#define max_length 8388608 /* ==> 2 x 32 MB per process */
MCS_Mutex hdl_binary;      /* Mutex handle */
int BinaryTreeBcast(buf_dtype *origin_addr, buf_dtype *rcv_buf, int my_rank,
                    int root, int nproc,
                    MPI_Win win, MPI_Comm comm)
{
    int result;
    int rank = (my_rank - root + nproc) % nproc; // reorder, so fake root is zero

    int child1 = 2 * rank + 1;
    int child2 = 2 * rank + 2;
    MCS_Mutex_create(my_rank, comm, &hdl_binary);
    MCS_Mutex_lock(hdl_binary, my_rank);

    if (child1 < nproc)
    {
        child1 = (child1 + root) % nproc; // mapping to real
        result = MPI_Win_lock(MPI_LOCK_SHARED, child1, 0, win);

        *(rcv_buf + (child1 - my_rank)) = *(origin_addr);
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

        child2 = (child2 + root) % nproc; // mapping to real
        result = MPI_Win_lock(MPI_LOCK_SHARED, child2, 0, win);

        *(rcv_buf + (child2 - my_rank)) = *(origin_addr);
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
    MCS_Mutex_unlock(hdl_binary, my_rank);
    MCS_Mutex_free(&hdl_binary);
    return MPI_SUCCESS;
}
