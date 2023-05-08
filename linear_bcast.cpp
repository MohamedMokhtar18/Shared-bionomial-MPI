#include "linear_bcast.h"
#include "mcs_lock.h"
#define max_length 8388608 /* ==> 2 x 32 MB per process */
MCS_Mutex hdl;             /* Mutex handle */

int RMA_Bcast_Linear(buf_dtype *origin_addr, MPI_Datatype origin_datatype, MPI_Aint target_disp,
                     descr_t descr, int nproc, MPI_Win win, MPI_Comm comm)
{
    // MCS_Mutex_create(my_rank, comm, &hdl);
    // MCS_Mutex_lock(hdl, my_rank);
    int master_root = 0;
    //? declare arguments
    for (auto rank = 0; rank < nproc; rank++)
    {
        int result = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, 0, win);
        if (result != MPI_SUCCESS)
            MPI_Abort(comm, result);

        MPI_Put(origin_addr, descr.message_length, origin_datatype, rank, 0, descr.message_length, origin_datatype, win);
        result = MPI_Win_sync(win);
        if (result != MPI_SUCCESS)
            MPI_Abort(comm, result);

        result = MPI_Win_unlock(rank, win);
        if (result != MPI_SUCCESS)
            MPI_Abort(comm, result);
    }

    return MPI_SUCCESS;
}
