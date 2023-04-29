#include "linear_bcast.h"
#include "mcs_lock.h"
#define max_length 8388608 /* ==> 2 x 32 MB per process */
MCS_Mutex hdl;             /* Mutex handle */

int RMA_Bcast_Linear(buf_dtype *origin_addr, buf_dtype *rcv_buf, int my_rank, descr_t descr,
                     int nproc,
                     MPI_Win win, MPI_Comm comm)
{
    // MCS_Mutex_create(my_rank, comm, &hdl);
    // MCS_Mutex_lock(hdl, my_rank);
    int master_root = 0;
    int offset;
    // float rcv_buf[max_length];
    //? declare arguments
    for (auto rank = 0; rank < nproc; rank++)
    {
        int result = MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, win);
        if (result != MPI_SUCCESS)
        {
            MPI_Abort(comm, result);
        }
        offset = +(rank - master_root) * max_length;

        for (int i = 0; i < descr.message_length; i++)
        {
            rcv_buf[i + offset] = origin_addr[i];
        }
        // *(rcv_buf + (rank - my_rank)) = *(origin_addr);
        result = MPI_Win_sync(win);
        if (result != MPI_SUCCESS)
        {
            MPI_Abort(comm, result);
        }
        result = MPI_Win_unlock(rank, win);
        if (result != MPI_SUCCESS)
        {
            MPI_Abort(comm, result);
        }
    }
    // MCS_Mutex_unlock(hdl, my_rank);
    // MCS_Mutex_free(&hdl);
    return MPI_SUCCESS;
}

bool bcast_isdone(int comp_data, int myrank, MPI_Win &done_data_win)
{
    int read_data = 0;
    MPI_Fetch_and_op(NULL, &read_data, MPI_INT, myrank, 0, MPI_NO_OP,
                     done_data_win);

    if (read_data == comp_data - 1)
        return true;
    else
        return false;
}
