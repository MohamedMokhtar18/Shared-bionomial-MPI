#pragma once
using buf_dtype = float; /*type of transfered message*/

struct descr_t
{
    // Root (for collectives with root)
    int root = -1;

    // // Buf size for current message (in bytes)
    // int bufsize;
    // *message length for the current message
    int message_length;

    // RMA window's id
    // win_id_t wid = MPI_WIN_NO_ID;
};
inline int move_data(buf_dtype *source, buf_dtype *distination, int offset, int rank, int message_size, MPI_Win win)
{
    int result = MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, win);
    // ! assign values to rcv_buf pointer
    result = MPI_Win_sync(win);
    // int offset = (put_rank - master_root) * max_length;
    for (int k = 0; k < message_size; k++)
        distination[k + offset] = source[k];
    result = MPI_Win_sync(win);

    result = MPI_Win_unlock(rank, win);
    return result;
}