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

    if (child1 < nproc)
    {
        child1 = (child1 + descr.root) % nproc; // mapping to real
        offset = +(child1 - master_root) * max_length;
        if (rcv_buf[offset] != 0)
        {
            return MPI_SUCCESS;
        }
        result = move_data(origin_addr, rcv_buf, offset, child1, descr.message_length, win);
    }
    if (child2 < nproc)
    {

        child2 = (child2 + descr.root) % nproc; // mapping to real
        offset = +(child2 - master_root) * max_length;
        if (rcv_buf[offset] != 0)
        {
            return MPI_SUCCESS;
        }
        result = move_data(origin_addr, rcv_buf, offset, child2, descr.message_length, win);
    }
    return result;
}