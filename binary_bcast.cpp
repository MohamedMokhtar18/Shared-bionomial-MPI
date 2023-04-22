#include "binary_bcast.h"
#define max_length 8388608 /* ==> 2 x 32 MB per process */

int BinaryTreeBcast(buf_dtype *origin_addr, buf_dtype *rcv_buf, int my_rank,
                    int i,
                    int root, int nproc,
                    int j,
                    int mid,
                    int length, std::fstream &file,
                    MPI_Win win, MPI_Comm comm)
{
    int result;
    int rank = (my_rank - root + nproc) % nproc; // reorder, so fake root is zero

    int child1 = 2 * rank + 1;
    int child2 = 2 * rank + 2;

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
        file
            << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
            << " snd_buf[0," << mid << "," << (length - 1) << "]"
            << "=(" << origin_addr[0] << origin_addr[mid] << origin_addr[length - 1] << ")"
            << "rank " << child1
            << std::endl;
        file << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
             << " rcv_buf[0," << mid << "," << (length - 1) << "]"
             << "=(" << (rcv_buf + (child1 - my_rank))[0] << (rcv_buf + (child1 - my_rank))[mid] << (rcv_buf + (child1 - my_rank))[length - 1] << ")"
             << "rank " << child1
             << std::endl;
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
        file
            << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
            << " snd_buf[0," << mid << "," << (length - 1) << "]"
            << "=(" << origin_addr[0] << origin_addr[mid] << origin_addr[length - 1] << ")"
            << "rank " << child2
            << std::endl;
        file << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
             << " rcv_buf[0," << mid << "," << (length - 1) << "]"
             << "=(" << (rcv_buf + (child2 - my_rank))[0] << (rcv_buf + (child2 - my_rank))[mid] << (rcv_buf + (child2 - my_rank))[length - 1] << ")"
             << "rank " << child2
             << std::endl;
    }
    return MPI_SUCCESS;
}
