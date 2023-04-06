#include "binary_bcast.h"
#define max_length 8388608 /* ==> 2 x 32 MB per process */

void BinaryTreeBcast(buf_dtype *origin_addr, int my_rank,
                     int i,
                     int root, int nproc,
                     int j,
                     int mid,
                     int length, std::fstream &file,
                     MPI_Win win, MPI_Comm comm)
{
    buf_dtype *rcv_buf; // rcv_buf pointer type

    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm);

    MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm, &rcv_buf, &win);
    int rank = (my_rank - root + nproc) % nproc; // reorder, so fake root is zero

    int child1 = 2 * rank + 1;
    int child2 = 2 * rank + 2;

    if (child1 < nproc)
    {
        child1 = (child1 + root) % nproc; // mapping to real
        MPI_Win_lock(MPI_LOCK_SHARED, child1, 0, win);

        *(rcv_buf + (child1 - my_rank)) = *(origin_addr);
        MPI_Win_unlock(child1, win);
        ;

        file
            << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
            << " snd_buf[0," << mid << "," << (length - 1) << "]"
            << "=(" << origin_addr[0] << origin_addr[mid] << origin_addr[length - 1] << ")"
            << "rank " << child1
            << std::endl;
        file << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
             << " rcv_buf[0," << mid << "," << (length - 1) << "]"
             << "=(" << rcv_buf[0] << rcv_buf[mid] << rcv_buf[length - 1] << ")"
             << "rank " << child1
             << std::endl;
    }
    if (child2 < nproc)
    {

        child2 = (child2 + root) % nproc; // mapping to real
        MPI_Win_lock(MPI_LOCK_SHARED, child2, 0, win);

        *(rcv_buf + (child2 - my_rank)) = *(origin_addr);
        MPI_Win_sync(win);
        MPI_Win_unlock(child2, win);

        file
            << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
            << " snd_buf[0," << mid << "," << (length - 1) << "]"
            << "=(" << origin_addr[0] << origin_addr[mid] << origin_addr[length - 1] << ")"
            << "rank " << child2
            << std::endl;
        file << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
             << " rcv_buf[0," << mid << "," << (length - 1) << "]"
             << "=(" << rcv_buf[0] << rcv_buf[mid] << rcv_buf[length - 1] << ")"
             << "rank " << child2
             << std::endl;
    }
}
