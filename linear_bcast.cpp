#include "linear_bcast.h"
#define max_length 8388608 /* ==> 2 x 32 MB per process */

void RMA_Bcast_Linear(buf_dtype *origin_addr, int my_rank,
                      int i,
                      int nproc,
                      int j,
                      int mid,
                      int length, std::fstream &file,
                      MPI_Win win, MPI_Comm comm)
{

    // float rcv_buf[max_length];
    //? declare arguments
    buf_dtype *rcv_buf;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm);

    MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm, &rcv_buf, &win);
    for (auto rank = 0; rank < nproc; rank++)
    {
        MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, win);
        *(rcv_buf + (rank - my_rank)) = *(origin_addr);
        // MPI_Win_flush(rank, win);
        MPI_Win_sync(win);
        MPI_Win_unlock(rank, win);
        // MPI_Win_sync(win);
        file
            << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
            << " snd_buf[0," << mid << "," << (length - 1) << "]"
            << "=(" << origin_addr[0] << origin_addr[mid] << origin_addr[length - 1] << ")"
            << "rank " << rank
            << std::endl;
        file << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
             << " rcv_buf[0," << mid << "," << (length - 1) << "]"
             << "=(" << (rcv_buf + (rank - my_rank))[0] << (rcv_buf + (rank - my_rank))[mid] << (rcv_buf + (rank - my_rank))[length - 1] << ")"
             << "rank " << rank
             << std::endl;
    }
}
