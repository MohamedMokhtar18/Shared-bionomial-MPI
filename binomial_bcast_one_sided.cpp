#include "binomial_bcast_one_sided.h"
#define max_length 8388608 /* ==> 2 x 32 MB per process */

void RMA_Bcast_binomial_OneSide(buf_dtype *origin_addr, int my_rank,
                                int i,
                                const descr_t &descr, int nproc,
                                int j,
                                int mid,
                                int length, std::fstream &file,
                                MPI_Win win, MPI_Comm comm)
{
    //? declare arguments
    buf_dtype *rcv_buf;
    MPI_Win_allocate((MPI_Aint)(max_length * sizeof(buf_dtype)), sizeof(buf_dtype), MPI_INFO_NULL, MPI_COMM_WORLD, &rcv_buf, &win);

    int srank = comp_srank(my_rank, descr.root, nproc); // Compute rank relative to root
    auto mask = 1;
    while (mask < nproc)
    {
        if ((srank & mask) == 0)
        { // send data to the next process if bit is not set
            auto rank = srank | mask;
            if (rank < nproc)
            {
                rank = comp_rank(rank, descr.root, nproc); // Compute rank from srank
                                                           // MPI_Win_fence(MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE, win);
                MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, win);

                // ! assign values to rcv_buf pointer
                MPI_Get(rcv_buf, length, MPI_FLOAT, rank, (MPI_Aint)0, length, MPI_FLOAT, win);

                MPI_Win_unlock(rank, win);
                MPI_Win_sync(win);
                // ! add results to the file
                file
                    << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
                    << " snd_buf[0," << mid << "," << (length - 1) << "]"
                    << "=(" << origin_addr[0] << origin_addr[mid] << origin_addr[length - 1] << ")"
                    << "rank " << rank
                    << std::endl;
                file << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
                     << " rcv_buf[0," << mid << "," << (length - 1) << "]"
                     << "=(" << rcv_buf[0] << rcv_buf[mid] << rcv_buf[length - 1] << ")"
                     << "rank " << rank
                     << std::endl;
            }
            else
            {

                break;
            }

            mask = mask << 1;
        }
    }
    MPI_Win_free(&win);
}
// comp_srank: Compute rank relative to root
int comp_srank(int myrank, int root, int nproc)
{
    return (myrank - root + nproc) % nproc;
}

// comp_rank: Compute rank from srank
int comp_rank(int srank, int root, int nproc)
{
    return (srank + root) % nproc;
}