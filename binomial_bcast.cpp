#include "binomial_bcast.h"
#define max_length 8388608 /* ==> 2 x 32 MB per process */

void RMA_Bcast_binomial(buf_dtype *origin_addr, int my_rank,
                        int i,
                        const descr_t &descr, int nproc,
                        int j,
                        int mid,
                        int length, std::fstream &file,
                        MPI_Win win, MPI_Comm comm)
{
    //? declare arguments
    buf_dtype *rcv_buf; // rcv_buf pointer type

    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm);

    MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm, &rcv_buf, &win);

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
                MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, win);
                // ! assign values to rcv_buf pointer

                *(rcv_buf + (rank - my_rank)) = *(origin_addr);
                MPI_Win_sync(win);
                MPI_Win_unlock(rank, win);
                // ! add results to the file
                file << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
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
            else
            {
                // If bit is set, break
                // (in original non-RMA algorithm it's the receive phase)
                break;
            }

            mask = mask << 1;
        }
    }
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