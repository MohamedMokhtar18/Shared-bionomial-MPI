#include "linear_bcast.h"
#define max_length 8388608 /* ==> 2 x 32 MB per process */

int RMA_Bcast_Linear(buf_dtype *origin_addr, buf_dtype *rcv_buf, int my_rank,
                     int i,
                     int nproc,
                     int j,
                     int mid,
                     int length, std::fstream &file,
                     MPI_Win win, MPI_Comm comm)
{

    // float rcv_buf[max_length];
    //? declare arguments
    for (auto rank = 0; rank < nproc; rank++)
    {
        int result = MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, win);
        if (result != MPI_SUCCESS)
        {
            MPI_Abort(comm, result);
        }
        *(rcv_buf + (rank - my_rank)) = *(origin_addr);
        // MPI_Win_flush(rank, win);
        result = Compare_data(rcv_buf, origin_addr, rank, my_rank, file, win);
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

        // MPI_Win_sync(win);
        // file
        //     << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
        //     << " snd_buf[0," << mid << "," << (length - 1) << "]"
        //     << "=(" << origin_addr[0] << origin_addr[mid] << origin_addr[length - 1] << ")"
        //     << "rank " << rank
        //     << std::endl;

        // file << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
        //      << " rcv_buf[0," << mid << "," << (length - 1) << "]"
        //      << "=(" << (rcv_buf + (rank - my_rank))[0] << (rcv_buf + (rank - my_rank))[mid] << (rcv_buf + (rank - my_rank))[length - 1] << ")"
        //      << "rank " << rank
        //      << std::endl;
    }
    return MPI_SUCCESS;
}

static int Compare_data(buf_dtype *buf, buf_dtype *origin_addr, int put_rank, int myrank, std::fstream &file, MPI_Win win)
{
    MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, MPI_COMM_WORLD, &buf, &win);

    int result = MPI_Compare_and_swap(&buf, origin_addr, origin_addr, MPI_FLOAT, put_rank, (MPI_Aint)max_length * sizeof(buf_dtype), win);
    file << "output " << buf[0] << " orginal " << origin_addr[0] << std::endl;
    //  printf("result orginal %f result %f", origin_addr[0], buf[0]);
    return result;
}
