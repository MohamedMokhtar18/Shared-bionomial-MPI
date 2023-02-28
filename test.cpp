
#include <stdio.h>
#include <mpi.h>
#include "test.h"
int main(int argc, char *argv[])
{
    int comp_srank(int myrank, int root, int nproc);
    int comp_rank(int srank, int root, int nproc);
    void RMA_Bcast_binomial(const void *origin_addr, int my_rank,
                            MPI_Datatype origin_datatype,
                            MPI_Aint target_disp, int nproc,
                            MPI_Datatype target_datatype,
                            MPI_Win win, MPI_Comm comm);

    // auto type_size = 0;
    // MPI_Type_size(origin_datatype, &type_size);
    int my_rank,
        size;
    int srank, rank;
    int sum, i;
    int provided = MPI_THREAD_MULTIPLE;
    if (argv[1] == NULL)
    {
        throw std::runtime_error(" invalid argument error must be number to be send ");
    }

    int snd_buf = std::stoi(std::string(argv[1]));
    MPI_Win win;
    MPI_Comm comm_sm;
    setbuf(stdout, NULL);
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    RMA_Bcast_binomial(&snd_buf, my_rank, MPI_INT, 0, size, MPI_INT, win, comm_sm);
    MPI_Finalize();
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
void RMA_Bcast_binomial(const void *origin_addr, int my_rank,
                        MPI_Datatype origin_datatype,
                        MPI_Aint target_disp, int nproc,
                        MPI_Datatype target_datatype,
                        MPI_Win win, MPI_Comm comm)
{
    int *rcv_buf_ptr;
    // int rank;
    descr_t descr;
    descr.root = my_rank;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm);

    MPI_Win_allocate_shared((MPI_Aint)sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &rcv_buf_ptr, &win);

    // descr.bufsize = origin_count * type_size;
    // descr.wid = wid;
    int srank = comp_srank(my_rank, descr.root, nproc);
    auto mask = 1;
    while (mask < nproc)
    {
        if ((srank & mask) == 0)
        { // send data to the next process if bit is not set
            auto rank = srank | mask;
            if (rank < nproc)
            {
                rank = comp_rank(rank, descr.root, nproc);
                MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, win);
                *rcv_buf_ptr = *(int *)(origin_addr);

                MPI_Win_unlock(rank, win);
                MPI_Win_sync(win);

                printf("%i Root %i\t with value %i\n", my_rank, rank, *rcv_buf_ptr);
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