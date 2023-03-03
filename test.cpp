
#include <stdio.h>
#include <mpi.h>
#include "test.h"
#include <bits/stdc++.h>

#define number_of_messages 50 /*number of messages transfered per package */
#define start_length 4        /*length of the array */
#define length_factor 8
#define max_length 8388608           /* ==> 2 x 32 MB per process */
#define number_package_sizes 8       /*number of package sizes for each size 16,128,...*/
/* #define max_length 67108864    */ /* ==> 2 x 0.5 GB per process */
/* #define number_package_sizes 9 */
using buf_dtype = float; /*type of transfered message*/

int main(int argc, char *argv[])
{ // ? variable declaration
    int my_rank,
        size;
    int srank, rank;
    int i, mid, length, test_value;
    double start, finish, transfer_time;
    int provided = MPI_THREAD_MULTIPLE;
    float snd_buf[max_length];
    MPI_Win win;
    MPI_Comm comm_sm;
    std::fstream file; /* value for result file*/

    // ? Methid declaration
    int comp_srank(int myrank, int root, int nproc);
    int comp_rank(int srank, int root, int nproc);
    void RMA_Bcast_binomial(buf_dtype * origin_addr, int my_rank,
                            int i,
                            MPI_Aint target_disp, int nproc,
                            int j,
                            int mid,
                            int length,
                            std::fstream &file,
                            MPI_Win win, MPI_Comm comm);

    // auto type_size = 0;
    // MPI_Type_size(origin_datatype, &type_size);
    file.open("results.dat", std::ios::out); /*create file and open it*/
    setbuf(stdout, NULL);
    // Todo  create passing value for tybe of broadcast
    // if (argv[1] == NULL)
    // {
    //     throw std::runtime_error(" invalid argument error must be number to be send ");
    // }
    // ? MPI Intialization
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (my_rank == 0)
        printf("    message size      transfertime  duplex bandwidth per process and neighbor\n");
    length = start_length;
    for (int j = 1; j <= number_package_sizes; j++)
    {
        for (i = 0; i <= number_of_messages; i++)
        {
            if (i == 1)
                start = MPI_Wtime(); // start the timer
            test_value = j * 1000000 + i * 10000 + my_rank * 10;
            mid = (length - 1) / number_of_messages * i;
            snd_buf[0] = test_value + 1;
            snd_buf[mid] = test_value + 2;
            snd_buf[length - 1] = test_value + 3;
            // Todo make a generic method for each type to be compared
            RMA_Bcast_binomial((buf_dtype *)snd_buf, my_rank, i, 0, size, j, mid, length, file, win, comm_sm);
        }
        finish = MPI_Wtime();
        if (my_rank == 0)
        {
            transfer_time = (finish - start) / number_of_messages; // calculate transfer message for each pacakge
            file << std::setw(10) << length * sizeof(float) << " bytes " << std::setw(12) << transfer_time * 1e6 << " usec " << std::setw(13) << 1.0e-6 * 2 * length * sizeof(float) / transfer_time << " MB/s" << std::endl;
            printf("%10i bytes %12.3f usec %13.3f MB/s\n",
                   length * (int)sizeof(float), transfer_time * 1e6, 1.0e-6 * 2 * length * sizeof(float) / transfer_time);
        }
        length = length * length_factor;
    }

    //  MPI_Win_free(&win);
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
void RMA_Bcast_binomial(buf_dtype *origin_addr, int my_rank,
                        int i,
                        MPI_Aint target_disp, int nproc,
                        int j,
                        int mid,
                        int length, std::fstream &file,
                        MPI_Win win, MPI_Comm comm)
{
    // float rcv_buf[max_length];
    //? declare arguments
    buf_dtype *rcv_buf; // rcv_buf pointer type
    // int rank;
    descr_t descr;
    descr.root = my_rank;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm);

    MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm, &rcv_buf, &win);

    // descr.bufsize = origin_count * type_size;
    // descr.wid = wid;
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
                // ! assign values to rcv_buf pointer for each
                rcv_buf[0] = (buf_dtype)(origin_addr[0]);
                rcv_buf[mid] = (buf_dtype)(origin_addr[mid]);
                rcv_buf[length - 1] = (buf_dtype)(origin_addr[length - 1]);

                MPI_Win_unlock(rank, win);
                MPI_Win_sync(win);
                // ! add results to the file
                file << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
                     << " snd_buf[0," << mid << "," << (length - 1) << "]"
                     << "=(" << origin_addr[0] << origin_addr[mid] << origin_addr[length - 1] << ")"
                     << std::endl;
                file << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
                     << " rcv_buf[0," << mid << "," << (length - 1) << "]"
                     << "=(" << rcv_buf[0] << rcv_buf[mid] << rcv_buf[length - 1] << ")"
                     << std::endl;

                // printf("%d: j=%d, i=%d --> snd_buf[0,%d,%d]=(%f,%f,%f)\n",
                //        my_rank, j, i, mid, length - 1, origin_addr[0], origin_addr[mid], origin_addr[length - 1]);

                // printf("%d: j=%d, i=%d --> rcv_buf[0,%d,%d]=(%f,%f,%f)\n",
                //        my_rank, j, i, mid, length - 1, rcv_buf[0], rcv_buf[mid], rcv_buf[length - 1]);
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