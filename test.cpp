
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

    int i, mid, length, test_value;
    double start, finish, transfer_time;
    int provided = MPI_THREAD_MULTIPLE;
    float snd_buf[max_length];
    descr_t descr;

    MPI_Comm tree_comm;
    MPI_Win win;
    MPI_Comm comm_sm;
    std::fstream file; /* value for result file*/
    std::fstream fileBench;

    // ? Methid declaration
    int comp_srank(int myrank, int root, int nproc);
    int comp_rank(int srank, int root, int nproc);
    void RMA_Bcast_binomial(buf_dtype * origin_addr, int my_rank,
                            int i,
                            const descr_t &descr, int nproc,
                            int j,
                            int mid,
                            int length,
                            std::fstream &file,
                            MPI_Win win, MPI_Comm comm);
    void RMA_Bcast_binomial_OneSide(buf_dtype * origin_addr, int my_rank,
                                    int i,
                                    const descr_t &descr, int nproc,
                                    int j,
                                    int mid,
                                    int length, std::fstream &file,
                                    MPI_Win win, MPI_Comm comm);
    void RMA_Bcast_Linear(buf_dtype * origin_addr, int my_rank,
                          int i,
                          const descr_t &descr, int nproc,
                          int j,
                          int mid,
                          int length, std::fstream &file,
                          MPI_Win win, MPI_Comm comm);
    void BinaryTreeBcast(buf_dtype * origin_addr, int my_rank,
                         int i,
                         const descr_t &descr, int nproc,
                         int j,
                         int mid,
                         int length, std::fstream &file,
                         MPI_Win win, MPI_Comm comm);
    void RingBcast(buf_dtype * origin_addr, int my_rank,
                   int i,
                   const descr_t &descr, int nproc,
                   int j,
                   int mid,
                   int length, std::fstream &file,
                   MPI_Win win, MPI_Comm comm);

    // auto type_size = 0;
    // MPI_Type_size(origin_datatype, &type_size);
    file.open("resultRing200.dat", std::ios::out); /*create file and open it*/
    fileBench.open("resultRingBench200.dat", std::ios::out);
    setbuf(stdout, NULL);
    // Todo  create passing value for tybe of broadcast
    // if (argv[1] == NULL)
    // {
    //     throw std::runtime_error(" invalid argument error must be number to be send ");
    // }
    // ? MPI Intialization
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    // MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    descr.root = my_rank;
    // Create tree topology
    // int dims[1] = {size};
    // int periods[1] = {1};
    // MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 0, &tree_comm); // create tree topology

    if (my_rank == 0)
    {
        printf("    message size      transfertime  duplex bandwidth per process and neighbor\n");
        fileBench << "    message size      transfertime  duplex bandwidth per process and neighbor" << std::endl;
    }
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
            // RMA_Bcast_binomial((buf_dtype *)snd_buf, my_rank, i, descr, size, j, mid, length, file, win, comm_sm);
            RingBcast((buf_dtype *)snd_buf, my_rank, i, descr, size, j, mid, length, file, win, comm_sm);
            //  MPI_Win_allocate((MPI_Aint)(max_length * sizeof(buf_dtype)), sizeof(buf_dtype), MPI_INFO_NULL, MPI_COMM_WORLD, &snd_buf, &win);
            // BinaryTreeBcast((buf_dtype *)snd_buf, my_rank, i, descr, size, j, mid, length, file, win, comm_sm);

            // RMA_Bcast_Linear((buf_dtype *)snd_buf, my_rank, i, descr, size, j, mid, length, file, win, comm_sm);
            // MPI_Win_free(&win);
        }
        finish = MPI_Wtime();
        if (my_rank == 0)
        {
            transfer_time = (finish - start) / number_of_messages; // calculate transfer message for each pacakge
            fileBench << std::setw(10) << length * sizeof(float) << " bytes " << std::setw(12) << transfer_time * 1e6 << " usec " << std::setw(13) << 1.0e-6 * 2 * length * sizeof(float) / transfer_time << " MB/s" << std::endl;
            printf("%10i bytes %12.3f usec %13.3f MB/s\n",
                   length * (int)sizeof(float), transfer_time * 1e6, 1.0e-6 * 2 * length * sizeof(float) / transfer_time);
        }
        length = length * length_factor;
    }
    // delete[] snd_buf;

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
                        const descr_t &descr, int nproc,
                        int j,
                        int mid,
                        int length, std::fstream &file,
                        MPI_Win win, MPI_Comm comm)
{
    // float rcv_buf[max_length];
    //? declare arguments
    buf_dtype *rcv_buf; // rcv_buf pointer type
    // int rank;

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
                // ! assign values to rcv_buf pointer

                *(rcv_buf + (rank - my_rank)) = *(origin_addr);
                // rcv_buf[0] = (buf_dtype)(origin_addr[0]);
                // rcv_buf[mid] = (buf_dtype)(origin_addr[mid]);
                // rcv_buf[length - 1] = (buf_dtype)(origin_addr[length - 1]);

                MPI_Win_unlock(rank, win);
                MPI_Win_sync(win);
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
void RMA_Bcast_binomial_OneSide(buf_dtype *origin_addr, int my_rank,
                                int i,
                                const descr_t &descr, int nproc,
                                int j,
                                int mid,
                                int length, std::fstream &file,
                                MPI_Win win, MPI_Comm comm)
{
    // float rcv_buf[max_length];
    //? declare arguments
    buf_dtype *rcv_buf; // rcv_buf pointer type
                        // int rank;
                        // MPI_Win_create_dynamic(MPI_INFO_NULL, MPI_COMM_WORLD, &win);
                        // MPI_Win_attach(win, rcv_buf, (MPI_Aint)max_length * sizeof(buf_dtype));
                        // MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm);
    MPI_Win_allocate((MPI_Aint)(max_length * sizeof(buf_dtype)), sizeof(buf_dtype), MPI_INFO_NULL, MPI_COMM_WORLD, &rcv_buf, &win);
    //  MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm, &rcv_buf, &win);
    //  descr.bufsize = origin_count * type_size;
    //  descr.wid = wid;
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
                // MPI_Win_fence(MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE, win);
                // rcv_buf[0] = (buf_dtype)(origin_addr[0]);
                // rcv_buf[mid] = (buf_dtype)(origin_addr[mid]);
                // rcv_buf[length - 1] = (buf_dtype)(origin_addr[length - 1]);
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
    MPI_Win_free(&win);
}

void RMA_Bcast_Linear(buf_dtype *origin_addr, int my_rank,
                      int i,
                      const descr_t &descr, int nproc,
                      int j,
                      int mid,
                      int length, std::fstream &file,
                      MPI_Win win, MPI_Comm comm)
{

    // float rcv_buf[max_length];
    //? declare arguments
    buf_dtype *rcv_buf; // rcv_buf pointer type
                        // int rank;
                        // MPI_Win_create_dynamic(MPI_INFO_NULL, MPI_COMM_WORLD, &win);
                        // MPI_Win_attach(win, rcv_buf, (MPI_Aint)max_length * sizeof(buf_dtype));
                        // MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm);
                        //  MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm, &rcv_buf, &win);
                        //  descr.bufsize = origin_count * type_size;
                        //  descr.wid = wid;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm);

    MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm, &rcv_buf, &win);
    for (auto rank = 0; rank < nproc; rank++)
    {
        MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, win);
        *(rcv_buf + (rank - my_rank)) = *(origin_addr);

        MPI_Win_unlock(rank, win);
        MPI_Win_sync(win);
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

void BinaryTreeBcast(buf_dtype *origin_addr, int my_rank,
                     int i,
                     const descr_t &descr, int nproc,
                     int j,
                     int mid,
                     int length, std::fstream &file,
                     MPI_Win win, MPI_Comm comm)
{
    buf_dtype *rcv_buf; // rcv_buf pointer type

    // int size, real_rank;

    // MPI_Comm_size(comm, &size); // get num of processes
    // if (root < 0 || root >= size) {
    //     return MPI_ERR_ROOT;
    // }
    // // otherwise - let's do some send-recv
    // MPI_Comm_rank(comm, &real_rank);
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm);

    MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm, &rcv_buf, &win);
    int rank = (my_rank - descr.root + nproc) % nproc; // reorder, so fake root is zero

    // if (rank != 0) { // not root
    //     // recv
    //     int from = ((rank - 1) / 2 + descr.root ) % nproc; // fake rank, map to real
    //     MPI_Status status;
    //     MPI_Recv(buffer, count, datatype, from, TAG + real_rank, comm, &status);

    // }
    int child1 = 2 * rank + 1;
    int child2 = 2 * rank + 2;

    if (child1 < nproc)
    {
        child1 = (child1 + descr.root) % nproc; // mapping to real
        MPI_Win_lock(MPI_LOCK_SHARED, child1, 0, win);

        // MPI_Put(origin_addr, length, MPI_FLOAT, child1, (MPI_Aint)0, length, MPI_FLOAT, win);
        *(rcv_buf + (child1 - my_rank)) = *(origin_addr);
        MPI_Win_unlock(child1, win);
        MPI_Win_sync(win);
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

        child2 = (child2 + descr.root) % nproc; // mapping to real
        MPI_Win_lock(MPI_LOCK_SHARED, child2, 0, win);

        // MPI_Put(origin_addr, length, MPI_FLOAT, child1, (MPI_Aint)0, length, MPI_FLOAT, win);
        *(rcv_buf + (child2 - my_rank)) = *(origin_addr);
        MPI_Win_unlock(child2, win);
        MPI_Win_sync(win);

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
void RingBcast(buf_dtype *origin_addr, int my_rank,
               int i,
               const descr_t &descr, int nproc,
               int j,
               int mid,
               int length, std::fstream &file,
               MPI_Win win, MPI_Comm comm)
{
    buf_dtype *rcv_buf; // rcv_buf pointer type

    int right = (my_rank + 1) % nproc;
    int left = (my_rank - 1 + nproc) % nproc;
    /* ... this SPMD-style neighbor computation with modulo has the same meaning as: */
    /* right = my_rank + 1;          */
    /* if (right == size) right = 0; */
    /* left = my_rank - 1;           */
    /* if (left == -1) left = size-1;*/

    /* Create the window. */

    /* ParaStation MPI may no allow MPI_Win_allocate_shared on MPI_COMM_WORLD. Workaround: Substitute MPI_COMM_WORLD by comm_sm: */
    MPI_Win_allocate_shared((MPI_Aint)(1 * sizeof(int)), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &rcv_buf, &win);

    for (i = 0; i < nproc; i++)
    {
        MPI_Win_lock(MPI_LOCK_SHARED, right, 0, win);
        /* MPI_Put(&snd_buf, 1, MPI_INT, right, (MPI_Aint) 0, 1, MPI_INT, win); */
        /*   ... is substited by (with offset "right-my_rank" to store into right neigbor's rcv_buf): */
        *(rcv_buf + (right - my_rank)) = *origin_addr;

        MPI_Win_unlock(right, win);
        MPI_Win_sync(win);
        file << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
             << " snd_buf[0," << mid << "," << (length - 1) << "]"
             << "=(" << origin_addr[0] << origin_addr[mid] << origin_addr[length - 1] << ")"
             << "rank " << right
             << std::endl;
        file << " " << my_rank << ": j=" << j << ", i=" << i << " --> "
             << " rcv_buf[0," << mid << "," << (length - 1) << "]"
             << "=(" << (rcv_buf + (right - my_rank))[0] << (rcv_buf + (right - my_rank))[mid] << (rcv_buf + (right - my_rank))[length - 1] << ")"
             << "rank " << right
             << std::endl;
    }

    // MPI_Win_free(&win);
}