#include "gtest/gtest.h"
#include "binomial_bcast.h"
#include "linear_bcast.h"
#include "binary_bcast.h"
#define start_length 4     /*length of the array */
#define max_length 8388608 /* ==> 2 x 32 MB per process */

std::fstream filetestbinomial, filetestbinary, filetestlinear; /* value for result file*/
TEST(BinaryTreeBcast, BroadcastToOtherProcesses)
{
    descr_t descr;
    int length = start_length;
    descr.message_length = length;
    int mid = 2;
    int my_rank, size;
    MPI_Win win;
    descr.root = 0;
    MPI_Comm comm_sm;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int result = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_sm);

    buf_dtype snd_buf[max_length];
    buf_dtype *rcv_buf;
    snd_buf[0] = 4.2 + 1;
    snd_buf[mid] = 4.2 + 2;
    snd_buf[4 - 1] = 4.2 + 3;
    int res = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm_sm, &rcv_buf, &win);
    if (my_rank == 0)
    {
        res = BinaryTreeBcast(snd_buf, rcv_buf, my_rank, descr, size, win, comm_sm);
        MPI_Win_flush(my_rank, win);
        MPI_Win_sync(win);
        filetestbinary << size << " processes" << std::endl;
        filetestbinary << "rank 0 Bcast([" << snd_buf[0] << "," << snd_buf[mid] << "," << snd_buf[4 - 1] << "])" << std::endl;
        filetestbinary << "rank " << my_rank << "-" << size << ":[" << rcv_buf[0] << "," << rcv_buf[mid] << "," << rcv_buf[4 - 1] << "]" << std::endl;
    }
    else
    {
        // usleep((my_rank + 1) * 50000);
        MPI_Win_flush(my_rank, win);
        MPI_Win_sync(win);
        filetestbinary << "rank " << my_rank << "-" << size << ":[" << rcv_buf[0] << "," << rcv_buf[mid] << "," << rcv_buf[4 - 1] << "]" << std::endl;
    }

    MPI_Win_free(&win);
    // MPI_Finalize();
    EXPECT_EQ(res, MPI_SUCCESS);
}
TEST(RMA_Bcast_binomial, BroadcastToOtherProcesses)
{
    descr_t descr;
    int length = start_length;
    descr.message_length = length;
    int mid = 2;
    int my_rank, size, disp_unit;
    MPI_Win win;
    MPI_Aint buf_size;
    descr.root = 0;
    MPI_Comm comm_sm;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int result = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_sm);

    buf_dtype snd_buf[max_length];
    buf_dtype *rcv_buf;
    snd_buf[0] = 4.2 + 1;
    snd_buf[mid] = 4.2 + 2;
    snd_buf[4 - 1] = 4.2 + 3;
    int res = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm_sm, &rcv_buf, &win);
    MPI_Win_shared_query(win, my_rank, &buf_size, &disp_unit, &rcv_buf);
    if (my_rank == 0)
    {
        res = RMA_Bcast_binomial(snd_buf, rcv_buf, my_rank, descr, size, win, comm_sm);
        filetestbinomial << size << " processes" << std::endl;
        filetestbinomial << "rank 0 Bcast([" << snd_buf[0] << "," << snd_buf[mid] << "," << snd_buf[4 - 1] << "])" << std::endl;
        filetestbinomial << "rank " << my_rank << "-" << size << ":[" << rcv_buf[0] << "," << rcv_buf[mid] << "," << rcv_buf[4 - 1] << "]" << std::endl;
        MPI_Win_flush(my_rank, win);
        MPI_Win_sync(win);
    }
    else
    {
        // usleep((my_rank + 1) * 50000);
        MPI_Win_flush(my_rank, win);
        MPI_Win_sync(win);
        filetestbinomial << "rank " << my_rank << "-" << size << ":[" << rcv_buf[0] << "," << rcv_buf[mid] << "," << rcv_buf[4 - 1] << "]" << std::endl;
    }

    MPI_Win_free(&win);
    //  MPI_Finalize();
    EXPECT_EQ(res, MPI_SUCCESS);
}
TEST(RMA_Bcast_Linear, BroadcastToOtherProcesses)
{
    descr_t descr;
    int length = start_length;
    descr.message_length = length;
    int mid = 2;
    int my_rank, size, disp_unit;
    MPI_Win win;
    descr.root = 0;
    MPI_Comm comm_sm;
    MPI_Aint buf_size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int result = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_sm);

    buf_dtype snd_buf[max_length];
    buf_dtype *rcv_buf;
    snd_buf[0] = 4.2 + 1;
    snd_buf[mid] = 4.2 + 2;
    snd_buf[4 - 1] = 4.2 + 3;
    int res = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm_sm, &rcv_buf, &win);
    MPI_Win_shared_query(win, my_rank, &buf_size, &disp_unit, &rcv_buf);

    if (my_rank == 0)
    {
        res = RMA_Bcast_Linear((buf_dtype *)snd_buf, MPI_FLOAT, buf_size, descr, size, win, comm_sm);
        ;
        MPI_Win_flush(my_rank, win);
        MPI_Win_sync(win);

        filetestlinear << size << " processes" << std::endl;
        filetestlinear << "rank 0 Bcast([" << snd_buf[0] << "," << snd_buf[mid] << "," << snd_buf[4 - 1] << "])" << std::endl;
        filetestlinear << "rank " << my_rank << "-" << size << ":[" << rcv_buf[0] << "," << rcv_buf[mid] << "," << rcv_buf[4 - 1] << "]" << std::endl;
    }
    else
    {
        // usleep((my_rank + 1) * 50000);
        MPI_Win_flush(my_rank, win);
        MPI_Win_sync(win);
        filetestlinear << "rank " << my_rank << "-" << size << ":[" << rcv_buf[0] << "," << rcv_buf[mid] << "," << rcv_buf[4 - 1] << "]" << std::endl;
    }

    MPI_Win_free(&win);
    // MPI_Finalize();
    EXPECT_EQ(res, MPI_SUCCESS);
}