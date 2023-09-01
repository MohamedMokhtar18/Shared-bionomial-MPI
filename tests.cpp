// #include "gtest/gtest.h"
// #include "binomial_bcast.h"
// #include "linear_bcast.h"
// #include "binary_bcast.h"
// #include <random>
// #define start_length 4     /*length of the array */
// #define max_length 8388608 /* ==> 2 x 32 MB per process */

// std::fstream filetestbinomial, filetestbinary, filetestlinear; /* value for result file*/
// TEST(BinaryTreeBcast, BroadcastToOtherProcesses)
// {

//     //?intiallize rank and world
//     int my_rank, size;
//     int result = MPI_Comm_size(MPI_COMM_WORLD, &size);
//     result = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//     //?allocation window
//     MPI_Comm comm_sm;
//     MPI_Win win;
//     MPI_Aint buf_size;
//     buf_dtype *rcv_buf;
//     int disp_unit;
//     result = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_sm);
//     result = MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm_sm, &rcv_buf, &win);
//     result = MPI_Win_shared_query(win, 0, &buf_size, &disp_unit, &rcv_buf);

//     //? send buffer intialization with values
//     buf_dtype snd_buf[max_length];
//     int mid = 2;
//     snd_buf[0] = 4.2 + 1;
//     snd_buf[mid] = 4.2 + 2;
//     snd_buf[4 - 1] = 4.2 + 3;
//     descr_t descr;
//     int length = start_length;
//     descr.message_length = length;
//     descr.root = 0;
//     if (my_rank == 0)
//     {
//         result = BinaryTreeBcast(snd_buf, rcv_buf, my_rank, descr, size, win, comm_sm);
//         //! flush memory
//         MPI_Win_flush(my_rank, win);
//         MPI_Win_sync(win);
//         //! put result into file for root
//         filetestbinary << size << " processes" << std::endl;
//         filetestbinary << "rank 0 Bcast([" << snd_buf[0] << "," << snd_buf[mid] << "," << snd_buf[4 - 1] << "])" << std::endl;
//         filetestbinary << "rank " << my_rank << "-" << size << ":[" << rcv_buf[0] << "," << rcv_buf[mid] << "," << rcv_buf[4 - 1] << "]" << std::endl;
//     }
//     else
//     {
//         // usleep((my_rank + 1) * 50000);
//         //! flush memory
//         MPI_Win_flush(my_rank, win);
//         MPI_Win_sync(win);
//         //! put result into file for recived
//         filetestbinary << "rank " << my_rank << "-" << size << ":[" << rcv_buf[0] << "," << rcv_buf[mid] << "," << rcv_buf[4 - 1] << "]" << std::endl;
//     }

//     MPI_Win_free(&win);
//     // MPI_Finalize();
//     EXPECT_EQ(result, MPI_SUCCESS);
// }
// TEST(RMA_Bcast_binomial, BroadcastToOtherProcesses)
// {
//     //?intiallize rank and world
//     int my_rank, size;
//     int result = MPI_Comm_size(MPI_COMM_WORLD, &size);
//     result = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//     //?allocation window
//     MPI_Win win;
//     MPI_Aint buf_size;
//     MPI_Comm comm_sm;
//     int disp_unit;
//     buf_dtype *rcv_buf;
//     result = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_sm);
//     MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm_sm, &rcv_buf, &win);
//     MPI_Win_shared_query(win, 0, &buf_size, &disp_unit, &rcv_buf);
//     //? send buffer intialization with values
//     int length = start_length;
//     descr_t descr;
//     descr.root = 0;
//     descr.message_length = length;
//     int mid = 2;
//     buf_dtype snd_buf[max_length];
//     snd_buf[0] = 4.2 + 1;
//     snd_buf[mid] = 4.2 + 2;
//     snd_buf[4 - 1] = 4.2 + 3;
//     if (my_rank == 0)
//     {
//         result = RMA_Bcast_binomial(snd_buf, rcv_buf, my_rank, descr, size, win, comm_sm);
//         //! flush memory
//         MPI_Win_flush(my_rank, win);
//         MPI_Win_sync(win);
//         //! put result into file for root
//         filetestbinomial << size << " processes" << std::endl;
//         filetestbinomial << "rank 0 Bcast([" << snd_buf[0] << "," << snd_buf[mid] << "," << snd_buf[4 - 1] << "])" << std::endl;
//         filetestbinomial << "rank " << my_rank << "-" << size << ":[" << rcv_buf[0] << "," << rcv_buf[mid] << "," << rcv_buf[4 - 1] << "]" << std::endl;
//     }
//     else
//     {
//         //! flush memory
//         // usleep((my_rank + 1) * 50000);
//         MPI_Win_flush(my_rank, win);
//         MPI_Win_sync(win);
//         //! put result into file for root

//         filetestbinomial << "rank " << my_rank << "-" << size << ":[" << rcv_buf[0] << "," << rcv_buf[mid] << "," << rcv_buf[4 - 1] << "]" << std::endl;
//     }

//     MPI_Win_free(&win);
//     //  MPI_Finalize();
//     EXPECT_EQ(result, MPI_SUCCESS);
// }
// TEST(RMA_Bcast_Linear, BroadcastToOtherProcesses)
// {
//     //?intiallize rank and world
//     int my_rank, size;
//     int result = MPI_Comm_size(MPI_COMM_WORLD, &size);
//     result = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//     //?allocation window
//     MPI_Win win;
//     MPI_Aint buf_size;
//     MPI_Comm comm_sm;
//     int disp_unit;
//     buf_dtype *rcv_buf;
//     result = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_sm);
//     MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm_sm, &rcv_buf, &win);
//     MPI_Win_shared_query(win, 0, &buf_size, &disp_unit, &rcv_buf);
//     //? send buffer intialization with values
//     int length = start_length;
//     descr_t descr;
//     descr.root = 0;
//     descr.message_length = length;
//     int mid = 2;
//     buf_dtype snd_buf[max_length];
//     snd_buf[0] = 4.2 + 1;
//     snd_buf[mid] = 4.2 + 2;
//     snd_buf[4 - 1] = 4.2 + 3;
//     if (my_rank == 0)
//     {
//         result = RMA_Bcast_Linear((buf_dtype *)snd_buf, MPI_FLOAT, buf_size, descr, size, win, comm_sm);
//         //! flush memory
//         MPI_Win_flush(my_rank, win);
//         MPI_Win_sync(win);
//         //! put result into file for root
//         filetestlinear << size << " processes" << std::endl;
//         filetestlinear << "rank 0 Bcast([" << snd_buf[0] << "," << snd_buf[mid] << "," << snd_buf[4 - 1] << "])" << std::endl;
//         filetestlinear << "rank " << my_rank << "-" << size << ":[" << rcv_buf[0] << "," << rcv_buf[mid] << "," << rcv_buf[4 - 1] << "]" << std::endl;
//     }
//     else
//     {
//         // usleep((my_rank + 1) * 50000);
//         //! flush memory
//         MPI_Win_flush(my_rank, win);
//         MPI_Win_sync(win);
//         //! put result into file for root
//         filetestlinear << "rank " << my_rank << "-" << size << ":[" << rcv_buf[0] << "," << rcv_buf[mid] << "," << rcv_buf[4 - 1] << "]" << std::endl;
//     }

//     MPI_Win_free(&win);
//     // MPI_Finalize();
//     EXPECT_EQ(result, MPI_SUCCESS);
// }