
#include <stdio.h>
#include "binomial_bcast.h"
#include "linear_bcast.h"
#include "binary_bcast.h"
#include "gtest/gtest.h"

#define number_of_messages 50 /*number of messages transfered per package */
#define start_length 4        /*length of the array */
#define length_factor 8
#define max_length 8388608           /* ==> 2 x 32 MB per process */
#define number_package_sizes 8       /*number of package sizes for each size 16,128,...*/
/* #define max_length 67108864    */ /* ==> 2 x 0.5 GB per process */
/* #define number_package_sizes 9 */
using buf_dtype = float; /*type of transfered message*/
enum bcast_types_t
{
    linear = 1,
    binomial = 2,
    binary = 3
};
bcast_types_t bcast_type = linear;
int main(int argc, char *argv[])
{ // ? variable declaration
    int my_rank,
        size;
    std::string algname;
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

    setbuf(stdout, NULL);
    // Todo  create passing value for tybe of broadcast
    // ? MPI Intialization
    //::testing::InitGoogleTest(&argc, argv);
    // boost::mpi::environment env(argc, argv);

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    descr.root = my_rank;
    if (argc == 2)
    {
        if (std::string(argv[1]) == "linear")
            bcast_type = linear;
        else if (std::string(argv[1]) == "binomial")
            bcast_type = binomial;
        else if (std::string(argv[1]) == "binary")
            bcast_type = binary;
        else
            throw std::runtime_error("Invalid argument");
    }
    file.open("results/result" + std::string(argv[1]) + std::to_string(size) + ".dat", std::ios::out); /*create file and open it*/
    fileBench.open("results/result" + std::string(argv[1]) + std::to_string(size) + "Bench.dat", std::ios::out);
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
            if (bcast_type == binomial)
                RMA_Bcast_binomial((buf_dtype *)snd_buf, my_rank, i, descr, size, j, mid, length, file, win, comm_sm);
            else if (bcast_type == linear)
                RMA_Bcast_Linear((buf_dtype *)snd_buf, my_rank, i, size, j, mid, length, file, win, comm_sm);
            else if (bcast_type == binary)
                BinaryTreeBcast((buf_dtype *)snd_buf, my_rank, i, descr.root, size, j, mid, length, file, win, comm_sm);
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
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
    MPI_Finalize();
}
TEST(BinaryTreeBcast, BroadcastToOtherProcesses)
{
    int mid = (4 - 1) / number_of_messages * 1;
    int my_rank, size;
    MPI_Win win;
    MPI_Comm comm_sm;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::fstream file;                                                              /* value for result file*/
    file.open("results/resultTest" + std::to_string(size) + ".dat", std::ios::out); /*create file and open it*/

    float snd_buf[max_length];
    snd_buf[0] = 4.2 + 1;
    snd_buf[mid] = 4.2 + 2;
    snd_buf[4 - 1] = 4.2 + 3;
    int res = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    res = BinaryTreeBcast(snd_buf, my_rank, 1, 0, 1, 1, mid, 4, file, win, comm_sm);
    EXPECT_EQ(res, MPI_SUCCESS);
}
TEST(RMA_Bcast_binomial, BroadcastToOtherProcesses)
{
    int mid = (4 - 1) / number_of_messages * 1;
    int my_rank, size;
    MPI_Win win;
    descr_t descr;
    descr.root = 0;
    MPI_Comm comm_sm;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::fstream file;                                                              /* value for result file*/
    file.open("results/resultTest" + std::to_string(size) + ".dat", std::ios::out); /*create file and open it*/

    float snd_buf[max_length];
    snd_buf[0] = 4.2 + 1;
    snd_buf[mid] = 4.2 + 2;
    snd_buf[4 - 1] = 4.2 + 3;
    int res = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    res = RMA_Bcast_binomial(snd_buf, my_rank, 1, descr, 1, 1, mid, 4, file, win, comm_sm);
    EXPECT_EQ(res, MPI_SUCCESS);
}
TEST(RMA_Bcast_Linear, BroadcastToOtherProcesses)
{
    int mid = (4 - 1) / number_of_messages * 1;
    int my_rank, size;
    MPI_Win win;
    descr_t descr;
    descr.root = 0;
    MPI_Comm comm_sm;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::fstream file;                                                              /* value for result file*/
    file.open("results/resultTest" + std::to_string(size) + ".dat", std::ios::out); /*create file and open it*/

    float snd_buf[max_length];
    snd_buf[0] = 4.2 + 1;
    snd_buf[mid] = 4.2 + 2;
    snd_buf[4 - 1] = 4.2 + 3;
    int res = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    res = RMA_Bcast_Linear(snd_buf, my_rank, 1, size, 1, mid, 4, file, win, comm_sm);
    EXPECT_EQ(res, MPI_SUCCESS);
}