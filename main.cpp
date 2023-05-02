
#include <stdio.h>
#include "binomial_bcast.h"
#include "linear_bcast.h"
#include "binary_bcast.h"
#include "gtest/gtest.h"
#include "tests.cpp"
#define number_of_messages 50 /*number of messages transfered per package */
#define start_length 4        /*length of the array */
#define length_factor 8
#define max_length 8388608           /* ==> 2 x 32 MB per process */
#define number_package_sizes 8       /*number of package sizes for each size 16,128,...*/
/* #define max_length 67108864    */ /* ==> 2 x 0.5 GB per process */
/* #define number_package_sizes 9 */
// using buf_dtype = int; /*type of transfered message*/
enum bcast_types_t
{
    linear = 1,
    binomial = 2,
    binary = 3,
    test = 4
};
bcast_types_t bcast_type = linear;

int main(int argc, char *argv[])
{ // ? variable declaration
    std::string algname;
    int provided = MPI_THREAD_MULTIPLE;
    setbuf(stdout, NULL);
    // // Todo  create passing value for tybe of broadcast
    // ? MPI Intialization
    // boost::mpi::environment env(argc, argv);
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm comm_sm;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int result = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_sm);
    if (result != MPI_SUCCESS)
        MPI_Abort(comm_sm, result);

    //? choose the type of broadcast
    if (argc == 2)
    {
        if (std::string(argv[1]) == "linear")
            bcast_type = linear;
        else if (std::string(argv[1]) == "binomial")
            bcast_type = binomial;
        else if (std::string(argv[1]) == "binary")
            bcast_type = binary;
        else if (std::string(argv[1]) == "test")
        {
            filetestbinomial.open("results/resultTestBinomial" + std::to_string(size) + ".dat", std::ios::app); /*create file and open it*/
            filetestbinary.open("results/resultTestBinary" + std::to_string(size) + ".dat", std::ios::app);     /*create file and open it*/
            filetestlinear.open("results/resultTestLinear" + std::to_string(size) + ".dat", std::ios::app);     /*create file and open it*/
            ::testing::InitGoogleTest(&argc, argv);
            bcast_type = test;
            int res = RUN_ALL_TESTS();
        }
        else
            throw std::runtime_error("Invalid argument");
    }
    //? Shared memory allocate region
    MPI_Aint buf_size;
    MPI_Win win;
    int disp_unit;
    buf_dtype *rcv_buf; // rcv_buf pointer type
    result = MPI_Win_allocate_shared((MPI_Aint)max_length * sizeof(buf_dtype), sizeof(buf_dtype), MPI_INFO_NULL, comm_sm, &rcv_buf, &win);
    if (result != MPI_SUCCESS)
        MPI_Abort(comm_sm, result);

    MPI_Win_shared_query(win, my_rank, &buf_size, &disp_unit, &rcv_buf);

    //? File declaration
    std::fstream file; /* value for result file*/
    std::fstream fileBench;
    file.open("results/result" + std::string(argv[1]) + std::to_string(size) + ".dat", std::ios::out); /*create file and open it*/
    fileBench.open("results/result" + std::string(argv[1]) + std::to_string(size) + "Bench.dat", std::ios::out);
    if (my_rank == 0 && bcast_type != test)
    {
        printf("    message size      transfertime  duplex bandwidth per process and neighbor\n");
        fileBench << "    message size      transfertime  duplex bandwidth per process and neighbor" << std::endl;
    }
    //! test variables for Broadcast

    if (my_rank == 0 && bcast_type != test)
    {
        double start, finish, transfer_time;
        int i, mid, length, test_value;
        length = start_length;
        descr_t descr;
        descr.root = 0;
        buf_dtype snd_buf[max_length]; //[max_length];

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
                descr.message_length = length;
                // Todo make a generic method for each type to be compared
                if (bcast_type == binomial)
                    RMA_Bcast_binomial((buf_dtype *)snd_buf, rcv_buf, my_rank, descr, size, win, comm_sm);
                else if (bcast_type == linear)
                    RMA_Bcast_Linear((buf_dtype *)snd_buf, MPI_FLOAT, buf_size, descr, size, win, comm_sm);
                else if (bcast_type == binary)
                    BinaryTreeBcast((buf_dtype *)snd_buf, rcv_buf, my_rank, descr, size, win, comm_sm);
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
        MPI_Win_flush(my_rank, win);
    }
    else
    {
        if (bcast_type != test)
            MPI_Win_flush(my_rank, win);
    }

    MPI_Win_free(&win);
    MPI_Finalize();
}
