/**
 * Student Name: Şükrü Can Gültop
 * Student Number: 2014400201
 * Compile Status: Compiling
 * Program Status: Working
 * Notes: Please compile the code with mpiCC not mpicc because of c++(I also tested with mpiCC).
 */

#include <mpi.h>
#include <fstream>
#include <cmath>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <iterator>

using namespace std;

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Find out rank, size
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    //the number of slave processes.
    int N = world_size-1;
    ifstream inputFile(argv[1]);        // Input file stream object

    //if process is master process. It reads the input file broadcasts subarrays to slaves and recieves finished subarrays and joins them.
    if (world_rank == 0) {
        //number of slave processes check.
        if(200%N!=0){
            cout << "please give proper numbers for 200 pixels to divide equally." << endl;
            exit(EXIT_FAILURE);
        }
        //input array.
        int X[200][200];
        //reading the file and assigning corresponding array element.
        int current_number = 0;
        for (int i = 0; i < 200; i++) {
            for (int j = 0; j < 200; ++j) {
                inputFile >> current_number;
                X[i][j] = (current_number);
            }

        }

        //broadcasting to slave processes.
        for(int i = 1 ; i <= N ; i++)
            MPI_Send(X[(200/N)*(i-1)], 200*200/N, MPI_INT, i, 0, MPI_COMM_WORLD);

        // Close the file.
        inputFile.close();
        // output array after processing.
        int outputArr[200][200];
        //recieves subarrays and join them according to their process ids.
        for (int k = 1; k <= N; ++k) {
            int* subarr = NULL;
            subarr = (int *)malloc(sizeof(int) * 200*200/N);
            MPI_Recv(subarr, 200*200/N, MPI_INT, k, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < 200 / N; ++i) {
                for (int j = 0; j < 200; ++j) {
                    outputArr[(k-1)*(200/N)+i][j] = subarr[i*200+j];
                }
            }
        }
        // output file for the final array(denoised image).
        ofstream outputfile;
        outputfile.open(argv[2]);
        // writing to output file.
        for (int n = 0; n < 200; ++n) {
            for (int i = 0; i < 200; ++i) {
                outputfile << outputArr[n][i] << " ";
            }
            outputfile << endl;
        }
        outputfile.close();
    } else{
        //subarray for receiving from master process.
        int* subarr = NULL;
        //this array is one row but it corresponds to whole subarray just flatted.
        subarr = (int *)malloc(sizeof(int) * 200*200/N);
        //gets the corresponding subarray from master process.
        MPI_Recv(subarr, 200*200/N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //deflated subarray.
        int subArray[200/N][200];
        //deflating received array to double dimentional array to work with.
        for (int i = 0; i < 200 / N; ++i) {
            for (int j = 0; j < 200; ++j) {
                subArray[i][j] = subarr[i*200+j];
            }
        }
        //mathematical constants we get from user to use in formula.
        double beta = stod(argv[3]);
        double pi = stod(argv[4]);
        double gamma = 0.5 * log((1 - pi) / pi);
        //Z array which has processed pixels and it is copy of X(subarray) at first.
        int Z[200/N][200];
        copy(&subArray[0][0], &subArray[0][0] + 200/N * 200, &Z[0][0]);
        //iteration number.
        int T = 1000000/N;
        for (int k = 0; k < T; ++k) {
            //upper row from upper process' array.
            int* upperRow = NULL;
            //lower row from lower process' array
            int* lowerRow = NULL;
            //they are only one rows so it is enough only 200.
            upperRow = (int *)malloc(sizeof(int) * 200);
            lowerRow = (int *)malloc(sizeof(int) * 200);

            //messaging part
            //if the process is not first it can get upper row from previous process. And it check whether there is only master process or not. If only master no messaging. Also tag is iteration number so there wont be data race.
            if (world_rank != 1){
                MPI_Recv(upperRow, 200, MPI_INT, world_rank-1, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else if(world_size > 2){
                MPI_Send(subArray[200/N-1], 200, MPI_INT, world_rank+1, k, MPI_COMM_WORLD);
            }
            //If process is not last and first it can send to next process.
            if (world_rank != N && world_rank != 1) {
                MPI_Send(subArray[200/N-1], 200, MPI_INT, world_rank+1, k, MPI_COMM_WORLD);
            }
            //If process is not last, it can recieve from next process to lower row..
            if (world_rank != N){
                MPI_Recv(lowerRow, 200, MPI_INT, world_rank+1, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }else if(world_size > 2) {
                MPI_Send(subArray[0], 200, MPI_INT, world_rank-1, k, MPI_COMM_WORLD);
            }
            //If process is not last, it can send to previous process to upper row..
            if (world_rank != N && world_rank != 1){
                MPI_Send(subArray[0], 200, MPI_INT, world_rank-1, k, MPI_COMM_WORLD);
            }

            //choosing random pixel
            int i = rand() % (200/N);
            int j = rand() % 200;


            //sum of all neighbours and itself(the pixel itself will be substracted in formula.).
            int sum = 0;
            for (int l = max(i - 1, 0); l < min(i + 2, 200/N); ++l) {
                for (int m = max(j - 1, 0); m < min(j + 2, 200); ++m) {
                    sum += Z[l][m];
                }
            }
            //if choosen pixel is an bound pixel it checks upper and lower rows for sum.
            if(world_rank!=1 && i == 0){
                for (int l = max(j-1,0); l < min(j+2,200); ++l) {
                    sum += upperRow[l];
                }
            }
            if(world_rank!=N && i == 200/N-1){
                for (int l = max(j-1,0); l < min(j+2, 200); ++l) {
                    sum += lowerRow[l];
                }
            }

            //calculating formula and reverting the pixel.
            double delta_E = -2 * gamma * subArray[i][j] * Z[i][j] - 2 * beta * Z[i][j] * (sum - Z[i][j]);
            double prob = log((double) rand() / RAND_MAX);
            if (prob < delta_E) {
                Z[i][j] = -Z[i][j];
            }
        }
        //sends finalized subarray to masterprocess.
        MPI_Send(Z[0], 200*200/N, MPI_INT, 0, world_rank, MPI_COMM_WORLD);
    }
    MPI_Finalize();
}