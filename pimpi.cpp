/*
Pi by Monte Carlo (MPI edition)
Example code by R. K. Lindsey (2024)

Generated via UM-GPT based on my pi_by_MC.py code, with the following prompt, in the same chat used to request the C++ translation from python:

Please re-write this code so that it can run via MPI. Please keep it simple and explain all the additions necessary to make it MPI parallel in the comments.

Compile code with: mpic++ -o pi_by_MC_mpi pi_by_MC_mpi.cpp -o pi_by_MC-mpi

Run code with: mpirun -np <num_processes> ./pi_by_MC-mpi <nsamples>

*/

// Include necessary header files

#include <iostream>      // For input and output functions
#include <cstdlib>       // For standard library functions like atoi
#include <cmath>         // For mathematical functions like pow
#include <random>        // For random number generation
#include <mpi.h>         // For MPI functions

int main(int argc, char *argv[]) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Check to make sure the user has entered the number of samples
    if (argc != 2) {
        if (world_rank == 0) { // Only the master node should print the error
            std::cerr << "Usage: " << argv[0] << " <nsamples>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    // Convert the command-line argument from string to integer
    int nsamples = std::atoi(argv[1]);

    // Calculate the number of samples per process
    int samples_per_process = nsamples / world_size;

    // Initialize the number of samples that fall within the circle
    int ncircle_local = 0;

    // Create a random number generator unique to each process
    std::random_device rd;
    std::mt19937 gen(rd() + world_rank);
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Loop over the number of samples for this process
    for (int i = 0; i < samples_per_process; ++i) {
        double x = dis(gen);
        double y = dis(gen);

        if (std::pow(x, 2.0) + std::pow(y, 2.0) <= 1.0) {
            ncircle_local += 1;
        }
    }

    // Reduce all of the local ncircle values to a global ncircle value
    int ncircle_global = 0;
    MPI_Reduce(&ncircle_local, &ncircle_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Only the master process should do the final calculation and print the result
    if (world_rank == 0) {
        double pi_estimate = 4.0 * (double(ncircle_global) / double(nsamples));
        std::cout << pi_estimate << std::endl;
    }

    // Finalize the MPI environment
    MPI_Finalize();
    return 0;
}
