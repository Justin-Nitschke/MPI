#include <mpi.h>
#include <iostream>

#include "Body.hpp"
#include "Node.hpp"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int num_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    std::cout << "My rank is " << rank << "\n";
    std::cout << "Num processes is " << num_processes << "\n";
}