#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "Body.hpp"
#include "Node.hpp"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int num_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    std::string input_file = "";

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-i" || arg == "-I") {
            if (i + 1 < argc) {
                input_file = argv[i + 1];
                i++;
            } else {
                if (rank == 0) {
                    std::cerr << "Error: " << arg << " requires a filename argument." << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
        }
    }

    std::ifstream infile(input_file);

    if (!infile.is_open()) {
        if (rank == 0) {
            std::cerr << "Error: Can't open input file " << input_file << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    int num_bodies;
    infile >> num_bodies;

    std::vector<Body> bodies(num_bodies);

    for (int i = 0; i < num_bodies; i++) {
        infile >> bodies[i].index
               >> bodies[i].x_position
               >> bodies[i].y_position
               >> bodies[i].mass
               >> bodies[i].x_velocity
               >> bodies[i].y_velocity;
    }
    
    MPI_Finalize();
}