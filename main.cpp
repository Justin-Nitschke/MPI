#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <GL/glu.h>
#include <GL/glut.h>

#include "Node.hpp"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int num_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    std::string input_file = "";
    std::string output_file = "";
    int timesteps = 100;
    float theta = 0.5;
    float dt = 0.005;
    bool visualize = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-i" || arg == "-I") {
            if (i + 1 < argc) {
                input_file = argv[i + 1];
                i++;
            } else {
                if (rank == 0) {
                    std::cerr << "Error: " << arg << " requires an input filename argument." << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
        } else if (arg == "-s" || arg == "-S") {
            if (i + 1 < argc) {
                timesteps = atoi(argv[i + 1]);
                i++;
            } else {
                if (rank == 0) {
                    std::cerr << "Error: " << arg << " requires a number of timesteps argument." << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
        } else if (arg == "-t" || arg == "-T") {
            if (i + 1 < argc) {
                theta = atof(argv[i + 1]);
                i++;
            } else {
                if (rank == 0) {
                    std::cerr << "Error: " << arg << " requires a theta argument." << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
        } else if (arg == "-d" || arg == "-D") {
            if (i + 1 < argc) {
                dt = atof(argv[i + 1]);
                i++;
            } else {
                if (rank == 0) {
                    std::cerr << "Error: " << arg << " requires a timestep argument." << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
        } else if (arg == "-o" || arg == "-O") {
            if (i + 1 < argc) {
                output_file = argv[i + 1];
                i++;
            } else {
                if (rank == 0) {
                    std::cerr << "Error: " << arg << " requires an output filename argument." << std::endl;
                }
                MPI_Finalize();
                return 1;
            }
        } else if (arg == "-v" || arg == "-V") {
            visualize = true;
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

    std::vector<int> recvcounts(num_processes);
    std::vector<int> displs(num_processes);

    int remainder = num_bodies % num_processes;
    int sum = 0;

    for (int i = 0; i < num_processes; i++) {
        int count = num_bodies / num_processes;
        if (i < remainder) {
            count++;
        }

        recvcounts[i] = count;
        displs[i] = sum;
        sum += count;
    }

    int start_index = displs[rank];
    int end_index = start_index + recvcounts[rank];
    int local_count = recvcounts[rank];

    std::vector<float> local_x(local_count);
    std::vector<float> local_y(local_count);
    std::vector<float> local_vx(local_count);
    std::vector<float> local_vy(local_count);

    std::vector<float> global_x(num_bodies);
    std::vector<float> global_y(num_bodies);
    std::vector<float> global_vx(num_bodies);
    std::vector<float> global_vy(num_bodies);

    /*********************************** SIMULATION LOOP ***********************************/
    for (int t = 0; t < timesteps; t++) {
        //Set the domain to (0,0) to (4,4)
        Node *root = new Node(0, 0, 2);

        //Build tree
        for(int i = 0; i < num_bodies; i++) {
            root->insert(&(bodies[i]));
            bodies[i].x_force = 0;
            bodies[i].y_force = 0;
        }

        //Compute centers of mass
        root->center_of_mass();

        //calculate forces
        for(int i = start_index; i < end_index; i++) {
            root->force(&bodies[i], theta);
        }

        //move bodies
        for(int i = start_index; i < end_index; i++) {
            float ax = bodies[i].x_force / bodies[i].mass;
            float ay = bodies[i].y_force / bodies[i].mass;

            bodies[i].x_position = bodies[i].x_position + bodies[i].x_velocity * dt + 0.5 * ax * dt * dt;
            bodies[i].y_position = bodies[i].y_position + bodies[i].y_velocity * dt + 0.5 * ay * dt * dt;

            bodies[i].x_velocity = bodies[i].x_velocity + ax * dt;
            bodies[i].y_velocity = bodies[i].y_velocity + ay * dt;
        }

        for(int i = 0; i < local_count; i++) {
            local_x[i] = bodies[start_index + i].x_position;
            local_y[i] = bodies[start_index + i].y_position;
            local_vx[i] = bodies[start_index + i].x_velocity;
            local_vy[i] = bodies[start_index + i].y_velocity;
        }

        MPI_Allgatherv(
            local_x.data(),
            local_count,
            MPI_FLOAT,
            global_x.data(),
            recvcounts.data(),
            displs.data(),
            MPI_FLOAT,
            MPI_COMM_WORLD
        );

        MPI_Allgatherv(
            local_y.data(),
            local_count,
            MPI_FLOAT,
            global_y.data(),
            recvcounts.data(),
            displs.data(),
            MPI_FLOAT,
            MPI_COMM_WORLD
        );

        MPI_Allgatherv(
            local_vx.data(),
            local_count,
            MPI_FLOAT,
            global_vx.data(),
            recvcounts.data(),
            displs.data(),
            MPI_FLOAT,
            MPI_COMM_WORLD
        );

        MPI_Allgatherv(
            local_vy.data(),
            local_count,
            MPI_FLOAT,
            global_vy.data(),
            recvcounts.data(),
            displs.data(),
            MPI_FLOAT,
            MPI_COMM_WORLD
        );

        for(int i = 0; i < num_bodies; i++) {
            bodies[i].x_position = global_x[i];
            bodies[i].y_position = global_y[i];
            bodies[i].x_velocity = global_vx[i];
            bodies[i].y_velocity = global_vy[i];
        }

        //cleanup
        delete root;
    }

    if (rank == 0) {
        std::ofstream outfile(output_file);

        if(!outfile.is_open()) {
            std::cerr << "Error: could not open output file" << output_file << std::endl;
        } else {
            outfile << num_bodies << "\n";
            outfile << std::fixed << std::setprecision(6);
            
            for(int i = 0; i < num_bodies; i++) {
                outfile << bodies[i].index << " "
                        << bodies[i].x_position << " "
                        << bodies[i].y_position << " "
                        << bodies[i].mass << " "
                        << bodies[i].x_velocity << " "
                        << bodies[i].y_velocity << "\n";
            }

            outfile.close();
        }
    }

    MPI_Finalize();
    return 0;
}