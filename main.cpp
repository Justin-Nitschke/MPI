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

float map_coord(float value, float min, float total_width) {
    return -1.0 + 2.0 * (value - min) / total_width;
}

void drawParticle2D(float x_sim, float y_sim, float radius, float r, float g, float b, float domain_min_x, float domain_min_y, float domain_width) {
    float x_gl = map_coord(x_sim, domain_min_x, domain_width);
    float y_gl = map_coord(y_sim, domain_min_y, domain_width);

    float rad_gl = radius * (2.0 / domain_width);

    int k = 0;
    float angle = 0.0;
    glBegin(GL_TRIANGLE_FAN);
    glColor3f(r, g, b);
    glVertex2f(x_gl, y_gl);
    for(k = 0; k <= 20; k++) {
        angle = (float)(k) / 20 * 2 * 3.141592;
        glVertex2f(x_gl + rad_gl * cos(angle),
                   y_gl + rad_gl * sin(angle));
    }
    glEnd();
}

void drawOctreeBounds2D(Node *node, float domain_min_x, float domain_min_y, float domain_width) {
    if (!node || !node->is_internal) return;

    float center_x_sim = node->square_x + node->half_width;
    float center_y_sim = node->square_y + node->half_width;

    float cx = map_coord(center_x_sim, domain_min_x, domain_width);
    float cy = map_coord(center_y_sim, domain_min_y, domain_width);

    float left = map_coord(node->square_x, domain_min_x, domain_width);
    float right = map_coord(node->square_x + 2 * node->half_width, domain_min_x, domain_width);
    float bottom = map_coord(node->square_y, domain_min_y, domain_width);
    float top = map_coord(node->square_y + 2 * node->half_width, domain_min_y, domain_width);

    glBegin(GL_LINES);
    glColor3f(1, 1, 1);

    glVertex2f(left, cy);
    glVertex2f(right, cy);

    glVertex2f(cx, bottom);
    glVertex2f(cx, top);
    glEnd();

    for(int i = 0; i < 4; i++) {
        drawOctreeBounds2D(node->children[i], domain_min_x, domain_min_y, domain_width);
    }
}

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

    /*********************************** VISUALIZATION ***********************************/
    /* OpenGL window dims */
    int width=600, height=600;
    GLFWwindow* window = nullptr;
    if (visualize && rank == 0) {
        if( !glfwInit() ){
            fprintf( stderr, "Failed to initialize GLFW\n" );
            return -1;
        }
        // Open a window and create its OpenGL context
        window = glfwCreateWindow( width, height, "Simulation", NULL, NULL);
        if( window == NULL ){
            fprintf( stderr, "Failed to open GLFW window.\n" );
            glfwTerminate();
            return -1;
        }
        glfwMakeContextCurrent(window); // Initialize GLEW
        if (glewInit() != GLEW_OK) {
            fprintf(stderr, "Failed to initialize GLEW\n");
            return -1;
        }
        // Ensure we can capture the escape key being pressed below
        glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    }

    /*********************************** SIMULATION LOOP ***********************************/
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    for (int t = 0; t < timesteps; t++) {
        //Set the domain to (0,0) to (4,4)
        float min_bound = 0;
        float max_bound = 4;
        Node *root = new Node(0, 0, 2);

        for(int i = 0; i < num_bodies; i++) {
            if (bodies[i].mass != -1) {
                bool out_of_bounds = 
                    bodies[i].x_position <= min_bound ||
                    bodies[i].x_position >= max_bound ||
                    bodies[i].y_position <= min_bound ||
                    bodies[i].y_position >= max_bound;

                if (out_of_bounds) {
                    bodies[i].mass = -1;
                }
            }
        }

        //Build tree
        for(int i = 0; i < num_bodies; i++) {
            bodies[i].x_force = 0;
            bodies[i].y_force = 0;

            if (bodies[i].mass != -1) {
                root->insert(&(bodies[i]));
            }
        }

        //Compute centers of mass
        root->center_of_mass();

        //calculate forces
        for(int i = start_index; i < end_index; i++) {
            if (bodies[i].mass != -1) {
                root->force(&bodies[i], theta);
            }
        }

        //move bodies
        for(int i = start_index; i < end_index; i++) {
            if (bodies[i].mass != -1) {
                float ax = bodies[i].x_force / bodies[i].mass;
                float ay = bodies[i].y_force / bodies[i].mass;

                bodies[i].x_position = bodies[i].x_position + bodies[i].x_velocity * dt + 0.5 * ax * dt * dt;
                bodies[i].y_position = bodies[i].y_position + bodies[i].y_velocity * dt + 0.5 * ay * dt * dt;

                bodies[i].x_velocity = bodies[i].x_velocity + ax * dt;
                bodies[i].y_velocity = bodies[i].y_velocity + ay * dt;
            }
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

        if (visualize && rank == 0) {
            float d_min_x = root->square_x;
            float d_min_y = root->square_y;
            float d_width = root->half_width * 2;

            glClear(GL_COLOR_BUFFER_BIT);

            drawOctreeBounds2D(root, d_min_x, d_min_y, d_width);

            for(int i = 0; i < num_bodies; i++) {
                if (bodies[i].mass != -1) {
                    float r = 0.2, g = 0.6, b = 1.0;
                    drawParticle2D(bodies[i].x_position, bodies[i].y_position, 0.05, r, g, b, d_min_x, d_min_y, d_width);
                }
            }

            glfwSwapBuffers(window);
            glfwPollEvents();
        }

        //cleanup
        delete root;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();

    if (rank == 0) {
        double elapsed = end_time - start_time;
        std::cout << elapsed << std::endl;

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

    if (visualize && rank == 0) {
        glfwTerminate();
    }

    MPI_Finalize();
    return 0;
}