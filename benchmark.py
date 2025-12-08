import subprocess
import csv
import statistics
import time
import os
import math
import random

# ==========================================
# CONFIGURATION
# ==========================================

# Path to your C++ executable
EXECUTABLE = "./main"

# Parameters to test
# "N": Number of bodies (Script will generate temp input files for these)
BODY_COUNTS = [500, 1000, 2000] 

# "-n": Number of MPI Processes (mpiexec -n X)
MPI_RANKS = [1, 2, 4, 8]

# "-t": Theta (Accuracy)
THETAS = [0.2, 0.5, 0.8]

# "-s": Timesteps
TIMESTEPS = [20000]

# Fixed parameters
DT = 0.0005
USE_POOL = True # Add the -p flag?

# How many times to run each config to get an average
TRIALS = 5

# Output filename
CSV_FILENAME = "benchmark_results.csv"

# ==========================================
# INPUT GENERATOR (Galaxy Logic)
# ==========================================
def generate_galaxy_input(filename, num_bodies):
    """Generates a galaxy input file so we can test different N values."""
    with open(filename, 'w') as f:
        f.write(f"{num_bodies}\n")
        g1_center = (1.2, 2.0)
        g2_center = (2.8, 2.0)
        
        for i in range(num_bodies):
            if i < num_bodies // 2:
                cx, cy = g1_center
                vx_bulk, vy_bulk = 0.5, 0.2
                is_g1 = True
            else:
                cx, cy = g2_center
                vx_bulk, vy_bulk = -0.5, -0.2
                is_g1 = False
            
            angle = random.uniform(0, 2 * math.pi)
            dist = random.uniform(0.05, 0.6)
            
            x = cx + dist * math.cos(angle)
            y = cy + dist * math.sin(angle)
            
            # Orbital Velocity
            virtual_mass = 2.0
            velocity_magnitude = math.sqrt(virtual_mass / dist)
            vx_orb = -math.sin(angle) * velocity_magnitude
            vy_orb = math.cos(angle) * velocity_magnitude
            
            if not is_g1:
                vx_orb *= -1
                vy_orb *= -1
            
            vx = vx_bulk + vx_orb
            vy = vy_bulk + vy_orb
            mass = 0.01
            
            # Write line
            f.write(f"{i} {x:.6f} {y:.6f} {mass:.6f} {vx:.6f} {vy:.6f}\n")

# ==========================================
# BENCHMARK ENGINE
# ==========================================
def run_benchmark():
    # Prepare CSV
    with open(CSV_FILENAME, 'w', newline='') as csvfile:
        fieldnames = ['Num_Bodies', 'MPI_Ranks', 'Theta', 'Timesteps', 'Avg_Time_Sec']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        total_iterations = len(BODY_COUNTS) * len(MPI_RANKS) * len(THETAS) * len(TIMESTEPS)
        current_iter = 0

        print(f"Starting Benchmark... ({total_iterations} combinations)")

        for n_bodies in BODY_COUNTS:
            # Generate a temporary input file for this N
            input_file = f"temp_bench_{n_bodies}.txt"
            generate_galaxy_input(input_file, n_bodies)
            
            for ranks in MPI_RANKS:
                for theta in THETAS:
                    for steps in TIMESTEPS:
                        current_iter += 1
                        print(f"[{current_iter}/{total_iterations}] Running: N={n_bodies}, Ranks={ranks}, Theta={theta}, Steps={steps} ... ", end="", flush=True)
                        
                        run_times = []
                        
                        # Build Command
                        cmd = [
                            "mpiexec", "-n", str(ranks),
                            EXECUTABLE,
                            "-i", input_file,
                            "-o", "out/bench_trash.txt", # We don't care about output for timing
                            "-t", str(theta),
                            "-s", str(steps),
                            "-d", str(DT)
                        ]
                        
                        if USE_POOL:
                            cmd.append("-p")

                        # Run Trials
                        valid_run = True
                        for _ in range(TRIALS):
                            try:
                                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                                
                                # Parse the output (C++ prints only the float number now)
                                # We strip whitespace and convert to float
                                output_lines = result.stdout.strip().split('\n')
                                # Get the last line (in case MPI prints warnings)
                                time_str = output_lines[-1]
                                run_times.append(float(time_str))
                                
                            except subprocess.CalledProcessError as e:
                                print(f"FAILED (Exit Code {e.returncode})")
                                print(e.stderr)
                                valid_run = False
                                break
                            except ValueError:
                                print(f"FAILED (Could not parse output)")
                                print(f"Output was: {result.stdout}")
                                valid_run = False
                                break

                        if valid_run:
                            avg_time = statistics.mean(run_times)
                            print(f"Avg: {avg_time:.4f}s")
                            
                            writer.writerow({
                                'Num_Bodies': n_bodies,
                                'MPI_Ranks': ranks,
                                'Theta': theta,
                                'Timesteps': steps,
                                'Avg_Time_Sec': f"{avg_time:.6f}"
                            })
                            csvfile.flush() # Ensure data is written immediately

            # Cleanup temp file
            if os.path.exists(input_file):
                os.remove(input_file)

    print(f"\nBenchmark Complete! Results saved to {CSV_FILENAME}")

if __name__ == "__main__":
    # Ensure output directory exists
    if not os.path.exists("out"):
        os.makedirs("out")
        
    run_benchmark()