#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX_CMD_LEN 512
#define NUM_TESTS 12  // Number of parameter configurations

int main() {
    char galsim_cmd[MAX_CMD_LEN];
    char compare_cmd[MAX_CMD_LEN];
    char buffer[256];
    double execution_time, pos_maxdiff;
    
    double theta_values[] = {0, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.55, 0.4, 0.45,
                             0.5, 0.55, 0.6, 0.65, 0.7,
                             0.75, 0.85, 0.9, 0.95};  
                            

    double dt_values[] = {1e-5, 2e-5, 5e-5};
    printf("Starting Barnes-Hut Optimization Study...\n\n");

    FILE *file = fopen("optimization_results.csv", "w");
    if (!file) {
        fprintf(stderr, "Error: Could not open output file.\n");
        return 1;
    }

    fprintf(file, "theta_max,dt,nsteps,execution_time,pos_maxdiff\n");

    for (int t = 0; t < 21; t++) {
        for (int d = 0; d < 1; d++) {
            double theta = theta_values[t];
            double dt = dt_values[d];
            int nsteps = 200;  // Ensure fixed T = 1ms

            printf("Running θ = %.2f, Δt = %.1e, nsteps = %d...\n", theta, dt, nsteps);

            // Run simulation
            snprintf(galsim_cmd, MAX_CMD_LEN,
                     "/usr/bin/time -f '%%e' ./galsim 2000 input_data/ellipse_N_02000.gal %d %.1e %.2f 0 8 2> tmp_time.txt",
                     nsteps, dt, theta);

            int status = system(galsim_cmd);
            if (status != 0) {
                fprintf(stderr, "Error: galsim failed for θ = %.2f, Δt = %.1e\n", theta, dt);
                continue;
            }

            // Read execution time
            FILE *fp = fopen("tmp_time.txt", "r");
            if (!fp || !fgets(buffer, sizeof(buffer), fp)) {
                fprintf(stderr, "Error: Failed to read execution time.\n");
                continue;
            }
            sscanf(buffer, "%lf", &execution_time);
            fclose(fp);

            // Compare results with reference
            snprintf(compare_cmd, MAX_CMD_LEN, "./compare_gal_files 2000 result.gal reference_solution.gal > tmp_diff.txt");
            // snprintf(compare_cmd, MAX_CMD_LEN, "./compare_gal_files 2000 result.gal ref_output_data/ellipse_N_02000_after200steps.gal > tmp_diff.txt");

            system(compare_cmd);
            // Read position difference
            fp = fopen("tmp_diff.txt", "r");
            if (!fp) {
                fprintf(stderr, "Error: tmp_diff.txt could not be opened.\n");
                pos_maxdiff = -1.0;
            } else {
                while (fgets(buffer, sizeof(buffer), fp)) {
                    if (sscanf(buffer, "pos_maxdiff = %lf", &pos_maxdiff) == 1) {
                        break; // Found the correct line, exit loop
                    }
                }
                fclose(fp);
            }

            // Store results
            fprintf(file, "%.2f,%.1e,%d,%.5f,%.5e\n", theta, dt, nsteps, execution_time, pos_maxdiff);

            // Check if within tolerance
            if (pos_maxdiff > 0 && pos_maxdiff < 1e-2) {
                printf("✅ θ = %.2f, Δt = %.1e is within tolerance (%.5e).\n", theta, dt, pos_maxdiff);
            } else {
                printf("❌ θ = %.2f, Δt = %.1e EXCEEDS tolerance (%.5e).\n", theta, dt, pos_maxdiff);
            }

        }
    }

    fclose(file);
    printf("\nOptimization study complete. Results saved to optimization_results.csv\n");
    return 0;
}
