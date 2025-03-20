/*
 * Barnes-Hut Scaling Study
 *
 * This program automates the process of evaluating how the execution time of
 * the Barnes-Hut algorithm scales with the number of particles N.
 *
 * The simulation is run for different values of N while keeping:
 *   - θ_max = 0.21
 *   - nsteps = 200
 *   - delta_t = 1e-5
 *   - graphics = 0 (off)
 *   - nthreads = 8 (for parallelization)
 *
 * The execution time is recorded and written to `scaling_results.csv`, allowing
 * further analysis to determine how the performance scales with N.
 *
 * To compile:
 *     gcc -o scaling_study scaling_study.c -lm
 *
 * Run:
 *     ./scaling_study
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX_CMD_LEN 512
#define NUM_TESTS 10  // Number of different N values to test

int main() {
    char galsim_cmd[MAX_CMD_LEN];
    char buffer[256];
    double execution_time;
    int N_values[NUM_TESTS] = {1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000};
    double times[NUM_TESTS];

    printf("Starting Barnes-Hut Scaling Study (θ_max = 0.21)...\n\n");

    FILE *file = fopen("scaling_results.csv", "w");
    if (!file) {
        fprintf(stderr, "Error: Could not open output file.\n");
        return 1;
    }

    for (int i = 0; i < NUM_TESTS; i++) {
        int N = N_values[i];
        printf("Running N = %d...\n", N);

        // Build command to execute galsim and measure time
        snprintf(galsim_cmd, MAX_CMD_LEN,
                 "/usr/bin/time -f '%%e' ./galsim %d input_data/ellipse_N_%05d.gal 200 1e-5 0.21 0 8 2> tmp_time.txt",
                 N, N);
        
        int status = system(galsim_cmd);
        if (status != 0) {
            fprintf(stderr, "Error: galsim failed for N = %d\n", N);
            continue;
        }

        // Read execution time from tmp_time.txt
        FILE *fp = fopen("tmp_time.txt", "r");
        if (!fp) {
            fprintf(stderr, "Error: Failed to read execution time.\n");
            continue;
        }
        if (fgets(buffer, sizeof(buffer), fp)) {
            sscanf(buffer, "%lf", &execution_time);
            times[i] = execution_time;
            fprintf(file, "%d,%.5f\n", N, execution_time);
        }
        fclose(fp);
    }

    fclose(file);
    printf("\nScaling study complete. Results saved to scaling_results.csv\n");
    return 0;
}
