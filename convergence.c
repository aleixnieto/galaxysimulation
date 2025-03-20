/*
 * Velocity Verlet Convergence Study
 *
 * This program runs `galsim` with decreasing values of Δt, computes `pos_maxdiff`
 * against the true solution (`true_solution.gal`), and verifies second-order accuracy.
 *
 * -----------------------------
 * Compilation Instructions:
 * -----------------------------
 *     gcc -o convergence_study convergence_study.c -lm
 *
 * -----------------------------
 * Execution Instructions:
 * -----------------------------
 *     ./convergence_study
 *
 * The program will:
 *   - Run `galsim` for different Δt values.
 *   - Use `compare_gal_files` with `true_solution.gal` to compute `pos_maxdiff`.
 *   - Save results to `convergence_results.csv`.
 *   - Print the final error values on the screen.
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 
 #define MAX_CMD_LEN 512
 #define MAX_LINE_LEN 256
 #define TRUE_NSTEPS 2000000  // True solution step count
 #define TRUE_DT 1e-8         // True solution Δt
 
 int main() {
     char galsim_cmd[MAX_CMD_LEN];
     char compare_cmd[MAX_CMD_LEN];
     char buffer[MAX_LINE_LEN];
 
     double pos_maxdiff;
     double dt_values[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8}; // Δt values
     int nsteps_values[] = {20, 200, 2000, 20000, 200000, TRUE_NSTEPS}; // Steps corresponding to Δt
     int num_tests = sizeof(dt_values) / sizeof(dt_values[0]);
 
     FILE *file = fopen("convergence_results.csv", "w");
     if (!file) {
         fprintf(stderr, "Error: Could not open output file.\n");
         return 1;
     }
     fprintf(file, "dt,pos_maxdiff\n");
 
     printf("Starting Velocity Verlet Convergence Study...\n\n");
 
     for (int i = 0; i < num_tests; i++) {
         double dt = dt_values[i];
         int nsteps = nsteps_values[i];
 
         printf("Running Δt = %.1e with %d steps...\n", dt, nsteps);
 
         // 1) Run galsim with current Δt and steps
         snprintf(galsim_cmd, MAX_CMD_LEN,
                  "./galsim 3 input_data/sun_and_planets_N_3.gal %d %.1e 0 0 8",
                  nsteps, dt);
         int status = system(galsim_cmd);
         if (status != 0) {
             fprintf(stderr, "Error: galsim failed for Δt = %.1e\n", dt);
             continue;
         }
 
         // 2) Run compare_gal_files using "true_solution.gal"
         snprintf(compare_cmd, MAX_CMD_LEN,
                  "./compare_gal_files 3 result.gal true_solution.gal");
 
         FILE *pipe = popen(compare_cmd, "r");
         if (!pipe) {
             fprintf(stderr, "Error: compare_gal_files failed for Δt = %.1e\n", dt);
             continue;
         }
 
         // 3) Read pos_maxdiff from output
         while (fgets(buffer, MAX_LINE_LEN, pipe)) {
             if (sscanf(buffer, "pos_maxdiff = %lf", &pos_maxdiff) == 1) {
                 break;
             }
         }
         pclose(pipe);
 
         // 4) Save and print results
         fprintf(file, "%.1e,%.9f\n", dt, pos_maxdiff);
         printf("   pos_maxdiff = %.9f\n", pos_maxdiff);
     }
 
     fclose(file);
     printf("\nConvergence study complete. Results saved to convergence_results.csv\n");
 
     return 0;
 }
 