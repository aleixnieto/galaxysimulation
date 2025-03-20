/*
 * θ_max Tuning for Barnes-Hut Simulation
 *
 * This program was generated with the help of ChatGPT.
 *
 * Prompt used:
 * "Create a C program that runs `galsim` with increasing values of θ_max from 0.02 to 0.5 
 * (incrementing by 0.01) and measures `pos_maxdiff` using `compare_gal_files`. 
 * The program should find the largest θ_max where pos_maxdiff < 1e-3 and print the results. 
 * It should not write to a temporary file but directly display results on screen. 
 *
 * To compile:
 *     gcc -o theta_tuning theta_tuning.c -lm
 *
 * Run:
 *     ./theta_tuning
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 
 #define MAX_CMD_LEN 512
 #define MAX_LINE_LEN 256
 #define THETA_START 0.02
 #define THETA_END 0.5
 #define THETA_STEP 0.01
 #define THRESHOLD 1e-3
 
 int main() {
     char galsim_cmd[MAX_CMD_LEN];
     char compare_cmd[MAX_CMD_LEN];
 
     double best_theta = -1;
     double pos_maxdiff = 0.0;
 
     printf("Starting θ_max tuning...\n\n");
 
     // Loop over θ_max values from 0.02 to 0.5, increasing by 0.01
     for (double theta = THETA_START; theta <= THETA_END; theta += THETA_STEP) {
         printf("Testing θ_max = %.3f\n", theta);
 
         // 1) Run galsim with current θ_max
         snprintf(galsim_cmd, MAX_CMD_LEN,
                  "./galsim 2000 input_data/ellipse_N_02000.gal 200 1e-5 %.3f 0 8",
                  theta);
         int status = system(galsim_cmd);
         if (status != 0) {
             fprintf(stderr, "Error: galsim failed for θ_max = %.3f\n", theta);
             continue;
         }
 
         // 2) Run compare_gal_files and capture output
         snprintf(compare_cmd, MAX_CMD_LEN,
                  "./compare_gal_files 2000 result.gal ref_output_data/ellipse_N_02000_after200steps.gal");
 
         FILE *pipe = popen(compare_cmd, "r");
         if (!pipe) {
             fprintf(stderr, "Error: Failed to execute compare_gal_files for θ_max = %.3f\n", theta);
             continue;
         }
 
         // 3) Read and extract pos_maxdiff directly from output
         char buffer[MAX_LINE_LEN];
         while (fgets(buffer, MAX_LINE_LEN, pipe)) {
             if (sscanf(buffer, "pos_maxdiff = %lf", &pos_maxdiff) == 1) {
                 break;
             }
         }
         pclose(pipe);
 
         // 4) Print and check threshold
         printf("   pos_maxdiff = %.9f\n", pos_maxdiff);
         if (pos_maxdiff < THRESHOLD) {
             best_theta = theta;
         } else {
             printf("   -> pos_maxdiff exceeded threshold! Stopping.\n");
             break;
         }
     }
 
     // Print final result
     if (best_theta > 0) {
         printf("\nOptimal θ_max for pos_maxdiff < 1e-3: θ_max = %.3f\n", best_theta);
     } else {
         printf("\nNo suitable θ_max found within the given range.\n");
     }
 
     return 0;
 }
 