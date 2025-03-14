#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "graphics.h"

#define EPSILON 1e-3

typedef struct {
    double *pos_x;      // x-coordinates of all particles
    double *pos_y;      // y-coordinates of all particles
    double *mass;       // mass of all particles
    double *vel_x;      // x-velocities of all particles
    double *vel_y;      // y-velocities of all particles
    double *brightness; // brightness of all particles
} ParticleData;

double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

void allocateParticleData(ParticleData *p, int N) {
    p->pos_x      = (double *)malloc(N * sizeof(double));
    p->pos_y      = (double *)malloc(N * sizeof(double));
    p->mass       = (double *)malloc(N * sizeof(double));
    p->vel_x      = (double *)malloc(N * sizeof(double));
    p->vel_y      = (double *)malloc(N * sizeof(double));
    p->brightness = (double *)malloc(N * sizeof(double));
}

void freeParticleData(ParticleData *p) {
    free(p->pos_x);
    free(p->pos_y);
    free(p->mass);
    free(p->vel_x);
    free(p->vel_y);
    free(p->brightness);
}

/* Read data from binary file into the SoA structure */
int readParticleData(const char *filename, ParticleData *p, int N) {
    FILE *file = fopen(filename, "rb");
    if (!file) {
        fprintf(stderr, "Error: Unable to open file %s\n", filename);
        return 1;
    }

    for (int i = 0; i < N; i++) {
        if (fread(&p->pos_x[i],      sizeof(double), 1, file) != 1 ||
            fread(&p->pos_y[i],      sizeof(double), 1, file) != 1 ||
            fread(&p->mass[i],       sizeof(double), 1, file) != 1 ||
            fread(&p->vel_x[i],      sizeof(double), 1, file) != 1 ||
            fread(&p->vel_y[i],      sizeof(double), 1, file) != 1 ||
            fread(&p->brightness[i], sizeof(double), 1, file) != 1)
        {
            fprintf(stderr, "Error: reading particle %d from file %s\n", i, filename);
            fclose(file);
            return 1;
        }
    }
    fclose(file);
    return 0;
}

/* Write final data to "result.gal" from the SoA structure */
void writeParticleData(const char *filename, ParticleData *p, int N) {
    FILE *binary_file = fopen(filename, "wb");
    if (!binary_file) {
        fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
        return;
    }

    for (int i = 0; i < N; i++) {
        fwrite(&p->pos_x[i],      sizeof(double), 1, binary_file);
        fwrite(&p->pos_y[i],      sizeof(double), 1, binary_file);
        fwrite(&p->mass[i],       sizeof(double), 1, binary_file);
        fwrite(&p->vel_x[i],      sizeof(double), 1, binary_file);
        fwrite(&p->vel_y[i],      sizeof(double), 1, binary_file);
        fwrite(&p->brightness[i], sizeof(double), 1, binary_file);
    }
    fclose(binary_file);
}

void updateForce(ParticleData *p, int N, double G, double time_step) {
    double *force_x = (double *)calloc(N, sizeof(double));
    double *force_y = (double *)calloc(N, sizeof(double));

    for (int i = 0; i < N; i++) {
        double x_i    = p->pos_x[i];
        double y_i    = p->pos_y[i];
        double mass_i = p->mass[i];

        double fxi = 0.0, fyi = 0.0;  // Accumulate force contributions for particle i

        for (int j = i + 1; j < N; j++) {
            double dx = x_i - p->pos_x[j];
            double dy = y_i - p->pos_y[j];

            double r = sqrt(dx * dx + dy * dy) + EPSILON;
            double force = (mass_i * p->mass[j]) / (r * r * r);

            fxi -= force * dx;
            fyi -= force * dy;

            // Update forces for particle j
            force_x[j] += force * dx;
            force_y[j] += force * dy;
        }

        // Store accumulated forces for particle i
        force_x[i] += fxi;
        force_y[i] += fyi;
    }

    // Update velocities and positions (symplectic Euler)
    for (int i = 0; i < N; i++) {
        p->vel_x[i] += ((G * force_x[i]) / p->mass[i]) * time_step;
        p->vel_y[i] += ((G * force_y[i]) / p->mass[i]) * time_step;

        p->pos_x[i] += p->vel_x[i] * time_step;
        p->pos_y[i] += p->vel_y[i] * time_step;
    }

    free(force_x);
    free(force_y);
}


int main(int argc, char* argv[])
{
    if (argc != 6) {
        printf("Usage: %s N input_file nsteps delta_t graphics\n", argv[0]);
        return 0;
    }

    int    N_bodies = atoi(argv[1]);
    char  *filename = argv[2];
    int    N_iter   = atoi(argv[3]);
    double delta_t  = atof(argv[4]);
    int    graphics = atoi(argv[5]);

    ParticleData particles;
    allocateParticleData(&particles, N_bodies);

    // Read initial data
    if (readParticleData(filename, &particles, N_bodies) != 0) {
        freeParticleData(&particles);
        return 1;
    }

    // Gravitational constant scaled by 1/N
    double G = 100.0 / (double)N_bodies;

    // Initialize graphics if needed
    if (graphics) {
        InitializeGraphics(argv[0]);
    }

    double start_time = get_wall_seconds();

    // Main simulation loop
    for (int step = 0; step < N_iter; step++) {
        updateForce(&particles, N_bodies, G, delta_t);

        // If graphics is on, draw the current state
        if (graphics) {
            ClearScreen();
            for (int i = 0; i < N_bodies; i++) {
                DrawCircle((float)particles.pos_x[i], (float)particles.pos_y[i]);
            }
            RefreshDisplay();
            usleep(3000);
            CheckForQuit();
            if (quit_flag) break;
        }
    }

    double end_time = get_wall_seconds();
    printf("Simulation completed in %f seconds.\n", end_time - start_time);

    // Write final state to file
    writeParticleData("result.gal", &particles, N_bodies);

    if (graphics) {
        CloseGraphics();
    }
    freeParticleData(&particles);

    return 0;
}