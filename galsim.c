/*
 * BARNES-HUT N-BODY SIMULATION WITH OPENMP PARALLELIZATION
 *
 * This program simulates an N-body system using the Barnes-Hut algorithm, 
 * which improves efficiency by approximating long-range gravitational interactions.
 *
 * The simulation organizes particles into a quadtree to reduce the computational 
 * complexity of force calculations from O(N²) to approximately O(N log N). 
 * It supports both simple explicit Euler integration and the more accurate 
 * Velocity Verlet method for updating particle positions. The code also takes 
 * advantage of OpenMP to parallelize force calculations and particle updates, 
 * allowing for faster execution on multi-core systems.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include "graphics.h"

#define EPSILON 1e-3
#define VERLET 0  // Set to 1 for Verlet integration, 0 for simple integration
double THETA_MAX;

typedef struct {
    double *pos_x;      // x-coordinates of all particles
    double *pos_y;      // y-coordinates of all particles
    double *mass;       // mass of all particles
    double *vel_x;      // x-velocities of all particles
    double *vel_y;      // y-velocities of all particles
    double *brightness; // brightness of all particles
} ParticleData;

// Quadtree node structure 
typedef struct QuadNode {
    struct QuadNode *child[4];

    // Bounding box for this node
    double x_min, x_max;
    double y_min, y_max;

    // Accumulated mass and center of mass
    double mass;
    double center_x, center_y;

    // If  it is a leaf, westore index of a single particle (or -1 if internal node)
    int particle_index;
} QuadNode;

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

// Quadtree
QuadNode* createNode(double x_min, double x_max, double y_min, double y_max) {
    QuadNode* node = (QuadNode*)calloc(1, sizeof(QuadNode));
    node->x_min = x_min;  node->x_max = x_max;
    node->y_min = y_min;  node->y_max = y_max;
    node->mass  = 0.0;
    node->center_x = 0.0;
    node->center_y = 0.0;
    node->particle_index = -1;

    for (int i = 0; i < 4; i++) {
        node->child[i] = NULL;
    }
    return node;
}

// Determine which child quadrant (0,1,2,3) a particle belongs to
int getQuadrant(QuadNode *node, double px, double py) {
    double mid_x = 0.5 * (node->x_min + node->x_max);
    double mid_y = 0.5 * (node->y_min + node->y_max);
    int quadrant = 0;
    // child[0] = bottom-left, [1] = bottom-right, [2] = top-left, [3] = top-right
    int left  = (px < mid_x);
    int bot   = (py < mid_y);

    if ( left &&  bot ) quadrant = 0; // BL
    if (!left &&  bot ) quadrant = 1; // BR
    if ( left && !bot ) quadrant = 2; // TL
    if (!left && !bot ) quadrant = 3; // TR
    return quadrant;
}

// We insert a particle index into the quadtree node recursively
void insertParticle(QuadNode *node, ParticleData *p, int idx) {
    // If the node is empty, i.e., it is a leaf with no particle
    if (node->particle_index == -1 && 
        !node->child[0] && !node->child[1] && !node->child[2] && !node->child[3]) 
    {
        node->particle_index = idx;
        return;
    }

    // If the node already holds one particle, we have to split
    if (node->particle_index != -1) {
        int existing_idx = node->particle_index;
        node->particle_index = -1;  // Now the node becomes an internal node

        // Insert the old particle into one of the children
        double px_old = p->pos_x[existing_idx];
        double py_old = p->pos_y[existing_idx];
        int q_old = getQuadrant(node, px_old, py_old);

        // If that child does no exist, we create it
        if (!node->child[q_old]) {
            // THe bounding box for that child
            double mid_x = 0.5 * (node->x_min + node->x_max);
            double mid_y = 0.5 * (node->y_min + node->y_max);
            switch (q_old) {
                case 0: // BL
                    node->child[0] = createNode(node->x_min, mid_x, node->y_min, mid_y);
                    break;
                case 1: // BR
                    node->child[1] = createNode(mid_x, node->x_max, node->y_min, mid_y);
                    break;
                case 2: // TL
                    node->child[2] = createNode(node->x_min, mid_x, mid_y, node->y_max);
                    break;
                case 3: // TR
                    node->child[3] = createNode(mid_x, node->x_max, mid_y, node->y_max);
                    break;
            }
        }
        insertParticle(node->child[q_old], p, existing_idx);
    }

    // Now insert the new particle
    double px = p->pos_x[idx];
    double py = p->pos_y[idx];
    int q_new = getQuadrant(node, px, py);

    if (!node->child[q_new]) {
        double mid_x = 0.5 * (node->x_min + node->x_max);
        double mid_y = 0.5 * (node->y_min + node->y_max);
        switch (q_new) {
            case 0: // BL
                node->child[0] = createNode(node->x_min, mid_x, node->y_min, mid_y);
                break;
            case 1: // BR
                node->child[1] = createNode(mid_x, node->x_max, node->y_min, mid_y);
                break;
            case 2: // TL
                node->child[2] = createNode(node->x_min, mid_x, mid_y, node->y_max);
                break;
            case 3: // TR
                node->child[3] = createNode(mid_x, node->x_max, mid_y, node->y_max);
                break;
        }
    }
    insertParticle(node->child[q_new], p, idx);
}

// Build the entire quadtree for all particles. We assume the bounding box is [0,1] x [0,1]
QuadNode* buildTree(ParticleData *p, int N) {
    QuadNode *root = createNode(0.0, 1.0, 0.0, 1.0);
    for (int i = 0; i < N; i++) {
        insertParticle(root, p, i);
    }
    return root;
}

// Recursively compute mass and center of mass in each node
void computeCenters(QuadNode *node, ParticleData *p) {
    if (!node) return;

    // BAE CASE: If it is a leaf with a single particle
    if (node->particle_index != -1) {
        int idx = node->particle_index;
        node->mass     = p->mass[idx];
        node->center_x = p->pos_x[idx];
        node->center_y = p->pos_y[idx];
        return;
    }

    // Internal node: recursively compute for children and sum up
    node->mass = 0.0;
    node->center_x = 0.0;
    node->center_y = 0.0;
    for (int i = 0; i < 4; i++) {
        if (node->child[i]) {
            computeCenters(node->child[i], p);
            node->mass     += node->child[i]->mass;
            node->center_x += node->child[i]->mass * node->child[i]->center_x;
            node->center_y += node->child[i]->mass * node->child[i]->center_y;
        }
    }
    if (node->mass > 0.0) {
        node->center_x /= node->mass;
        node->center_y /= node->mass;
    }
}

// Free the quadtree recursively
void freeTree(QuadNode *node) {
    if (!node) return;
    for (int i = 0; i < 4; i++) {
        freeTree(node->child[i]);
    }
    free(node);
}

// Barnes-Hut force computation for a single particle
void computeForceFromNode(QuadNode *node, ParticleData *p, int i, double theta, double G,
                          double *fx, double *fy)
{
    if (!node) return;

    // If it is a leaf with a single particle
    if (node->particle_index != -1 && node->particle_index != i) {
        // Compute direct force from that single particle
        double dx = p->pos_x[i] - node->center_x;
        double dy = p->pos_y[i] - node->center_y;
        double dist = sqrt(dx*dx + dy*dy) + EPSILON;
        double F = (p->mass[i] * node->mass) / (dist * dist * dist);
        *fx -= F * dx;
        *fy -= F * dy;
        return;
    }

    // If it is an internal node
    double dx = p->pos_x[i] - node->center_x;
    double dy = p->pos_y[i] - node->center_y;
    double dist = sqrt(dx*dx + dy*dy);
    double width = node->x_max - node->x_min; 
    // Barnes-Hut acceptance criterion
    if ((width / dist) < theta) {
        // Approximate entire node as single mass
        dist += EPSILON;
        double F = (p->mass[i] * node->mass) / (dist * dist * dist);
        *fx -= F * dx;
        *fy -= F * dy;
    } else {
        // Recurse on children
        for (int c = 0; c < 4; c++) {
            computeForceFromNode(node->child[c], p, i, theta, G, fx, fy);
        }
    }
}

// Verlet integration for updating particle positions and velocities
 void updateParticlesVerlet(ParticleData *p, int N, double G, double time_step, QuadNode *root)
{
    double *force_x_old = (double *)calloc(N, sizeof(double));
    double *force_y_old = (double *)calloc(N, sizeof(double));
    double *force_x_new = (double *)calloc(N, sizeof(double));
    double *force_y_new = (double *)calloc(N, sizeof(double));

    if (!force_x_old || !force_y_old || !force_x_new || !force_y_new) {
        fprintf(stderr, "Error: Memory allocation failed in updateParticlesVerlet\n");
        exit(1);
    }

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        double fx = 0.0, fy = 0.0;
        computeForceFromNode(root, p, i, THETA_MAX, G, &fx, &fy);
        force_x_old[i] = fx;
        force_y_old[i] = fy;
    }

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        double ax = (G * force_x_old[i]) / p->mass[i];
        double ay = (G * force_y_old[i]) / p->mass[i];

        p->pos_x[i] += p->vel_x[i] * time_step + 0.5 * ax * time_step * time_step;
        p->pos_y[i] += p->vel_y[i] * time_step + 0.5 * ay * time_step * time_step;
    }

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        double fx = 0.0, fy = 0.0;
        computeForceFromNode(root, p, i, THETA_MAX, G, &fx, &fy);
        force_x_new[i] = fx;
        force_y_new[i] = fy;
    }

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        double ax_old = (G * force_x_old[i]) / p->mass[i];
        double ay_old = (G * force_y_old[i]) / p->mass[i];
        double ax_new = (G * force_x_new[i]) / p->mass[i];
        double ay_new = (G * force_y_new[i]) / p->mass[i];

        p->vel_x[i] += 0.5 * (ax_old + ax_new) * time_step;
        p->vel_y[i] += 0.5 * (ay_old + ay_new) * time_step;
    }

    free(force_x_old);
    free(force_y_old);
    free(force_x_new);
    free(force_y_new);
}

// Symplectic Euler implementation
void updateParticlesSimple(ParticleData *p, int N, double G, double time_step, QuadNode *root)
{
    double *force_x = (double *)calloc(N, sizeof(double));
    double *force_y = (double *)calloc(N, sizeof(double));

    if (!force_x || !force_y) {
        fprintf(stderr, "Error: Memory allocation failed in updateParticlesSimple\n");
        exit(1);
    }

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; i++)
    {
        double fx = 0.0, fy = 0.0;
        computeForceFromNode(root, p, i, THETA_MAX, G, &fx, &fy);
        force_x[i] = fx;
        force_y[i] = fy;
    }

    #pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        p->vel_x[i] += ((G * force_x[i]) / p->mass[i]) * time_step;
        p->vel_y[i] += ((G * force_y[i]) / p->mass[i]) * time_step;

        p->pos_x[i] += p->vel_x[i] * time_step;
        p->pos_y[i] += p->vel_y[i] * time_step;
    }

    free(force_x);
    free(force_y);
}

void updateParticlesBarnesHut(ParticleData *p, int N, double G, double time_step, QuadNode *root)
{
#if VERLET
    updateParticlesVerlet(p, N, G, time_step, root);
#else
    updateParticlesSimple(p, N, G, time_step, root);
#endif
}


int main(int argc, char* argv[])
{
    if (argc != 7 && argc != 8) {
        printf("Usage: %s [N] [input file to read] [nsteps] [delta t] [theta_max] [graphics on/off] [nthreads]", argv[0]);
        return 0;
    }

    int    N_bodies = atoi(argv[1]);
    char  *filename = argv[2];
    int    N_iter   = atoi(argv[3]);
    double delta_t  = atof(argv[4]);
    THETA_MAX       = atof(argv[5]);
    int    graphics = atoi(argv[6]);
    int    nthreads = 1;
    if (argc == 8) {
        nthreads = atoi(argv[7]);
    }

    // Set the number of OpenMP threads if provided
    omp_set_num_threads(nthreads);

    ParticleData particles;
    allocateParticleData(&particles, N_bodies);

    // Read initial data
    if (readParticleData(filename, &particles, N_bodies) != 0) {
        freeParticleData(&particles);
        return 1;
    }

    double G = 100.0 / (double)N_bodies;  // scaled G

    // Initialize graphics if needed
    if (graphics) {
        InitializeGraphics(argv[0]);
    }

    double start_time = get_wall_seconds();

    // Main simulation loop
    for (int step = 0; step < N_iter; step++) {

        // 1) Build quadtree over [0,1]x[0,1]
        QuadNode *root = buildTree(&particles, N_bodies);
        
        // 2) Compute centers of mass
        computeCenters(root, &particles);

        // 3) Update particle velocities and positions using Barnes-Hut
        updateParticlesBarnesHut(&particles, N_bodies, G, delta_t, root);

        // 4) Clean up quadtree
        freeTree(root);

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
    printf("Simulation completed in %f seconds (using %d threads).\n", 
           end_time - start_time, nthreads);

    // Write final positions to  the results file
    writeParticleData("result.gal", &particles, N_bodies);

    if (graphics) {
        CloseGraphics();
    }
    freeParticleData(&particles);

    return 0;
}