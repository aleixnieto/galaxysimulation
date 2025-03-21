# High-Performance Programming in C
## Individual Final Project: Barnes-Hut Method

This repository contains my individual final project for the course "High-Performance Programming" at Uppsala University, Spring 2023. The project involves implementing and optimizing the Barnes-Hut algorithm, a divide-and-conquer approach for efficiently simulating N-body gravitational interactions.

### Project Overview

The goal of this project is to simulate the gravitational forces between particles (or stars) using Newton's law of gravitation and Plummer spheres, where the Barnes-Hut algorithm significantly reduces the computational complexity from O(N²) to approximately O(N log N). The implementation employs the Velocity Verlet method for time integration, ensuring second-order accuracy.

### Key Components:
- **Barnes-Hut Algorithm:**
  - Quad-tree data structure to efficiently group distant particles.
  - Recursive computation of mass and center of mass for particle groups.
  - Approximation threshold parameter (θ_max) tuning for accuracy and performance trade-off.

- **Velocity Verlet Integration:**
  - Accurate numerical method for particle position and velocity updates.
  - Stability analysis and parameter tuning for time steps (Δt).

- **Performance Optimization:**
  - Parallelization using OpenMP.
  - Runtime and complexity analysis comparing Barnes-Hut against direct computation.
  - Code optimization techniques learned throughout the course to achieve maximal performance.

This project showcases efficient implementation practices in C, with a strong emphasis on numerical accuracy, computational efficiency, and parallel computing.

