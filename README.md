# Binary Star Comet Dynamics

This repository contains a MATLAB simulation of a binary star system, where two stars orbit their common barycenter. A planet orbits one of the stars, and a comet enters the system from a distant location. The model uses Newtonian gravitation and the Euler integration method to update the positions and velocities of all bodies in discrete time steps.

## Features
- Full-body gravitational interaction between:
  - Star 1 (M1)
  - Star 2 (M2)
  - Planet (Mp) orbiting Star 2
  - Comet (Mc) approaching from a distant region
- Velocity and position updates using Euler's method
- Real-time collision detection between the comet and other bodies
- Velocity magnitude plotting for all objects
- 2D animated trajectory visualization

## Method
The simulation solves the discrete-time form of Newton's second law for each body:
v[k+1] = v[k] + (F_total[k] / m) * Δt

r[k+1] = r[k] + v[k+1] * Δt

The gravitational forces are computed at each time step and applied accordingly to each mass in the system.

## License
This project is for academic and educational use. Please cite appropriately if referenced.

