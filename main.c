#include <math.h>
#include <sys/stat.h>
#include "read_params.h"
#include "read_params.c"

// Define Simulation Parameters
// Global Variables used here, 
// but a struct also contains
// each of these fields
int COSMOLOGY;
int SAVE_OUTPUT;
int HALO_FINDER;

int N_PARTICLES;
int N_STEPS;
float M_PARTICLES;
float V_PARTICLES;
int N_BINS;
float L_BOX;

int halo_counts[1] = {0};

float M_MIN = 1e10;
float M_MAX = 1e15;

float G = 6.67430e-11;

float h;
float HUBBLE_CONSTANT;
float OMEGA_M;
float OMEGA_LAMBDA;
//
//
//

// Particle attributes
struct Particle {
    double mass;
    double position[3];
    double velocity[3];
    double acceleration[3];
};

struct Particle particles[1]; // Initialized particles var

void initialize_particles(struct Particle *particles) {
    // Initialize particles with random positions, masses, and velocities
    // according to some global minima and maxima
    for (int i = 0; i < N_PARTICLES; i++) {
        particles[i].mass = rand() / (RAND_MAX + 1.0) * M_PARTICLES;

        particles[i].position[0] = rand() / (RAND_MAX + 1.0) * L_BOX;
        particles[i].position[1] = rand() / (RAND_MAX + 1.0) * L_BOX;
        particles[i].position[2] = rand() / (RAND_MAX + 1.0) * L_BOX;
        particles[i].velocity[0] = rand() / (RAND_MAX + 1.0) * V_PARTICLES - 10.0;
        particles[i].velocity[1] = rand() / (RAND_MAX + 1.0) * V_PARTICLES - 10.0;
        particles[i].velocity[2] = rand() / (RAND_MAX + 1.0) * V_PARTICLES - 10.0;

//        printf("Particle %d Position: (%lf, %lf, %lf)\n", i, particles[i].position[0], particles[i].position[1], particles[i].position[2]);

    }
}

void save_particle_positions(float step, char* filename, struct Particle *particles) {
  printf("\nSaving file at step %f\n", step);
  char* dot_position = strrchr(filename, '.');
  if (dot_position != NULL) {
      *dot_position = '\0'; // Replace the dot with null terminator to remove the extension
  }
    
  // Create directory path
  char directory[100];
  sprintf(directory, "output/%s", filename);
    
  // Create the directory if it doesn't exist
  mkdir(directory, 0777);

  char out_file[100];
  sprintf(out_file, "%s/particle_positions_%f.csv", directory, step);
  FILE* fp = fopen(out_file, "w");
  if (fp == NULL) {
      printf("Error saving file.\n");
      return;
  }

  fprintf(fp, "X,Y,Z\n");
  for (int i = 0; i < N_PARTICLES; i++) {
      fprintf(fp, "%f,%f,%f\n", particles[i].position[0], particles[i].position[1], particles[i].position[2]);
  }

  fclose(fp);
}

double calculate_hubble_parameter(double a) {
  if ((int)COSMOLOGY == 0) {
    // Calculate the Hubble parameter 
    // H(a) = H0 * sqrt(OMEGA_M / a^3 + OMEGA_LAMBDA)
    return HUBBLE_CONSTANT * sqrt(OMEGA_M / pow(a, 3) + OMEGA_LAMBDA);
  } else {
    return 1.0;  
  }
}

void calculate_friedmann_equations(double a, double* hubble_parameter, double* acceleration_parameter) {
  if ((int)COSMOLOGY == 0) {
    // Calculate the Hubble parameter and acceleration 
    // parameter using Friedmann equations
    *hubble_parameter = calculate_hubble_parameter(a);
    *acceleration_parameter = -4 * M_PI * G * (OMEGA_M / (a * a * a) + 3 * OMEGA_LAMBDA) / 3;
  }
}

void calculate_forces(double a) {
    // Calculate gravitational forces between particles
    // and account for cosmological expansion
//    double hubble_parameter, acceleration_parameter;
//    calculate_friedmann_equations(
//	a, 
//	&hubble_parameter, 
//	&acceleration_parameter
//   );

    int acceleration_parameter = 1;

    for (int i = 0; i < N_PARTICLES; i++) {
        double ax = 0.0, ay = 0.0, az = 0.0;  // Acceleration components for particle i

        for (int j = 0; j < N_PARTICLES; j++) {
            if (i != j) {
                double dx = particles[j].position[0] - particles[i].position[0];
                double dy = particles[j].position[1] - particles[i].position[1];
                double dz = particles[j].position[2] - particles[i].position[2];
                double r = sqrt(dx*dx + dy*dy + dz*dz);  // Distance between particles

                // Gravity
                double f = G * particles[i].mass * particles[j].mass / pow(r, 2);

                // Update acceleration components for particle i
                ax += f * dx / r;
                ay += f * dy / r;
                az += f * dz / r;
            }
        }

        if ((int)COSMOLOGY == 0) {
          // Account for cosmic acceleration
          ax += acceleration_parameter * particles[i].position[0];
          ay += acceleration_parameter * particles[i].position[1];
          az += acceleration_parameter * particles[i].position[2];
        }

        // Update velocities of particle i
        particles[i].velocity[0] += ax;
        particles[i].velocity[1] += ay;
        particles[i].velocity[2] += az;
    }
}

void update_positions(double a) {
    // Update positions of particles and account for cosmological expansion
    for (int i = 0; i < N_PARTICLES; i++) {
        particles[i].position[0] += particles[i].velocity[0];
        particles[i].position[1] += particles[i].velocity[1];
        particles[i].position[2] += particles[i].velocity[2];

        printf("Particle %d Position: (%lf, %lf, %lf)\n", i, particles[i].position[0], particles[i].position[1], particles[i].position[2]);
    }
}

void compute_halo_mass_function(int N_BINS) {
  // Compute the halo mass function from the simulation data
  int halo_counts[N_BINS];
  double bin_width = log10(M_MAX / M_MIN) / N_BINS;

  // Loop through particles and identify halos
  for (int i = 0; i < N_PARTICLES; i++) {
      double mass = particles[i].mass;
      if (mass >= M_MIN && mass <= M_MAX) {
          int bin_index = (int)((log10(mass) - log10(M_MIN)) / bin_width);
          if (bin_index >= 0 && bin_index < N_BINS) {
              halo_counts[bin_index]++;
          }
      }
  }

  // Normalize halo counts by the volume of the simulation box
  double box_volume = 100.0 * 100.0 * 100.0;  // Example: 100 Mpc/h box size
  for (int i = 0; i < N_BINS; i++) {
      double mass = M_MIN * pow(10, (i + 0.5) * bin_width);
      double halo_density = halo_counts[i] / box_volume / bin_width;
      printf("Halo Mass: %.2e - %.2e, Halo Count: %d, Halo Density: %.2e\n",
             mass, mass * pow(10, bin_width), halo_counts[i], halo_density);
  }
}

void calculate_acceleration(struct Particle *particles, int num_particles) {
    for (int i = 0; i < num_particles; i++) {
        for (int j = 0; j < num_particles; j++) {
            if (i != j) {
                double dx = particles[j].position[0] - particles[i].position[0];
                double dy = particles[j].position[1] - particles[i].position[1];
                double dz = particles[j].position[2] - particles[i].position[2];
                double r = sqrt(dx*dx + dy*dy + dz*dz);
                double F = G * particles[i].mass * particles[j].mass / (r * r);
                particles[i].acceleration[0] = F * dx / r / particles[i].mass;
                particles[i].acceleration[1] = F * dy / r / particles[i].mass;
                particles[i].acceleration[2] = F * dz / r / particles[i].mass;
            }
        }
    }
}

void update_positions_velocities(struct Particle *particles, int num_particles, double dt) {
    for (int i = 0; i < num_particles; i++) {
        particles[i].velocity[0] += particles[i].acceleration[0] * dt;
        particles[i].velocity[1] += particles[i].acceleration[1] * dt;
        particles[i].velocity[2] += particles[i].acceleration[2] * dt;
        particles[i].position[0] += particles[i].velocity[0] * dt;
        particles[i].position[1] += particles[i].velocity[1] * dt;
        particles[i].position[2] += particles[i].velocity[2] * dt;
    }
    for (int i = 0; i < num_particles; i++) {
      printf("%f, %f, %f\n", particles[i].position[0], particles[i].position[1], particles[i].position[2]);
    }
}

void print_global_params(struct global_params params) {
  // Print each member of the struct
  // Used for debug of reading in file parameters
  printf("cosmology: %d\n", params.cosmology);
  printf("halo_finder: %d\n", params.halo_finder);
  printf("save_output: %d\n", params.save_output);
  printf("n_steps: %d\n", params.n_steps);
  printf("n_prts: %d\n", params.n_prts);
  printf("m_prts: %f\n", params.m_prts);
  printf("v_prts_max: %f\n", params.v_prts_max);
  printf("h: %f\n", params.h);
  printf("Omega_m: %f\n", params.Omega_m);
  printf("Omega_Lambda: %f\n", params.Omega_Lambda);
  printf("N_bins_hf: %d\n", params.N_bins_hf);
}

int main(int argc, char *argv[]) {
  char* filename = get_filename(argc, argv);
  print_params(filename);
  struct global_params simulation_params = get_params(filename);

  COSMOLOGY = simulation_params.cosmology;
  HALO_FINDER = simulation_params.halo_finder;
  SAVE_OUTPUT = simulation_params.save_output;

  N_STEPS = simulation_params.n_steps;
  N_PARTICLES = simulation_params.n_prts;
  M_PARTICLES = simulation_params.m_prts;
  V_PARTICLES = simulation_params.v_prts_max;
  L_BOX = simulation_params.l_box;

  h = simulation_params.h;
  HUBBLE_CONSTANT = h * 100;
  OMEGA_M = simulation_params.Omega_m;
  OMEGA_LAMBDA = simulation_params.Omega_Lambda;
  N_BINS = simulation_params.N_bins_hf;

  struct Particle particles[N_PARTICLES];

  // Initialize particles
  initialize_particles(particles);

  double dt = 0.01;
  double t_end = dt * N_STEPS;
  printf("%f", t_end);

  for (double t = 0.0; t < t_end; t += dt) {
        // Calculate accelerations
        calculate_acceleration(particles, N_PARTICLES);

        // Update positions and velocities
        update_positions_velocities(particles, N_PARTICLES, dt);

        // Output positions
        save_particle_positions(t, filename, particles);
    
  }

 
  return 0;

    double a = 0.0;  // Initial scale factor
    
    // Main simulation loop
    for (int step = 0; step < N_STEPS; step++) {
        // Calculate gravitational forces with cosmological expansion
        calculate_forces(a);

        // Update particle positions with cosmological expansion
        update_positions(a);

      if ((int)COSMOLOGY == 0) {
        // Update scale factor for next time step (using Euler integration for simplicity)
        double hubble_parameter, acceleration_parameter;
        calculate_friedmann_equations(a, &hubble_parameter, &acceleration_parameter);
        a += 0.01 * hubble_parameter;  // Increment scale factor
      }

      // Save particle positions to CSV file
//      if ((int)SAVE_OUTPUT == 0) {
//       save_particle_positions(step, filename);
//     }

      // Compute and analyze halo mass function
      if (step % 100 == 0) {
        if ((int)HALO_FINDER == 0) {
          printf("Finding Haloes at step %i\n\n", step);
          compute_halo_mass_function(N_BINS);
        }
      }
    }
    return 0;
}
