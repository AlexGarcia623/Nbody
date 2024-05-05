//#include <stdio.h>

//#include <stdlib.h>
#include <math.h>
#include "read_params.h"
#include "read_params.c"

//#define SAVE_OUTPUT
//#define COSMOLOGY
//#define HALO_FINDER

int COSMOLOGY;
int SAVE_OUTPUT;
int HALO_FINDER;

int N_PARTICLES;
int N_STEPS;
float M_PARTICLES;
float V_PARTICLES;
int N_BINS;

float M_MIN = 1e10;
float M_MAX = 1e15;

float G = 6.67430e-11;

float h;
float HUBBLE_CONSTANT;
float OMEGA_M;
float OMEGA_LAMBDA;


struct Particle {
    double mass;
    double position[3];
    double velocity[3];
};

void initialize_particles(struct Particle *particles) {
    // Initialize particles with random positions and velocities
    for (int i = 0; i < N_PARTICLES; i++) {
        particles[i].mass = rand() / (RAND_MAX + 1.0) * M_PARTICLES; 
        // (RAND_MAX + 1.0) * 1e15;

        particles[i].position[0] = rand() / (RAND_MAX + 1.0) * 100.0;
        particles[i].position[1] = rand() / (RAND_MAX + 1.0) * 100.0;
        particles[i].position[2] = rand() / (RAND_MAX + 1.0) * 100.0;
        particles[i].velocity[0] = rand() / (RAND_MAX + 1.0) * 20.0 - 10.0;
        particles[i].velocity[1] = rand() / (RAND_MAX + 1.0) * 20.0 - 10.0;
        particles[i].velocity[2] = rand() / (RAND_MAX + 1.0) * 20.0 - 10.0;

        printf("Particle %d Position: (%lf, %lf, %lf)\n", i, particles[i].position[0], particles[i].position[1], particles[i].position[2]);

    }
}

void save_particle_positions(int step) {
    char filename[50];
    sprintf(filename, "particle_positions_%d.csv", step);
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error opening file.\n");
        return;
    }

    fprintf(fp, "X,Y,Z\n");
    for (int i = 0; i < N_PARTICLES; i++) {
        fprintf(fp, "%lf,%lf,%lf\n", particles[i].position[0], particles[i].position[1], particles[i].position[2]);
    }

    fclose(fp);
}

double calculate_hubble_parameter(double a) {
#ifdef COSMOLOGY
    // Calculate the Hubble parameter 
    // H(a) = H0 * sqrt(OMEGA_M / a^3 + OMEGA_LAMBDA)
    return HUBBLE_CONSTANT * sqrt(OMEGA_M / pow(a, 3) + OMEGA_LAMBDA);
#else
    return 1.0;  
#endif
}

void calculate_friedmann_equations(double a, double* hubble_parameter, double* acceleration_parameter) {
#ifdef COSMOLOGY
    // Calculate the Hubble parameter and acceleration 
    // parameter using Friedmann equations
    *hubble_parameter = calculate_hubble_parameter(a);
    *acceleration_parameter = -4 * M_PI * G * (OMEGA_M / (a * a * a) + 3 * OMEGA_LAMBDA) / 3;
#endif
}

void calculate_forces(double a) {
    // Calculate gravitational forces between particles
    // and account for cosmological expansion
    double hubble_parameter, acceleration_parameter;
    calculate_friedmann_equations(
	a, 
	&hubble_parameter, 
	&acceleration_parameter
    );

    for (int i = 0; i < N_PARTICLES; i++) {
        double ax = 0.0, ay = 0.0, az = 0.0;  // Acceleration components for particle i

        for (int j = 0; j < N_PARTICLES; j++) {
            if (i != j) {
                double dx = particles[j].position[0] - particles[i].position[0];
                double dy = particles[j].position[1] - particles[i].position[1];
                double dz = particles[j].position[2] - particles[i].position[2];
                double r = sqrt(dx*dx + dy*dy + dz*dz);  // Distance between particles

                // Gravitational force calculation (Newton's law of gravitation) with cosmological expansion
                double f = G * particles[i].mass * particles[j].mass / pow((r * a), 2);

                // Update acceleration components for particle i
                ax += f * dx / r;
                ay += f * dy / r;
                az += f * dz / r;
            }
        }

#ifdef COSMOLOGY
        // Account for cosmic acceleration
        ax += acceleration_parameter * particles[i].position[0];
        ay += acceleration_parameter * particles[i].position[1];
        az += acceleration_parameter * particles[i].position[2];
#endif

        // Update velocities of particle i
        particles[i].velocity[0] += ax;
        particles[i].velocity[1] += ay;
        particles[i].velocity[2] += az;
    }
}

void update_positions(double a) {
    // Update positions of particles and account for cosmological expansion
    for (int i = 0; i < N_PARTICLES; i++) {
        particles[i].position[0] += particles[i].velocity[0] / a;
        particles[i].position[1] += particles[i].velocity[1] / a;
        particles[i].position[2] += particles[i].velocity[2] / a;

        printf("Particle %d Position: (%lf, %lf, %lf)\n", i, particles[i].position[0], particles[i].position[1], particles[i].position[2]);
    }
}

void compute_halo_mass_function() {
#ifdef HALO_FINDER
    // Compute the halo mass function from the simulation data
    int halo_counts[N_BINS] = {0};
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
#endif
}

void print_global_params(struct global_params params) {
    // Print each member of the struct
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
    char* filename = print_params(argc, argv);
    struct global_params simulation_params = get_params(filename);

    COSMOLOGY = simulation_params.cosmology;
    HALO_FINDER = simulation_params.halo_finder;
    SAVE_OUTPUT = simulation_params.save_output;

    N_STEPS = simulation_params.n_steps;
    N_PARTICLES = simulation_params.n_prts;
    M_PARTICLES = simulation_params.m_prts;
    V_PARTICLES = simulation_params.v_prts_max;

    h = simulation_params.h;
    HUBBLE_CONSTANT = h * 100;
    OMEGA_M = simulation_params.Omega_m;
    OMEGA_LAMBDA = simulation_params.Omega_Lambda;
    N_BINS = simulation_params.N_bins_hf;

    Particle particles[N_PARTICLES];

    if (COSMOLOGY == 0) {
      printf("Cosmology is enabled");
    } else {
      printf("No Cosmology");
    }
    
    //print_global_params(simulation_params);   
 
    return 0;

    double a = 1.0;  // Initial scale factor (normalized to 1)
    
    // Initialize particles
    initialize_particles(particle);

    // Main simulation loop
    for (int step = 0; step < N_STEPS; step++) {
        // Calculate gravitational forces with cosmological expansion
        calculate_forces(a);

        // Update particle positions with cosmological expansion
        update_positions(a);

#ifdef COSMOLOGY
        // Update scale factor for next time step (using Euler integration for simplicity)
        double hubble_parameter, acceleration_parameter;
        calculate_friedmann_equations(a, &hubble_parameter, &acceleration_parameter);
        a += 0.01 * hubble_parameter;  // Increment scale factor
#endif

        // Save particle positions to CSV file
#ifdef SAVE_OUTPUT
        save_particle_positions(step);
#endif        

        // Compute and analyze halo mass function
        if (step % 100 == 0) {
          compute_halo_mass_function();

          printf("Step %d: Particle positions updated.\n", step);
        }
    }
    return 0;
}
