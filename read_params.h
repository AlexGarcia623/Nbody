#ifndef READ_PARAMS_H
#define READ_PARAMS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <ctype.h>

typedef union {
  bool boolean_value;
  double float_value;
} Value;

struct global_params {
  bool merger;
  bool cosmology;
  bool halo_finder;
  bool save_output;
  bool edge;
  int snapshot_cadence;
  float delta_t;
  float fb_strength;
  int n_steps;
  int n_prts;
  float m_prts;
  float v_prts_max;
  float l_box;
  float grav_softening;
  float h;
  float Omega_m;
  float Omega_Lambda;
  int N_bins_hf;
};

int max_label_width(FILE *file) {
    char line[100];
    int max_width = 0;

    // Loop through each line in the file
    while (fgets(line, sizeof(line), file)) {
        char *colon_pos = strchr(line, ':');
        if (colon_pos != NULL) {
            *colon_pos = '\0'; // Replace colon with null terminator to split the string
            
            // Update max_width if needed
            int width = strlen(line);
            if (width > max_width) {
                max_width = width;
            }
        }
    }

    return max_width + 1;
}

void print_title() {
    printf("\n");
    printf("  /$$$$$$                                              \n");
    printf(" /$$__  $$                                             \n");
    printf("| $$  \\__/  /$$$$$$   /$$$$$$$ /$$$$$$/$$$$   /$$$$$$ \n");
    printf("| $$       /$$__  $$ /$$_____/| $$_  $$_  $$ /$$__  $$\n");
    printf("| $$      | $$  \\ $$|  $$$$$$ | $$ \\ $$ \\ $$| $$  \\ $$\n");
    printf("| $$    $$| $$  | $$ \\____  $$| $$ | $$ | $$| $$  | $$\n");
    printf("|  $$$$$$/|  $$$$$$/ /$$$$$$$/| $$ | $$ | $$|  $$$$$$/\n");
    printf(" \\______/  \\______/ |_______/ |__/ |__/ |__/ \\______/ \n");
    printf("                                                      \n");
    printf("                                                      \n");
    printf("                                                      \n");
    printf(" /$$   /$$ /$$                       /$$              \n");
    printf("| $$$ | $$| $$                      | $$              \n");
    printf("| $$$$| $$| $$$$$$$   /$$$$$$   /$$$$$$$ /$$   /$$    \n");
    printf("| $$ $$ $$| $$__  $$ /$$__  $$ /$$__  $$| $$  | $$    \n");
    printf("| $$  $$$$| $$  \\ $$| $$  \\ $$| $$  | $$| $$  | $$    \n");
    printf("| $$\\  $$$| $$  | $$| $$  | $$| $$  | $$| $$  | $$    \n");
    printf("| $$ \\  $$| $$$$$$$/|  $$$$$$/|  $$$$$$$|  $$$$$$$    \n");
    printf("|__/  \\__/|_______/  \\______/  \\_______/ \\____  $$    \n");
    printf("                                         /$$  | $$    \n");
    printf("                                        |  $$$$$$/    \n");
    printf("                                         \\______/     \n");
    printf("\n\n");

    sleep(1);

    printf("!!!Welcome to Cosmo Nbody!!!\n");
    printf("You're running an N-body Simulation with the following parameters:\n\n");
}


char* get_filename(int argc, char *argv[]) {
  print_title();

  char* fname = malloc(100* sizeof(char));

  struct global_params params = {false, false, false};

  FILE *file;

  char default_file[] = "default_simulation_params.txt";
  char *filename;

  if (argc == 1) {
    printf("No file provided, using %s\n\n", default_file);
    filename = default_file;
  } else if (argc == 2) {
    filename = argv[1];
  } else {
    printf("!!!This only takes one argument!!!");
    exit(1);
  }
  strcpy(fname, filename);

  fname[strlen(filename)] = '\0';

  return fname;
}

void print_params(char* filename) {
  FILE *file;
  file = fopen(filename,"r");

  char line[100]; // Max line length of 100

  if (file == NULL) {
    printf("NOT ABLE TO OPEN FILE");
    exit(1);
  }

  int label_width = max_label_width(file);

  fseek(file, 0, SEEK_SET);

  while(fgets(line, sizeof(line), file)) {
    char *colon_pos = strchr(line, ':');
    if (colon_pos != NULL) {
      *colon_pos = '\0';
      
      Value value;

      if (strstr(colon_pos + 1, "True") || strstr(colon_pos + 1, "False")) {
        value.boolean_value = (strstr(colon_pos + 1, "True") != NULL);
        printf("\t%-*s: %s\n", label_width, line, value.boolean_value ? "True": "False");
      } else {
        sscanf(colon_pos + 1, "%lf", &value.float_value);
        printf("\t%-*s: %.6lf\n", label_width, line, value.float_value);
      }
    }
  }

  printf("\n");

  fclose(file);
  sleep(2);
}

struct global_params get_params(char* filename);

#endif
