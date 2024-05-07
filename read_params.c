#include "read_params.h"

struct global_params get_params(char* filename) {

  struct global_params params;

  FILE *file;

  file = fopen(filename, "r");

  char line[100];
  float value;
  while(fgets(line, sizeof(line), file)) {
    if (strstr(line, "###") != NULL) {
      // Skip this line
      continue;
    }
    char *colon_pos = strchr(line, ':');
    if (colon_pos != NULL) {
      *colon_pos = '\0';
      
      char *key = line; // Key starts from the beginning of the line

      // Remove any trailing whitespace from the key
      while (isspace(*(colon_pos - 1))) {
          *(--colon_pos) = '\0';
      }
      // Skip leading whitespace after the colon
      char *value_str = colon_pos + 1;
      while (isspace(*value_str)) {
          value_str++;
      }
      // Convert value to the appropriate data type (float, int, etc.)
      Value value;

      if ((strcmp(key, "Cosmology") == 0) || (strcmp(key, "Halo_finder") == 0) ||( strcmp(key, "Save_output") == 0)) {
        value.boolean_value = (strstr(colon_pos + 1, "True") != NULL);
      } else {
        sscanf(colon_pos + 1, "%lf", &value.float_value);
      }
        
      // Handle the key and value as needed
      if (strcmp(key, "Cosmology") == 0) {
        params.cosmology = value.boolean_value ? false: true;
      } else if (strcmp(key, "Halo_finder") == 0) {
        params.halo_finder = value.boolean_value ? false: true;
      } else if (strcmp(key, "Save_output") == 0) {
        params.save_output =value.boolean_value ? false: true;
      } else if (strcmp(key, "N_steps") == 0) {
        params.n_steps = (int)value.float_value;
      } else if (strcmp(key, "N_particles") == 0) {
        params.n_prts = (int)value.float_value;
      } else if (strcmp(key, "M_particles") == 0) {
        params.m_prts = value.float_value;
      } else if (strcmp(key, "V_prt_init_max") == 0) {
        params.v_prts_max = value.float_value;
      } else if (strcmp(key, "L_box") == 0) {
        params.l_box = value.float_value;
      } else if (strcmp(key, "Hubble_Param") == 0) {
        params.h = value.float_value;
      } else if (strcmp(key, "Omega_m") == 0) {
        params.Omega_m = value.float_value;
      } else if (strcmp(key, "Omega_Lambda") == 0) {
        params.Omega_Lambda = value.float_value;
      } else if (strcmp(key, "N_bins_hf") == 0) {
        params.N_bins_hf = (int)value.float_value;
      }
    }
  }
  return params;
}

