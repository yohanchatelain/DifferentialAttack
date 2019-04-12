#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "inputs.h"

inputs_t *inputs = NULL;

void parse_inputs(char *argv[]) {
  
  inputs->random_couples = atoi(argv[1]);
  strcpy(inputs->cryptosystem,argv[2]);  
  inputs->print_key_step_by_step = argv[3][0];
  
  printf("Parsed inputs: %d\n", inputs->random_couples);
  printf("Parsed inputs: %s\n", inputs->cryptosystem);
  printf("Parsed inputs: %c\n", inputs->print_key_step_by_step);
}
