#include "utility.h"

void init_difference_table(difference_table_t diff_table) {
  byte_t col = 0, row = 0;
  for(byte_t i = 0 ; i < SBOX_SIZE ; i++){
    for (byte_t j = 0 ; j < SBOX_SIZE ; j++) {
      row = i^j;
      col = sbox[i]^sbox[j];
      diff_table[row * SBOX_SIZE + col]++;
    }
  }
}

void print_difference_table(difference_table_t diff_table) {
  for (int i = 0; i < SBOX_SIZE; i++) {
    for (int j = 0; j < SBOX_SIZE; j++) {
      printf("%02hhu ", diff_table[i*SBOX_SIZE+j]);
    }
    printf("\n\033[0m");
  }
}

byte_t previous_highest_value_difference_table(byte_t diff_table[16][16]) {
  byte_t res = 0;
  byte_t res_prev = 0;
  for (unsigned int i = 0; i < 16; i++) {
    for (unsigned int j = 0; j < 16; j++) {
      if (diff_table[i][j] != 16 && diff_table[i][j] > res) {
        res_prev = res;
        res = diff_table[i][j];
      }
    }
  }
  return res_prev;
}

#define NEQ16(X,Y) ((X) != 16) ? (X) : (Y)
#define MAX(X,Y) ((X) > (Y)) ? (X) : (Y)

byte_t previous_highest_value_difference_table_2(difference_table_t diff_table) {
  byte_t res = 0;
  byte_t res_prev = 0;
  for (unsigned int i = 0; i < SBOX_SIZE*SBOX_SIZE; i++) {
    if (diff_table[i] != 16 && diff_table[i] > res) {
      res_prev = res;
      res = diff_table[i];
    }
  }
  return res_prev;
}

void get_stats(difference_table_t diff_table) {

  typedef struct counter {
    int16_t value;
    uint16_t count;
  } counter_t;

  counter_t counter_array[SBOX_SIZE*SBOX_SIZE];
  for(int i = 0; i < SBOX_SIZE*SBOX_SIZE; i++) {
    counter_array[i].value = -1;
    counter_array[i].count = 0;
  }
  
  for(int i = 0; i < SBOX_SIZE*SBOX_SIZE; i++) {
    counter_array[i].value = -1;
    counter_array[i].count = 0;
  }
  
  void increment_counter(counter_t *counter_array, int value) {
    int pos = 0;
    int v = counter_array[pos].value;
    while (pos < SBOX_SIZE*SBOX_SIZE && v!=value && v!=-1 ) {
      pos++;
      v = counter_array[pos].value;
    }
    if (counter_array[pos].value == value) {
      counter_array[pos].count++;
    } else if (counter_array[pos].value == -1) {
      counter_array[pos].value = value;
      counter_array[pos].count = 1;
    }
  }

  byte_t value = 0;
  int pos = 0;
  
  for(int i = 0; i < SBOX_SIZE*SBOX_SIZE; i++) {
      value = diff_table[i];
      increment_counter(counter_array, value);    
  }

  int size_tab_diff = SBOX_SIZE*SBOX_SIZE;
  int count = 0;
  while (counter_array[pos].value != -1) {
    count = counter_array[pos].count;
    value = counter_array[pos].value;
    float freq = count / ((float)size_tab_diff);
    fprintf(stderr, "Frequency %02d = %02.3f%% (%d/%d)\n", value, 100*freq, count, size_tab_diff);
    pos++;
  }
}
