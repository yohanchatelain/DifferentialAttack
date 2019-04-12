#include "caracteristic.h"


static void
permute(diff_carac_tab_t res) {
  unsigned int cpt = NB_CARACTERISTICS;
  differential_caracterictic_t tmp;
  while (cpt != 0 && res[cpt].proba > res[cpt-1].proba) {
    tmp = res[cpt];
    res[cpt] = res[cpt-1];
    res[cpt-1] = tmp;
    cpt--;
  }
}

void
differential_caracteristic(diff_carac_tab_t res, fct_perm f, byte_t difference_table[16][16]) {

  const byte_t max = previous_highest_value_difference_table(difference_table);
  const unsigned int proba_min_round = max*max*max*max;

  unsigned long long int proba;

  unsigned long long int proba_i, proba_j, proba_k;
  for (block_t i = 1; i != 0; i++) {
    for (block_t i2 = 1; i2 != 0; i2++) {
      proba_i = difference_table[(i & 0xf)][(i2 & 0xf)];
      proba_i *= difference_table[((i & 0xf0) >> 4)][((i2 & 0xf0) >> 4)];
      proba_i *= difference_table[((i & 0xf00) >> 8)][((i2 & 0xf00) >> 8)];
      proba_i *= difference_table[((i & 0xf000) >> 12)][((i2 & 0xf000) >> 12)];
      if (proba_i == 0 || proba_i < proba_min_round) {
	continue;
      }
      block_t j = f(i2);
      for (block_t j2 = 1; j2 != 0; j2++) {
	proba_j = difference_table[(j & 0xf)][(j2 & 0xf)];
	proba_j *= difference_table[((j & 0xf0) >> 4)][((j2 & 0xf0) >> 4)];
	proba_j *= difference_table[((j & 0xf00) >> 8)][((j2 & 0xf00) >> 8)];
	proba_j *= difference_table[((j & 0xf000) >> 12)][((j2 & 0xf000) >> 12)];
	if (proba_j == 0 || proba_j < proba_min_round) {
	  continue;
	}
	block_t k = f(j2);
	for (block_t k2 = 1; k2 != 0; k2++) {
	  proba_k = difference_table[(k & 0xf)][(k2 & 0xf)];
	  proba_k *= difference_table[((k & 0xf0) >> 4)][((k2 & 0xf0) >> 4)];
	  proba_k *= difference_table[((k & 0xf00) >> 8)][((k2 & 0xf00) >> 8)];
	  proba_k *= difference_table[((k & 0xf000) >> 12)][((k2 & 0xf000) >> 12)];
	  if (proba_k == 0 || proba_k < proba_min_round) {
	    continue;
	  }
	  proba = proba_i * proba_j * proba_k;
	  if (proba > res[NB_CARACTERISTICS].proba) {
	    res[NB_CARACTERISTICS].input_diff = i;
	    res[NB_CARACTERISTICS].output_diff = f(k2);
	    res[NB_CARACTERISTICS].proba = proba;
#ifdef PATH_INFO
	    res[NB_CARACTERISTICS].output_diff_round_1 = i2;
	    res[NB_CARACTERISTICS].input_diff_round_2 = j;
	    res[NB_CARACTERISTICS].output_diff_round_2 = j2;
	    res[NB_CARACTERISTICS].input_diff_round_3 = k;
	    res[NB_CARACTERISTICS].output_diff_round_3 = k2;
#endif
	    permute(res);       
	  }
	}
      }
    }
  }
}

uint64_t medium_table[SBOX_SIZE*SBOX_SIZE*SBOX_SIZE*SBOX_SIZE] = {0};

void compute_medium_diff(difference_table_t difference_table, uint64_t *medium_table) {
  uint64_t proba;  
  uint32_t index[4];
  uint32_t index_i, index_j;
  
  for (block_t i = 0; i < 1<<8; i++) {
    for (block_t j = 0; j < 1<<8; j++) {
      index[0] = (i & 0x000f) * SBOX_SIZE + (j & 0x000f);
      index[1] = (i & 0x00f0) +  ((j & 0x00f0) >> 4);
     
      proba  = difference_table[index[0]];
      proba *= difference_table[index[1]];

      index_i = i & 0x00ff;
      index_j = j & 0x00ff;

      medium_table[index_i * SBOX_SIZE*SBOX_SIZE + index_j] = proba;
    }
  }
}
#ifndef DEBUG
inline
#endif
uint64_t compute_proba(block_t i, block_t j) {
  return medium_table[(i & 0x00ff) * SBOX_SIZE*SBOX_SIZE + (j&0x00ff)]*
    medium_table[(i & 0xff00) + ((j&0xff00) >> 8)];
}

void
differential_caracteristic_2(diff_carac_tab_t res, difference_table_t difference_table) {
  
  const byte_t max = previous_highest_value_difference_table_2(difference_table);
  const uint64_t proba_min_round = max*max*max*max;
  const uint32_t twoTo16Minus1 = (1 << 16) - 1;
  
  unsigned long long int proba;

  uint64_t proba_i, proba_j, proba_k;

  compute_medium_diff(difference_table, medium_table);
  
#pragma omp parallel for shared(proba)
  for (block_t i = 1; i < twoTo16Minus1; i++) {      	  
    for (block_t i2 = 1; i2 < twoTo16Minus1; i2++) {
      proba_i = compute_proba(i, i2);
      if (proba_i == 0 || proba_i < proba_min_round) {
	continue;
      }
      
      /* block_t j = f(i2); */
      block_t j = get_heys_perm(i2);
      for (block_t j2 = 1; j2 < twoTo16Minus1; j2++) {
	proba_j = compute_proba(j, j2);

	if (proba_j == 0 || proba_j < proba_min_round) {
	  continue;
	}
	
	/* block_t k = f(j2); */
	block_t k = get_heys_perm(j2);	
	for (block_t k2 = 1; k2 < twoTo16Minus1; k2++) {
	  proba_k = compute_proba(k, k2);	  
	  
	  if (proba_k == 0 || proba_k < proba_min_round) {
	    continue;
	  }
#pragma omp critical
	  {	    
	    proba = proba_i * proba_j * proba_k;
	    if (proba > res[NB_CARACTERISTICS].proba) {
	      res[NB_CARACTERISTICS].input_diff = i;
	      res[NB_CARACTERISTICS].output_diff = get_heys_perm(k2);
	      res[NB_CARACTERISTICS].proba = proba;
#ifdef PATH_INFO
	      res[NB_CARACTERISTICS].output_diff_round_1 = i2;
	      res[NB_CARACTERISTICS].input_diff_round_2 = j;
	      res[NB_CARACTERISTICS].output_diff_round_2 = j2;
	      res[NB_CARACTERISTICS].input_diff_round_3 = k;
	      res[NB_CARACTERISTICS].output_diff_round_3 = k2;
#endif
	      permute(res);
	    }
	  }
	}
      }
    }
  }
}
