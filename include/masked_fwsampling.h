#ifndef MASKED_FWSAMPLING__H
#define MASKED_FWSAMPLING__H

#include "masked_cmp.h"
#include "masked.h"
#include "masked_representations.h"

// BIN2TRI
void masked_binary_to_trinary(uint32_t signs[NSHARES][(N + 31) / 32],
                              const uint32_t poly[NSHARES][(N + 31) / 32]);
void masked_binary_to_trinary_fixed(uint32_t signs[NSHARES][(N + 31) / 32],
                                    const uint32_t poly[NSHARES][(N + 31) / 32],
                                    const uint8_t target_weight);

// I2C
void coeff_to_regular(uint32_t poly[NSHARES][N_PADDED / 32],
                              uint32_t index[NSHARES][W_PADDED]);
void bs_coeff_to_regular(uint32_t poly[NSHARES][N_PADDED / 32],
                              uint32_t index[NSHARES][W_PADDED]);

// CMP
uint32_t sample_cmp(uint32_t out[NSHARES][(N + 31) / 32], uint32_t randomness[NSHARES][N], uint32_t thresh);
uint8_t bs_sample_cmp(uint32_t out[NSHARES][(N + 31) / 32], uint32_t randomness[PAD32(N)/32][NSHARES][L]);

// FISHER-YATES
// expects randomness in p for first W entries
void bs_sample_fisher_yates(uint32_t p[NSHARES][W_PADDED]);
void sample_fisher_yates(uint32_t p[NSHARES][W_PADDED]);

// REJECT
// expects non bitsliced input
// expects rand filled with randomness
int bs_sample_reject(uint32_t p[NSHARES][W_PADDED], uint32_t rand[NSHARES][B_PADDED]);
int sample_reject(uint32_t p[NSHARES][W], uint32_t rand[NSHARES][B_PADDED]);

// REPAND
// expects A to be initialized to zero!
void bs_sample_repand(uint32_t A[NSHARES][(N + 31) / 32]); 

// SORT
// expects randomness in the upper 30bits for the first N entries in x
void bs_sample_sort(uint32_t x[NSHARES][PAD64(N)]);
void sample_sort(uint32_t x[NSHARES][N]);

#endif
