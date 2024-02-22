#include "gadgets.h"
#include "masked.h"
#include "masked_cmp.h"
#include "masked_fwsampling.h"
#include "core_ops.h"
#include <stddef.h>
#include <stdint.h>

static inline void sample_fixed_32(uint32_t *out, size_t stride,
                            uint32_t randomness[NSHARES][L]) {
#ifdef NTRULPR
#if N == 653 || N == 857 || N == 953 || N == 1013 // t = 3
    b_and(out, stride, &randomness[0][0], L, &randomness[0][1],
               L);
    b_not(out, stride, out, stride);
    for (size_t i = 2; i < L; i++) {
        b_and(out, stride, &randomness[0][i], L, out, stride);
    }
#elif N == 761  // t = 21
    b_and(out, stride, &randomness[0][0], L, &randomness[0][1],
               L);
    b_not(out, stride, out, stride);
    b_and(out, stride, &randomness[0][2], L, out, stride);
    b_not(out, stride, out, stride);
    b_and(out, stride, &randomness[0][3], L, out, stride);
    b_not(out, stride, out, stride);
    b_and(out, stride, &randomness[0][4], L, out, stride);
    b_not(out, stride, out, stride);
    b_and(out, stride, &randomness[0][5], L, out, stride);
#elif N == 1277 // t = 11
    b_and(out, stride, &randomness[0][0], L, &randomness[0][1],
               L);
    b_not(out, stride, out, stride);
    b_and(out, stride, &randomness[0][2], L, out, stride);
    b_not(out, stride, out, stride);
    b_and(out, stride, &randomness[0][3], L, out, stride);
    b_not(out, stride, out, stride);
    for (size_t i = 4; i < L; i++) {
        b_and(out, stride, &randomness[0][i], L, out, stride);
    }
#else
#error
#endif
#else
#if N == 12323 || N == 35851 || N == 677 || N == 761 || N == 857 ||            \
    N == 1277                               // t = 3
    b_and(out, stride, &randomness[0][0], L, &randomness[0][1],
               L);
    b_not(out, stride, out, stride);
    for (size_t i = 2; i < L; i++) {
        b_and(out, stride, &randomness[0][i], L, out, stride);
    }
#elif N == 24659 || N == 17669 || N == 8192 // t = 1
    b_and(out, stride, &randomness[0][0], L, &randomness[0][1],
               L);
    for (size_t i = 2; i < L; i++) {
        b_and(out, stride, &randomness[0][i], L, out, stride);
    }
#elif N == 40973 || N == 1013 || N == 653   // t = 7
    b_and(out, stride, &randomness[0][0], L, &randomness[0][1],
               L);
    b_and(out, stride, &randomness[0][2], L, out, stride);
    b_not(out, stride, out, stride);
    for (size_t i = 3; i < L; i++) {
        b_and(out, stride, &randomness[0][i], L, out, stride);
    }
#elif N == 57637 || N == 6960               // t = 9
    b_and(out, stride, &randomness[0][0], L, &randomness[0][1],
               L);
    b_and(out, stride, &randomness[0][2], L, out, stride);
    b_not(out, stride, out, stride);
    b_and(out, stride, &randomness[0][3], L, out, stride);
    b_not(out, stride, out, stride);
    for (size_t i = 4; i < L; i++) {
        b_and(out, stride, &randomness[0][i], L, out, stride);
    }
#elif N == 3488 || N == 6688 || N == 821    // t = 5
    b_and(out, stride, &randomness[0][0], L, &randomness[0][1],
               L);
    b_not(out, stride, out, stride);
    b_and(out, stride, &randomness[0][2], L, out, stride);
    b_not(out, stride, out, stride);
    for (size_t i = 3; i < L; i++) {
        b_and(out, stride, &randomness[0][i], L, out, stride);
    }
#elif N == 4608                             // t = 11
    b_and(out, stride, &randomness[0][0], L, &randomness[0][1],
               L);
    b_not(out, stride, out, stride);
    b_and(out, stride, &randomness[0][2], L, out, stride);
    b_not(out, stride, out, stride);
    b_and(out, stride, &randomness[0][3], L, out, stride);
    b_not(out, stride, out, stride);
    for (size_t i = 4; i < L; i++) {
        b_and(out, stride, &randomness[0][i], L, out, stride);
    }
#elif N == 509                              // t = 1
    for (size_t d = 0; d < NSHARES; d++) {
        out[d * stride] = randomness[d][0];
    }
#elif N == 953
    b_and(out, stride, &randomness[0][0], L, &randomness[0][1],
               L);
    b_not(out, stride, out, stride);
    b_and(out, stride, &randomness[0][2], L, out, stride);
    b_and(out, stride, &randomness[0][3], L, out, stride);
    b_not(out, stride, out, stride);
    b_and(out, stride, &randomness[0][4], L, out, stride);
#else
    // CAUTION: default 
    for (size_t d = 0; d < NSHARES; d++) {
        out[d * stride] = randomness[d][0];
    }
#endif
#endif
}

/**
 * Samples a polnomial candidate and returns its weight. The caller must check whether the weight is 
 * correct or not and restart the algorithm with new randomness accordingly (or not).
 *
 * The caller must ensure that each (unmasked) value in the randomness array does not exceed 2**L-1.
 */
uint32_t sample_cmp(uint32_t out[NSHARES][(N + 31) / 32], uint32_t randomness[NSHARES][N], uint32_t thresh) {
  uint32_t tmp[NSHARES], accum[NSHARES] = {0};
  uint32_t masked_thresh[NSHARES], masked_lsb[NSHARES];
  uint32_t weight;
  b_mask(masked_thresh, 1, thresh);
  b_mask(masked_lsb, 1, 1);
  for (size_t i = 0; i < N; i++)
  {
#if defined(USE_MASKED_ISA) && !defined(USE_MASKED_EXT)
    asm volatile (
    "lw a0, (%[in0])\n" //rand
    "lw a1, (%[in1])\n"
    "lw a2, (%[in2])\n" //thres
    "lw a3, (%[in3])\n"

    "mask.b.sub (a5,a4),(a3,a2),(a1,a0)\n"
    "mask.b.srli (a5,a4),(a5,a4),31\n"

    "lw a6, (%[ac0])\n"
    "lw a7, (%[ac1])\n"
    "mask.b.add (a7,a6),(a7,a6),(a5,a4)\n"
    "lw a0, (%[out0])\n"
    "lw a1, (%[out1])\n"
    "mask.b.slli (a1,a0),(a1,a0), 1\n"
    "lw a2, (%[lsb0])\n" 
    "lw a3, (%[lsb1])\n"
    "mask.b.and (a5,a4),(a5,a4),(a3,a2)\n"
    "mask.b.ior (a1,a0),(a1,a0),(a5,a4)\n"
    "sw a0, (%[out0])\n"
    "sw a1, (%[out1])\n"
    :
    : [in0] "r" (&randomness[0][i]), 
      [in1] "r" (&randomness[0][i]+N), 
      [in2] "r" (masked_thresh), 
      [in3] "r" (masked_thresh+1),
      [ac0] "r" (accum), 
      [ac1] "r" (accum+1),
      [out0] "r" (&out[0][i/32]),
      [out1] "r" (&out[0][i/32]+(N+31)/32),
      [lsb0] "r" (&masked_lsb),
      [lsb1] "r" (&masked_lsb+1)
    : "a7", "a6", "a5", "a4", "a3", "a2", "a1", "a0", "memory");
#elif defined USE_MASKED_EXT && defined USE_MASKED_ISA
    asm volatile (
    "lw a0, (%[in0])\n" //rand
    "lw a1, (%[in1])\n"
    "lw a2, (%[in2])\n" //thres
    "lw a3, (%[in3])\n"
    "mask.b.cmpgt (a5,a4),(a1,a0),(a3,a2)\n"
    "lw a6, (%[ac0])\n"
    "lw a7, (%[ac1])\n"
    "mask.b.add (a7,a6),(a7,a6),(a5,a4)\n"
    "lw a0, (%[out0])\n"
    "lw a1, (%[out1])\n"
    "mask.b.slli (a1,a0),(a1,a0), 1\n"
    "lw a2, (%[lsb0])\n" 
    "lw a3, (%[lsb1])\n"
    "mask.b.and (a5,a4),(a5,a4),(a3,a2)\n"
    "mask.b.ior (a1,a0),(a1,a0),(a5,a4)\n"
    "sw a0, (%[out0])\n"
    "sw a1, (%[out1])\n"
    :
    : [in0] "r" (&randomness[0][i]), 
      [in1] "r" (&randomness[0][i]+N), 
      [in2] "r" (masked_thresh), 
      [in3] "r" (masked_thresh+1),
      [ac0] "r" (accum), 
      [ac1] "r" (accum+1),
      [out0] "r" (&out[0][i/32]),
      [out1] "r" (&out[0][i/32]+(N+31)/32),
      [lsb0] "r" (&masked_lsb),
      [lsb1] "r" (&masked_lsb+1)
    : "a7", "a6", "a5", "a4", "a3", "a2", "a1", "a0", "memory");
#else
    b_cmpg(tmp, 1, &randomness[0][i], N, masked_thresh, 1);
    b_add(accum, 1, accum, 1, tmp, 1);
    b_sll1(&out[0][i/32], (N+31)/32, &out[0][i/32], (N+31)/32);
    b_and(tmp, 1, tmp, 1, masked_lsb, 1);
    b_or(&out[0][i/32], (N+31)/32, &out[0][i/32], (N+31)/32, tmp, 1);
#endif 
  }
  b_unmask(&weight, accum, 1);
  return weight;
}

/**
 * Samples a polnomial candidate and returns its weight. The caller must check whether the weight is 
 * correct or not and restart the algorithm with new randomness accordingly (or not).
 */
uint8_t bs_sample_cmp(uint32_t out[NSHARES][(N + 31) / 32], uint32_t randomness[PAD32(N)/32][NSHARES][L]) {
    // sample the candidate
    for (size_t i = 0; i < (N + 31) / 32; i++) {
        sample_fixed_32(&out[0][i], (N + 31) / 32, randomness[i]);
    }

#ifdef USE_MASKED_ISA
    uint32_t weight = 0;
    uint32_t out_weight[NSHARES] = {0};
    uint32_t one[NSHARES], coeffs_shifted[NSHARES];
    b_mask(one, 1, 1);
    for (size_t i = 0; i < N; i++)
    {
      uint32_t tmp[NSHARES];
      if ((i%32) == 0)
      {
        b_and(tmp, 1, &out[0][i/32], (N+31)/32, one, 1);
      } else if ((i%32) == 1) {
        b_srl1(coeffs_shifted, 1, &out[0][i/32], (N+31)/32);
        b_and(tmp, 1, coeffs_shifted, 1, one, 1);
      } else {
        b_srl1(coeffs_shifted, 1, coeffs_shifted, 1);
        b_and(tmp, 1, coeffs_shifted, 1, one, 1);
      }
      b_add(out_weight, 1, out_weight, 1, tmp, 1);
    }
    b_unmask(&weight, out_weight, 1);
#else
    uint32_t out_arith_bs[NSHARES][8], out_arith_bs_x4[32], negA0[32];
    size_t bsoffset = 0;
    uint16_t coefcnt[NSHARES] = {0};
    uint8_t weight = 0;
    uint8_t out_weight[NSHARES] = {0};
    // convert to arithmetic shares mod 256 and accumulate weight
    for (size_t i = 0; i < (N + 31) / 32; i++) {
        // B2A conversion
        uint32_t B1[2], carry[NSHARES];

        B1[1] = randomint();
        negA0[bsoffset + 0] = randomint();
        B1[0] = negA0[bsoffset + 0] ^ B1[1];
        sechalfadd(NSHARES, carry, 1, &out_arith_bs[0][0], 8, &out[0][i],
                   (N + 31) / 32, B1, 1);
        for (size_t j = 1; j < 7; j++) {
            B1[1] = randomint();
            negA0[bsoffset + j] = randomint();
            B1[0] = negA0[bsoffset + j] ^ B1[1];
            sechalfadd(NSHARES, carry, 1, &out_arith_bs[0][j], 8, carry, 1, B1,
                       1);
        }
        B1[1] = randomint();
        negA0[bsoffset + 7] = randomint();
        B1[0] = negA0[bsoffset + 7] ^ B1[1];
        b_xor(&out_arith_bs[0][7], 8, carry, 1, B1, 1);

        // now: each slice in out_arith_bs[0]^out_arith_bs[1] is one arithmetic
        // share; 256 minus each slice in negA0 is the other share
        for (size_t j = 0; j < 8; j++) {
            out_arith_bs_x4[bsoffset + j] = out_arith_bs[0][j];
            for (size_t n = 1; n < NSHARES; n++) {
                out_arith_bs_x4[bsoffset + j] ^= out_arith_bs[n][j];
            }
        }

        bsoffset += 8;
        if (bsoffset == 32 || i == ((N + 31) / 32 - 1)) {
            // un-bitslice
            bsoffset = 0;

            transpose32(out_arith_bs_x4);
            for (size_t j = 0; j < 4; j++) {
                for (size_t k = 0; k < 32; k++) {
                    if (coefcnt[0] == N)
                        break;
                    out_weight[0] += (out_arith_bs_x4[k] >> (j * 8)) & 0xff;
                    coefcnt[0] += 1;
                }
            }

            transpose32(negA0);
            for (size_t j = 0; j < 4; j++) {
                for (size_t k = 0; k < 32; k++) {
                    if (coefcnt[1] == N)
                        break;
                    out_weight[1] -= (negA0[k] >> (j * 8)) & 0xff;
                    coefcnt[1] += 1;
                }
            }
        }
    }

    // unmask accumulated weight and store to `weight` variable, then this is
    // repeated until the correct weight is found
    weight = (out_weight[0] + out_weight[1]) % 256;
#endif
    return weight;
}
