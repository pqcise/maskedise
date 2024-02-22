#include "frodo_sampler.h"
#include "core_ops.h"
#include "gadgets.h"

#if (FRODO_N == 1344)
static const uint16_t CDF_TABLE[7] = {9142, 23462, 30338, 32361, 32725, 32765, 32767};
#define CDF_TABLE_LEN 7


#elif (FRODO_N == 976)
static const uint16_t CDF_TABLE[11] = {5638,  15915, 23689, 28571, 31116, 32217,
                                32613, 32731, 32760, 32766, 32767};
#define CDF_TABLE_LEN 11


#elif (FRODO_N == 640)
static const uint16_t CDF_TABLE[13] = {4643,  13363, 20579, 25843, 29227, 31145, 32103,
                                32525, 32689, 32745, 32762, 32766, 32767};
#define CDF_TABLE_LEN 13

#else
#error No suitable FRODO_N defined
#endif

static const uint8_t LOGTABLE[] = {1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5};

/**
 * Mask stride is always FRODO_N, data stride is always 1!
 *
 * For each Boolean-masked entry in samples, we do:
 * 1. use the LSB as sign bit
 * 2. use all other bits to compare with the constants in CDF_TABLE
 * 3. sum all comparison results and write back to the same destination
 *
 * usage: write input randomness to input array, call the function.
 */
void frodo_sampler(uint16_t samples[NSHARES * FRODO_N]) {
    size_t i, j, k;
    uint32_t lowmask[NSHARES], lsbmask[NSHARES], masked_CDF_TABLE[NSHARES][CDF_TABLE_LEN];
    b_mask(lowmask, 1, 0xffff);
    b_mask(lsbmask, 1, 1);
    for (i = 0; i < CDF_TABLE_LEN; i++)
    {
      b_mask(&masked_CDF_TABLE[0][i], CDF_TABLE_LEN, CDF_TABLE[i]);
    }
    for (i = 0; i < FRODO_N; i += 2) // process two samples at once
    {
        uint32_t rand[2][NSHARES], result[2][NSHARES], signs[2][NSHARES];

        // extract the low and the high randomness
        b_and(rand[0], 1, lowmask, 1, (uint32_t *)&samples[i], FRODO_N / 2);
        masked_shiftr(rand[1], 1, (uint32_t *)&samples[i], FRODO_N / 2, 16);

        for (j = 0; j < 2; j++) // sample for both
        {
            uint32_t tmp[NSHARES];
            b_and(signs[j], 1, lsbmask, 1, rand[j], 1); // extract sign
            b_srl1(rand[j], 1, rand[j], 1); // shift right to obtain random
                                               // value that is compared to CDF
            b_cmpg(result[j], 1, rand[j], 1, &masked_CDF_TABLE[0][0],
                   CDF_TABLE_LEN); // compare to first value
            for (k = 1; k < CDF_TABLE_LEN; k++) {
#if defined(USE_MASKED_ISA) && !defined(USE_MASKED_EXT)
            asm volatile (
            "lw a0, (%[in0])\n"
            "lw a1, (%[in1])\n"
            "lw a2, (%[in2])\n"
            "lw a3, (%[in3])\n"

            "mask.b.sub (a5,a4),(a1,a0),(a3,a2)\n"
            "mask.b.srli (a5,a4),(a5,a4),31\n"
            
            "lw a0, (%[res0])\n"
            "lw a1, (%[res1])\n"
            "mask.b.add (a3,a2),(a5,a4),(a1,a0)\n"
            "sw a2, (%[res0])\n"
            "sw a3, (%[res1])\n"
            :
            : [in0] "r" (&masked_CDF_TABLE[0][k]), 
              [in1] "r" (&masked_CDF_TABLE[0][k]+CDF_TABLE_LEN), 
              [in2] "r" (&rand[j]), 
              [in3] "r" (&rand[j]+1),
              [res0] "r" (&result[j]),
              [res1] "r" (&result[j]+1)
            : "a5", "a4", "a3", "a2", "a1", "a0", "memory");
#elif defined USE_MASKED_EXT && defined USE_MASKED_ISA
            asm volatile (
            "lw a0, (%[in0])\n"
            "lw a1, (%[in1])\n"
            "lw a2, (%[in2])\n"
            "lw a3, (%[in3])\n"
            "mask.b.cmpgt (a5,a4),(a3,a2),(a1,a0)\n"
            "lw a0, (%[res0])\n"
            "lw a1, (%[res1])\n"
            "mask.b.add (a3,a2),(a5,a4),(a1,a0)\n"
            "sw a2, (%[res0])\n"
            "sw a3, (%[res1])\n"
            :
            : [in0] "r" (&masked_CDF_TABLE[0][k]), 
              [in1] "r" (&masked_CDF_TABLE[0][k]+CDF_TABLE_LEN), 
              [in2] "r" (&rand[j]), 
              [in3] "r" (&rand[j]+1),
              [res0] "r" (&result[j]),
              [res1] "r" (&result[j]+1)
            : "a5", "a4", "a3", "a2", "a1", "a0", "memory");
#else
                b_cmpg(tmp, 1, rand[j], 1, &masked_CDF_TABLE[0][k],
                       CDF_TABLE_LEN); // compare to all other values
                b_add(result[j], 1, result[j], 1, tmp, 1); // accumulate
#endif
            }
            // here, result[j] holds the sample value without sign
            b_not(tmp, 1, signs[j], 1);
            b_add(tmp, 1, tmp, 1, lsbmask, 1); // turn 0->0 and 1->all-ones
            // signed sample = ((-sign) ^ result) + sign = (tmp ^ result) + sign
            b_xor(result[j], 1, result[j], 1, tmp, 1);
            b_add(result[j], 1, result[j], 1, signs[j], 1);

            // copy to output, including conversion to 16 bit
            for (k = 0; k < NSHARES; k++) {
                samples[i + j + k * FRODO_N] = result[j][k];
            }
        }
    }
}

/**
 * 16 bit slices represent 32 input values, then we do the same as above
 *
 * Usage: put the input randomness bitsliced into the 3d-array samples, then call the function.
 * Output will be written to the same array.
 */
void frodo_sampler_bs(uint32_t samples[NSHARES][(FRODO_N + 31)/32][16]) {
    size_t i, j, k;
    uint32_t tmp[NSHARES] = {0};

    for (i = 0; i < (FRODO_N + 31) / 32; i++) 
    {
        uint32_t result[NSHARES][16] = {0};
        for (j = 0; j < CDF_TABLE_LEN; j++) {
            uint32_t cmpres[NSHARES] = {0};
            for (k = 0; k < 15;
                 k++) // iterate over the bits that are significant for
                      // comparison
            {
                // in the following, we access offset k+1 because the LSB is the
                // sign
                b_not(tmp, 1,
                      &samples[0][i][k+1],
                      ((FRODO_N+31)/32)*16);
                if (CDF_TABLE[j] & (1 << k)) {
                    b_or(cmpres, 1, cmpres, 1, tmp, 1);
                } else {
                    b_and(cmpres, 1, cmpres, 1, tmp, 1);
                }
            }
            if (j == 0) {
                b_not(result[0], 16, cmpres, 1);
            } else {
                b_not(cmpres, 1, cmpres, 1);
                // ripple-carry adder LOGTABLE[j] bit + 1 bit with tmp as carry
                b_and(tmp, 1, result[0], 16, cmpres, 1);
                b_xor(result[0], 16, result[0], 16, cmpres, 1);
                for (k = 1; k <= (size_t) LOGTABLE[j]+1; k++) {
                    uint32_t tmp2[2] = {0};
                    b_xor(tmp2, 1, &result[0][k], 16, tmp, 1);
                    b_and(tmp, 1, &result[0][k], 16, tmp, 1);
                    copy_sharing(NSHARES, &result[0][k], 16, tmp2, 1);
                }
            }
        }
        // negate conditionally (condition is samples[*][i][0])
        for (k = 0; k < 16; k++)
        {
          b_xor(&result[0][k], 16, &result[0][k], 16, &samples[0][i][0], ((FRODO_N+31)/32)*16);
        }
        b_and(tmp, 1, result[0], 16, &samples[0][i][0], ((FRODO_N+31)/32)*16);
        b_xor(&samples[0][i][0], ((FRODO_N+31)/32)*16, result[0], 16, &samples[0][i][0], ((FRODO_N+31)/32)*16);
        for (k = 1; k < 16; k++) {
            b_xor(&samples[0][i][k], ((FRODO_N+31)/32)*16, &result[0][k], 16, tmp, 1);
            b_and(tmp, 1, &result[0][k], 16, tmp, 1);
        }
    }
}
