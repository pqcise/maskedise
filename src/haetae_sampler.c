#include "haetae_sampler.h"
#include "core_ops.h"
#include "gadgets.h"

static const uint16_t CDF_TABLE[64] = {
 3266,  6520,  9748, 12938, 16079, 19159, 22168, 25096,
27934, 30674, 33309, 35833, 38241, 40531, 42698, 44742,
46663, 48460, 50135, 51690, 53128, 54454, 55670, 56781,
57794, 58712, 59541, 60287, 60956, 61554, 62085, 62556,
62972, 63337, 63657, 63936, 64178, 64388, 64569, 64724,
64857, 64970, 65066, 65148, 65216, 65273, 65321, 65361,
65394, 65422, 65444, 65463, 65478, 65490, 65500, 65508,
65514, 65519, 65523, 65527, 65529, 65531, 65533, 65534
};
#define CDF_TABLE_LEN 64
#define CDF_NUM_SIGNIFICANT_BITS 16

static const uint8_t LOGTABLE[] = {
  1, 2, 2, 3, 3, 3, 3, 4, 
  4, 4, 4, 4, 4, 4, 4, 5, 
  5, 5, 5, 5, 5, 5, 5, 5, 
  5, 5, 5, 5, 5, 5, 5, 6, 
  6, 6, 6, 6, 6, 6, 6, 6, 
  6, 6, 6, 6, 6, 6, 6, 6, 
  6, 6, 6, 6, 6, 6, 6, 6, 
  6, 6, 6, 6, 6, 6, 6, 7
};

/**
 * Mask stride is always HAETAE_N, data stride is always 1!
 *
 * For each Boolean-masked entry in samples, we do:
 * 1. use the LSB as sign bit
 * 2. use all other bits to compare with the constants in CDF_TABLE
 * 3. sum all comparison results and write back to the same destination
 *
 * usage: write input randomness to input array, call the function.
 */
void haetae_sampler(uint16_t samples[NSHARES * HAETAE_N]) {
    size_t i, j, k;
    uint32_t lowmask[NSHARES], masked_CDF_TABLE[NSHARES][CDF_TABLE_LEN];
    b_mask(lowmask, 1, 0xffff);
    for (i = 0; i < CDF_TABLE_LEN; i++)
    {
      b_mask(&masked_CDF_TABLE[0][i], CDF_TABLE_LEN, CDF_TABLE[i]);
    }
    for (i = 0; i < HAETAE_N; i += 2) // process two samples at once
    {
        uint32_t rand[2][NSHARES], result[2][NSHARES];

        // extract the low and the high randomness
        b_and(rand[0], 1, lowmask, 1, (uint32_t *)&samples[i], HAETAE_N / 2);
        masked_shiftr(rand[1], 1, (uint32_t *)&samples[i], HAETAE_N / 2, 16);

        for (j = 0; j < 2; j++) // sample for both
        {
            uint32_t tmp[NSHARES];
            b_cmpg(result[j], 1, rand[j], 1, &masked_CDF_TABLE[0][0],
                   CDF_TABLE_LEN); // compare to first value
            for (k = 1; k < CDF_TABLE_LEN; k++) {
                b_cmpg(tmp, 1, rand[j], 1, &masked_CDF_TABLE[0][k],
                       CDF_TABLE_LEN); // compare to all other values
                b_add(result[j], 1, result[j], 1, tmp, 1); // accumulate
            }

            // copy to output, including conversion to 16 bit
            for (k = 0; k < NSHARES; k++) {
                samples[i + j + k * HAETAE_N] = result[j][k];
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
void haetae_sampler_bs(uint32_t samples[NSHARES][(HAETAE_N + 31)/32][16]) {
    size_t i, j, k;
    uint32_t tmp[NSHARES] = {0};

    for (i = 0; i < (HAETAE_N + 31) / 32; i++) 
    {
        uint32_t result[NSHARES][7] = {0};
        for (j = 0; j < CDF_TABLE_LEN; j++) {
            uint32_t cmpres[NSHARES] = {0};
            for (k = 0; k < CDF_NUM_SIGNIFICANT_BITS;
                 k++) // iterate over the bits that are significant for
                      // comparison
            {
                b_not(tmp, 1,
                      &samples[0][i][k],
                      ((HAETAE_N+31)/32)*16);
                if (CDF_TABLE[j] & (1 << k)) {
                    b_or(cmpres, 1, cmpres, 1, tmp, 1);
                } else {
                    b_and(cmpres, 1, cmpres, 1, tmp, 1);
                }
            }
            if (j == 0) {
                b_not(result[0], 7, cmpres, 1);
            } else {
                b_not(cmpres, 1, cmpres, 1);
                // ripple-carry adder LOGTABLE[j] bit + 1 bit with tmp as carry
                b_and(tmp, 1, result[0], 7, cmpres, 1);
                b_xor(result[0], 7, result[0], 7, cmpres, 1);
                for (k = 1; k < (size_t) LOGTABLE[j]; k++) {
                    uint32_t tmp2[2] = {0};
                    b_xor(tmp2, 1, &result[0][k], 7, tmp, 1);
                    b_and(tmp, 1, &result[0][k], 7, tmp, 1);
                    copy_sharing(NSHARES, &result[0][k], 7, tmp2, 1);
                }
            }
        }
        for (j = 0; j < 7; j++) 
        {
          copy_sharing(NSHARES, &samples[0][i][j], ((HAETAE_N+31)/32)*16, &result[0][j], 7);
        }

        // zeroize remaining random input
        for (k = 0; k < NSHARES; k++)
        {
          for (j = 7; j < 16; j++)
          {
            samples[k][i][j] = 0;
          }
        }
    }
}
