#include "gadgets.h"
#include "masked.h"
#include "masked_fwsampling.h"
#include <stddef.h>
#include <stdint.h>

void masked_binary_to_trinary(uint32_t signs[NSHARES][(N + 31) / 32],
                              const uint32_t poly[NSHARES][(N + 31) / 32]) {
    size_t i, d;
    uint32_t tmp[NSHARES];
    for (i = 0; i < (N + 31) / 32; i++) {
        for (d = 0; d < NSHARES; d++) {
            tmp[d] = randomint();
        }
        masked_and(NSHARES, &signs[0][i], (N + 31) / 32, &poly[0][i],
                   (N + 31) / 32, tmp, 1);
    }
}

void masked_binary_to_trinary_fixed(uint32_t signs[NSHARES][(N + 31) / 32],
                                    const uint32_t poly[NSHARES][(N + 31) / 32],
                                    const uint8_t target_weight) {
    uint8_t weight = 0;
    do {
        uint32_t out_arith_bs[NSHARES][8], out_arith_bs_x4[32], negA0[32];
        size_t bsoffset = 0;
        uint8_t out_weight[NSHARES] = {0};
        uint16_t coefcnt[NSHARES] = {0};

        masked_binary_to_trinary(signs, poly);

        // convert to arithmetic shares mod 256 and accumulate weight
        for (size_t i = 0; i < (N + 31) / 32; i++) {
            // B2A conversion
            uint32_t B1[2], carry[NSHARES];

            B1[1] = randomint();
            negA0[bsoffset + 0] = randomint();
            B1[0] = negA0[bsoffset + 0] ^ B1[1];
            sechalfadd(NSHARES, carry, 1, &out_arith_bs[0][0], 8, &signs[0][i],
                       (N + 31) / 32, B1, 1);
            for (size_t j = 1; j < 7; j++) {
                B1[1] = randomint();
                negA0[bsoffset + j] = randomint();
                B1[0] = negA0[bsoffset + j] ^ B1[1];
                sechalfadd(NSHARES, carry, 1, &out_arith_bs[0][j], 8, carry, 1,
                           B1, 1);
            }
            B1[1] = randomint();
            negA0[bsoffset + 7] = randomint();
            B1[0] = negA0[bsoffset + 7] ^ B1[1];
            masked_xor(NSHARES, &out_arith_bs[0][7], 8, carry, 1, B1, 1);

            // now: each slice in out_arith_bs[0]^out_arith_bs[1] is one
            // arithmetic share; 256 minus each slice in negA0 is the other
            // share
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

        // unmask accumulated weight and store to `weight` variable, then this
        // is repeated until the correct weight is found
        weight = (out_weight[0] + out_weight[1]) % 256;
        bsoffset = 0;

    } while (weight != target_weight);
}
