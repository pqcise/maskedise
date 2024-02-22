#include "core_ops.h"
#include "gadgets.h"
#include "masked.h"
#include "masked_fwsampling.h"
#include <stddef.h>
#include <stdint.h>
#include <string.h>

// expects randomness in p for first W entries
static void bs_fy_part1(uint32_t p[NSHARES][W_PADDED]) {
    uint32_t t[NSHARES][64];
    union U64 x;

    memset((uint8_t*) &t[0][0], 0, sizeof(uint32_t) * NSHARES * 64);

    // component 1: sampling int out of range: p[i] = i + rand(n-1-i);
    for (int i = 0; i < W_PADDED; i += 32) {
        // assume we receive our random values boolean bitsliced masking domain
        for (int j = 0; j < 32; j++) {
            copy_sharing(NSHARES, &t[0][j], 64, &p[0][i + j], W_PADDED);
        }

        secb2a_32to48(NSHARES, &t[0][0], 64,
                      1); // transform to mod 2^48 arithmetic domain;
        transpose32_array(&t[0][0], 64); // unbitslice

        for (int j = 0; j < 32; j++) {
            for (int s = 0; s < NSHARES; s++) {
                x.v32[0] = t[s][j + 0];  // lower 32bit
                x.v32[1] = t[s][j + 32]; // upper 32bit
                x.v64 = x.v64 * ((uint64_t)N - 1 - i - j);
                t[s][j + 0] = x.v32[0];
                t[s][j + 32] = x.v32[1] & 0x0000FFFF; // mod 2^48
            }
            t[0][j + 32] =
                (t[0][j + 32] + i + j) &
                0x0000FFFF; // add i to the upper half, compute mod 2^48
        }
        // a2b:
        transpose32_array(&t[0][0], 64); // bitslice
        seca2b(NSHARES, 48, &t[0][0], 64, 1);

        for (int j = 0; j < 32; j++) {
            // copy back upper 32bits:
            copy_sharing(NSHARES, &p[0][i + j], W_PADDED, &t[0][j + 32], 64);
        }
    }
}

void bs_sample_fisher_yates(uint32_t p[NSHARES][W_PADDED]) {

    // component 1: sampling int out of range: p[i] = i + rand(n-1-i);
    bs_fy_part1(p);

    // component 2: if collision then exchange bitsliced
    int start_j32 = W_PADDED - 32; // init with rightmost index in 32er steps
    int i32 = W_PADDED - 32;
    uint32_t mask;

    uint32_t bitmask[NSHARES];
    uint32_t bs_i[NSHARES][32];  // 32 times i
    uint32_t bs_pi[NSHARES][32]; // 32 times p[i]
    uint32_t res[NSHARES];       // result of comparison

    for (int i = W - 2; i >= 0; i--) {
        start_j32 = start_j32 > i + 1
                        ? start_j32 - 32
                        : start_j32;    // index of 32er block containing j
        i32 = i32 > i ? i32 - 32 : i32; // index of 32er block containing i

        // prepare 32times index i:
        for (int j = 0; j < 32; j++) {
            b_mask(&bs_i[0][j], 32, i);
        }
        transpose32_shares(&bs_i[0][0], 32);

        // prepare 32times p[i]:
        // unbitslice to avoid bitfiddeling:
        transpose32_shares(&p[0][i32], W_PADDED);
        for (int j = 0; j < 32; j++) {
            // copy p[i] 32 times for comparison
            copy_sharing(NSHARES, &bs_pi[0][j], 32, &p[0][i], W_PADDED);
        }
        // bitslice again:
        transpose32_shares(&p[0][i32], W_PADDED);
        transpose32_shares(&bs_pi[0][0], 32);

        for (int j32 = start_j32; j32 + 31 < W_PADDED; j32 += 32) {
            masked_bs_eq(res, 1, &p[0][j32], W_PADDED, &bs_pi[0][0], 32, LOG_N);

            int j = i + 1;
            int diff = j - j32;
            mask = 0xFFFFFFFF;
            if (diff > 0) {
                // rightmost diff bits are zero -> donttouch these indizes
                mask = mask << diff;
                b_mask(bitmask, 1, mask);
                b_and(res, 1, res, 1, bitmask, 1); 
            }

            masked_bs_sel(res, 1, &p[0][j32], W_PADDED, &bs_i[0][0], 32, LOG_N);
        }
    }

    for (int i = 0; i < W_PADDED; i += 32) {
        transpose32_shares(&p[0][i], W_PADDED);
    }
}

// expects randomness in p for first W entries
void sample_fisher_yates(uint32_t p[NSHARES][W_PADDED]) {
    uint32_t res[NSHARES];
    uint32_t i_masked[NSHARES];

    // component 1: sampling int out of range: p[i] = i + rand(n-1-i);
    bs_fy_part1(p);

    // component 2: if collision then exchange
    for (int i = 0; i < W_PADDED; i += 32) {
        transpose32_shares(&p[0][i], W_PADDED); //unbitslice
    }

    for (int i = W - 2; i >= 0; i--) {
        b_mask(i_masked, 1, i);
        for(int j = i + 1; j < W; j++) {
#if defined(USE_MASKED_ISA) && !defined(USE_MASKED_EXT)
            asm volatile (
            "lw a0, (%[in0])\n"
            "lw a1, (%[in1])\n"
            "lw a2, (%[in2])\n"
            "lw a3, (%[in3])\n"

            "mask.b.sub (a5,a4),(a3,a2),(a1,a0)\n"
            "mask.b.srli (a5,a4),(a5,a4),31\n"
            "mask.b.sub (a7,a6), (a1,a0),(a3,a2)\n"
            "mask.b.srli (a7,a6), (a7,a6), 31\n"
            "mask.b.ior (a5,a4),(a5,a4),(a7,a6)\n"
            "mask.b.not (a5,a4),(a5,a4)\n"
            "li t0, 1\n"
            "mask.b.mask (a7,a6),t0\n"
            "mask.b.and (a5,a4),(a5,a4),(a7,a6)\n"

            "lw a0, (%[im0])\n"
            "lw a1, (%[im1])\n"

            "mask.b.xor (a7,a6),(a3,a2),(a1,a0)\n"
            "mask.b.and (t6,t5),(a7,a6),(a5,a4)\n"
            "mask.b.xor (a3,a2),(a3,a2),(t6,t5)\n"

            "sw a2, (%[in2])\n"
            "sw a3, (%[in3])\n"
            :
            : [in0] "r" (&p[0][j]), 
            [in1] "r" (&p[0][j]+W_PADDED), 
            [in2] "r" (&p[0][i]), 
            [in3] "r" (&p[0][i]+ W_PADDED),
            [im0] "r" (i_masked), 
            [im1] "r" (i_masked+1)
            : "t6", "t5", "t0", "a7", "a6", "a5", "a4", "a3", "a2", "a1", "a0", "memory");
#elif defined USE_MASKED_EXT && defined USE_MASKED_ISA
            asm volatile (
            "lw a0, (%[in0])\n"
            "lw a1, (%[in1])\n"
            "lw a2, (%[in2])\n"
            "lw a3, (%[in3])\n"
            "mask.b.cmpeq (a5,a4),(a3,a2),(a1,a0)\n"
            "lw a0, (%[im0])\n"
            "lw a1, (%[im1])\n"
            "mask.b.cmov (a3,a2),(a1,a0),(a5,a4)\n"
            "sw a2, (%[in2])\n"
            "sw a3, (%[in3])\n"
            :
            : [in0] "r" (&p[0][j]), 
            [in1] "r" (&p[0][j]+W_PADDED), 
            [in2] "r" (&p[0][i]), 
            [in3] "r" (&p[0][i]+ W_PADDED),
            [im0] "r" (i_masked), 
            [im1] "r" (i_masked+1)
            : "a5", "a4", "a3", "a2", "a1", "a0", "memory");
#else
           b_cmpeq(res, 1, &p[0][j], W_PADDED, &p[0][i], W_PADDED);
           b_cmov(&p[0][j], W_PADDED, i_masked, 1, res, 1); 
#endif
        }
    }
}

