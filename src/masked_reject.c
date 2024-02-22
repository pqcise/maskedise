#include "core_ops.h"
#include "gadgets.h"
#include "masked.h"
#include "masked_fwsampling.h"
#include <stddef.h>
#include <stdint.h>

/////////////////////////////////////////////////////////////////////////////////////////////
// Simple Rejection
// /////////////////////////////////////////////////////////////////////////
/*
void localsend(char byte) {
    volatile uint32_t* uart_tx = 0x40001004;
    *uart_tx = byte & 0xFF;
    for(int i = 0; i < 100; i++) {
        asm volatile ("nop");
    }
}
*/

// expects rand filled with randomness
int sample_reject(uint32_t p[NSHARES][W], uint32_t rand[NSHARES][B_PADDED]) {
    uint32_t i = 0;
    uint32_t r = 0;
    uint32_t j;
    uint32_t collision;
    uint32_t n[NSHARES];
    uint32_t res[NSHARES];
    uint32_t bitmask[NSHARES]; 
    

    b_mask(n, 1, N);
    b_mask(bitmask, 1, (1<<LOG_N)-1);

    while (i < W) {
        if (r >= B_PADDED)
            return 1; // very unlikely event of not enough randomness

        // zero out upper bits of rand:
        b_and(&rand[0][r], B_PADDED, &rand[0][r], B_PADDED, bitmask, 1);
        /*
        uint32_t testa;
        b_unmask(&testa, &rand[0][r], B_PADDED);
        */
        // check if N is greater than rand:
        b_cmpg(res, 1, &rand[0][r], B_PADDED, n, 1);
        b_unmask(res, res, 1);

        /*
        if(N > testa && res[0] != 1) {
            localsend(0xDE);
            localsend(0xA1);
        }
        */

        if((res[0] &1) != 1) {
            r++;
            /*
            for(int k = 0; k < 4; k++) {
                localsend((testa >> k*8) & 0xFF);
                for (int i = 0; i < 100; i++) { 
                        __asm__("nop");
                }
            }
            localsend(0xC0);
            localsend(0xCA);
            */
            continue;
        }

        // check for collisions:
        collision = 0;
        for(j=0; j < i; j++) {
            b_cmpeq(res, 1, &rand[0][r], B_PADDED, &p[0][j], W);
            b_unmask(res, res, 1);
            if(res[0] &1) {
                collision = 1;
                break; 
            }           
        }
        if(collision) {
            r++;
            continue;
        }

        // copy valid new value:
        copy_sharing(NSHARES, &p[0][i], W, &rand[0][r], B_PADDED);
        r++;
        i++;
    }
    return 0;
}

// expects rand filled with randomness
int bs_sample_reject(uint32_t p[NSHARES][W_PADDED], uint32_t rand[NSHARES][B_PADDED]) {
    const uint32_t n_not_pow2 = (1<<LOG_N != N);
    uint32_t i = 0;
    uint32_t c = 0;
    uint32_t r = 0;
    uint32_t index = 0;
    uint32_t r32[NSHARES][LOG_N];
    uint32_t n32[NSHARES][32];
    uint32_t valid[NSHARES];
    uint32_t res[NSHARES];

    // init N:
    for (int j = 0; j < 32; j++) {
        b_mask(&n32[0][j], 32, N);
    }
    transpose32_shares(&n32[0][0], 32);

    // init p:
    for (int j = 0; j < W_PADDED; j++) {
        b_mask(&p[0][j], W_PADDED, 0);
    }
    transpose32_array(&p[0][0], W_PADDED);

    while (i < W) {
        // compare next 32 random values
        if (r >= B_PADDED)
            return 1; // very unlikely event of not enough randomness
        transpose32_shares(&rand[0][r], B_PADDED);

        // If N is a power of two, a log(N) rand is always < N and thus valid
        // the comparison is thus not necessary and would actually require log(N)+1 bits.
        if(n_not_pow2) {
            bscmpg(valid, &rand[0][r], B_PADDED, &n32[0][0], 32);
            b_unmask(valid, valid, 1);
        }

        for (int j = 0; j < 32; j++) { // batch of 32 r
            if (((valid[0] >> j) & 1) == 0 && n_not_pow2)
                continue;              // r_j is >= N -> next r

            // prepare 32 times r in bitsliced domain for parallel collision check:
            for (int s = 0; s < NSHARES; s++) {
                for (int k = 0; k < LOG_N; k++) {
                    r32[s][k] = (rand[s][r + k] >> j) &
                                1; // copy r_j (bitsliced domain)
                    r32[s][k] =
                        ~(r32[s][k] + 0xffffffff); // extend single bit to 32bit
                }
            }

            for (c = 0; c < i; c += 32) {
                masked_bs_eq(res, 1, &r32[0][0], LOG_N, &p[0][c], W_PADDED,
                             LOG_N);  // check for collision
                b_unmask(res, res, 1);

                index = i & 31; // mod 32 -> bit within 32bit register where r
                                // gets stored bitsliced

                if (c + 31 > i - 1)         // we included comparisons with the
                                            // initialized 0 above i-1
                    res[0] &= 0xffffffff >>
                              (32 - index); // mask out these comparisons

                if (res[0] != 0)
                    break; // collision found
            }
            if (i > 0 && res[0] != 0)
                continue; // collision found -> next r

            // no collision:
            if (c > i)
                c -= 32; // c now points to nearest multiple of 32 <= i;

            for (int s = 0; s < NSHARES; s++) {
                for (int b = 0; b < LOG_N; b++) {
                    p[s][c + b] &= 0xffffffff ^ (1 << index); // set bit to zero
                    p[s][c + b] |= (r32[s][b] & 1)
                                   << index; // set bit according to r
                }
            }

            i++;
            if (i >= W)
                break; // stop 32er batch
        }
        r += 32;
    }
    transpose32_array(&p[0][0], W_PADDED); // unbitslice
    return 0;
}
