#include "bike_rotr.h"
#include "core_ops.h"
#include "gadgets.h"

// expects a 3*R_QWORDS sized array as input and qw_num < R_QWORDS
void rotr_big(uint32_t in[NSHARES][3 * R_QWORDS], uint32_t qw_num[NSHARES]) {
    uint32_t res[NSHARES];
    uint32_t shift_amount, tmp = R_QWORDS_HALF_LOG2;
    uint32_t idx = R_QWORDS_HALF_LOG2;
    uint32_t a[NSHARES];
    uint32_t b[NSHARES];

    for (shift_amount = 0; tmp != 1; tmp >>= 1, shift_amount++)
        ;

    for (idx = R_QWORDS_HALF_LOG2; idx >= 1; idx >>= 1) {
        // Get bit_i (i=shift_amount) from the secret.
        // We use a for loop, to be able to use the masked instruction for fixed shift amounts.
        // The inefficiency is not that relevant here.
        b_mask(res, 1, 0);
        b_xor(res, 1, res, 1, qw_num, 1);
        for(uint32_t i=0; i<shift_amount; i++)
            b_srl1(res, 1, res, 1); 

        // Rotate R_QWORDS words and another idx words,
        // as needed by the next iteration

        for (size_t offset = 0; offset < idx; offset++) {
            a[0] = in[0][0 + offset];
            a[1] = in[1][0 + offset];

            // keep this as c code:
            for (size_t i = 0; i + offset < (R_QWORDS + idx); i+=2*idx) {
            // inline asm again below:

                b[0] = in[0][i + offset + idx];
                b[1] = in[1][i + offset + idx]; 
                b_cmov(a, 1, b, 1, res, 1); 
                in[0][i + offset] = a[0];
                in[1][i + offset] = a[1]; 

                a[0] = in[0][i + offset + 2*idx];
                a[1] = in[1][i + offset + 2*idx];  
                b_cmov(b, 1, a, 1, res, 1); 
                in[0][i + offset + idx] = b[0];
                in[1][i + offset + idx] = b[1]; 
            }
        }
        shift_amount--;
    }
}