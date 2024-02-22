#include "core_ops.h"
#include "gadgets.h"
#include "masked.h"
#include "masked_fwsampling.h"
#include "masked_representations.h"
#include <stddef.h>
#include <stdint.h>

/*
 * a[0] holds the low bits of the first share
 * ...
 * a[31] holds the MSB of the first share
 * a[PAD64(N)] holds the LSBs of the second share
 * ...
 * a[PAD64(N)+31] holds the MSBs of the second share
 * and so on
 */
static void bscmp(uint32_t *res, const uint32_t *a, const uint32_t stride_a,
                  const uint32_t *b, const uint32_t stride_b) {
    size_t i;
    uint32_t xor [NSHARES], tmp[NSHARES];

    b_mask(res, 1, 0);

    for (i = 2; i < 32; i++) // lowest 2bit do not belong to sorted value
    {
        b_xor(xor, 1, a + i, stride_a, b + i, stride_b);
        b_xor(tmp, 1, res, 1, b + i, stride_b);
        b_and(tmp, 1, tmp, 1, xor, 1);
        b_xor(res, 1, res, 1, tmp, 1);
    }
}

static void bsint32_MINMAX(uint32_t *a, uint32_t *b, uint32_t direction,
                           uint32_t stride_b) {
    size_t i;
    uint32_t res[NSHARES], xor[NSHARES], mux[NSHARES][32];
    // res = 1 iff b > a
    bscmp(res, a, PAD64(N), b, stride_b);

    // flip res based on sorting direction
    res[0] = res[0] ^ direction;
    res[0] = ~res[0]; // flip res for faster cmov

    // now swap conditionally based on res
    for (i = 0; i < 32; i++) {
        b_xor(xor, 1, &a[i], PAD64(N), &b[i], stride_b);
        b_and(&mux[0][i], 32, xor, 1, res, 1);
    }

    // now, mux will contain a if b > a, else b
    for (i = 0; i < 32; i++) {
        b_xor(&a[i], PAD64(N), &a[i], PAD64(N), &mux[0][i], 32);
        b_xor(&b[i], stride_b, &b[i], stride_b, &mux[0][i], 32);
    }
}

static int isPowerOfTwo(unsigned int x) { return ((x != 0) && !(x & (x - 1))); }

static int greatestPowerOfTwoLessThan(int n) {
    int k = 1;
    while (k < n)
        k = k << 1;
    return k >> 1;
}

static int greatestPowerOfTwoLessThanEqual(int n) {
    if (isPowerOfTwo(n))
        return n;
    else
        return greatestPowerOfTwoLessThan(n);
}

static uint32_t rotr32a(uint32_t x, uint32_t n) {
    return (x << (32 - n)) | (x >> n);
}

static uint32_t rotl32a(uint32_t x, uint32_t n) {
    return (x >> (32 - n)) | (x << n);
}

// converts two consecutive 32bit words for a minmax step in bitonic sort with
// p<32 so that 32 minmax operations can be done in parallel via bitslicing
static void bsconvertsub32(uint32_t *x, int current_distance_log,
                           int target_distance_log) {
    uint32_t t;

    const uint32_t mask[5] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF,
                              0x0000FFFF};

    for (int c = 0; c < 32; c++) {
        for (int i = current_distance_log - 1; i >= target_distance_log; i--) {
            t = rotr32a(x[c] & ~mask[i], 1 << i);
            x[c] = (x[c] & mask[i]) | rotl32a(x[c + 32] & mask[i], 1 << i);
            x[c + 32] = t | (x[c + 32] & ~mask[i]);
        }
    }
}

static void bsbackconvertsub32(uint32_t *x, int current_distance_log,
                               int target_distance_log) {
    uint32_t t;

    const uint32_t mask[5] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF,
                              0x0000FFFF};

    for (int c = 0; c < 32; c++) {
        for (int i = current_distance_log; i < target_distance_log; i++) {
            t = rotr32a(x[c] & ~mask[i], 1 << i);
            x[c] = (x[c] & mask[i]) | rotl32a(x[c + 32] & mask[i], 1 << i);
            x[c + 32] = t | (x[c + 32] & ~mask[i]);
        }
    }
}

static void bsconvertsub32poly(uint32_t *x, int n, int current_distance,
                               int target_distance) {
    int current_distance_log = 0;
    int target_distance_log = 0;

    if (target_distance >= 32 && current_distance >= 32)
        return;
    if (target_distance > 32)
        target_distance = 32;
    if (current_distance > 32)
        current_distance = 32;

    while (1 << current_distance_log < current_distance)
        current_distance_log++;
    while (1 << target_distance_log < target_distance)
        target_distance_log++;

    if (target_distance < current_distance) {
        for (int s = 0; s < NSHARES; s++) {
            for (int i = 0; i < n - 63; i += 64) {
                bsconvertsub32(&x[s * PAD64(N) + i], current_distance_log,
                               target_distance_log);
            }
        }
    } else if (target_distance > current_distance) {
        for (int s = 0; s < NSHARES; s++) {
            for (int i = 0; i < n - 63; i += 64) {
                bsbackconvertsub32(&x[s * PAD64(N) + i], current_distance_log,
                                   target_distance_log);
            }
        }
    }
}

// expects already bitsliced polynomial
static void masked_bitonic_merge(uint32_t *x, int n, int k, int direction) {
    uint32_t directionmask;
    int i, j, old_j = 32;

    for (j = k; j > 0; j = j >> 1) {
        bsconvertsub32poly(x, n, old_j,
                           j); // from/to minmax within one register

        if (j >= 32) {         // pairs do not share registers
            for (i = 0; i + j + 31 < n; i += 32) {
                int ij = i ^ j;
                if ((ij) > i) {
                    directionmask = !(i & k << 1) ? 0x00000000 : 0xFFFFFFFF;
                    if (direction)
                        directionmask = ~directionmask;
                    bsint32_MINMAX(&x[i], &x[ij], directionmask, PAD64(N));
                }
            }
        } else if (j <
                   32) { // pairs do share registers -> seperate them with
                         // bsconvertsub32poly function -> thus we need to
                         // operate on 64 consecutive elements simulatenously
            switch (k) {
            case 16:
                directionmask = 0xFFFF0000;
                break;
            case 8:
                directionmask = 0xFF00FF00;
                break;
            case 4:
                directionmask = 0xF0F0F0F0;
                break;
            case 2:
                directionmask = 0xCCCCCCCC;
                break;
            case 1:
                directionmask = 0xAAAAAAAA;
                break;
            default:
                directionmask = 0x00000000;
                break;
            }

            for (i = 0; i + 63 < n; i += 64) {
                if (k > 16) {
                    directionmask = !(i & k << 1) ? 0x00000000 : 0xFFFFFFFF;
                }
                if (direction)
                    directionmask = ~directionmask;
                bsint32_MINMAX(&x[i], &x[i + 32], directionmask, PAD64(N));
            }
        }
        old_j = j;
    }
    bsconvertsub32poly(x, n, old_j, 32); // from/to minmax within one register
}

// expects already bitsliced polynomial
static void masked_bitonic_sort_pow2(uint32_t *x, int n, int direction) {
    for (int k = 2; k <= n; k = 2 * k) {
        masked_bitonic_merge(x, n, k >> 1, direction);
    }
}

// Very efficient if n is power of two.
// Not the most efficient approach for bitonic sort with an n that is not a
// power of two, we cant reach nlogn because n=64 is the lowest we can get down
// to while keeping the code somewhat simple. Slightly more efficient variant
// would be possible though.
static void masked_bitonic_sort(uint32_t x[NSHARES][PAD64(N)]) {
    transpose32_array(&x[0][0], PAD64(N)); // bitslicing

    int pows2[16] = {0}; // hardcoded indirect limit for size of N
    int n = PAD64(N);
    int m = 0;
    int i = 0;

    while (n > 0) {
        m = greatestPowerOfTwoLessThanEqual(n);
        pows2[i] = m;
        i++;
        n = n - m;
    }

    uint32_t *b = &x[0][0] + PAD64(N) - pows2[--i];
    masked_bitonic_sort_pow2(b, pows2[i], 0);
    int n_merge = pows2[i];

    while (i > 0) {
        b = b - pows2[--i];
        masked_bitonic_sort_pow2(b, pows2[i], 1);
        n_merge += pows2[i];
        masked_bitonic_merge(b, n_merge, greatestPowerOfTwoLessThan(n_merge),
                             0);
    }

    transpose32_array(&x[0][0], PAD64(N)); // unbitscliing
}

// expects randomness in the upper 30bits for the first N entries in x
void bs_sample_sort(uint32_t x[NSHARES][PAD64(N)]) {
    int i;
    uint32_t one[NSHARES];
    uint32_t bitmask[NSHARES];

    b_mask(one, 1, 1);
    b_mask(bitmask, 1, 0xFFFFFFFC);

    for (i = 0; i < W; i++) {
        // set first W coeffs to 1:
        b_or(&x[0][i], PAD64(N), &x[0][i], PAD64(N), one, 1);
    }
    for (; i < N; i++) {
        // set remaining coeffs to 0:
        b_and(&x[0][i], PAD64(N), &x[0][i], PAD64(N), bitmask, 1);
    }
    for (; i < PAD64(N); i++) {
        // set padded entries to max value:
        b_mask(&x[0][i], PAD64(N), 0xFFFFFFFF);
    }

    masked_bitonic_sort(x);
}

static void masked_minmax(uint32_t *a, size_t a_stride, uint32_t *b, size_t b_stride) {
    uint32_t tmp[NSHARES];
    uint32_t res[NSHARES];
    uint32_t ab[NSHARES];
    uint32_t a_shifted[NSHARES];
    uint32_t b_shifted[NSHARES];
#if defined(USE_MASKED_ISA) && !defined(USE_MASKED_EXT)
    asm volatile (
    "li a0, 1\n"
    "mask.b.mask (a1,a0), a0\n" //a1,a0 = 1
    "lw a2, (%[a0])\n" //a3,a2 = a
    "lw a3, (%[a1])\n"
    "lw a4, (%[b0])\n" //a4,a5 = b
    "lw a5, (%[b1])\n"
    "mask.b.xor (a7,a6), (a3,a2), (a5,a4)\n" //a7,a6 = ab
    "mask.b.srli (t2,t1),(a3,a2),1\n" //t6,t5 shifted a
    "mask.b.srli (t4,t3),(a5,a4),1\n" //t4,t3 shifted b
    
    "mask.b.sub (t2,t1),(t4,t3),(t2,t1)\n"
    "mask.b.srli (t2,t1),(t2,t1),31\n"
    
    "mask.b.and (t2,t1), (t2,t1), (a1,a0)\n"
    "mask.b.sub (t2,t1), (t2,t1), (a1,a0)\n"
    "mask.b.not (t2,t1), (t2,t1)\n"
    
    "mask.b.xor (t4,t3),(a7,a6),(a1,a0)\n"
    "mask.b.and (a7,a6),(t4,t3),(t2,t1)\n"
    "mask.b.xor (a1,a0),(a1,a0),(a7,a6)\n"

    "mask.b.xor (a3,a2),(a3,a2),(a1,a0)\n"
    "mask.b.xor (a5,a4),(a5,a4),(a1,a0)\n"

    "sw a2, (%[a0])\n" //a3,a2 = a
    "sw a3, (%[a1])\n"
    "sw a4, (%[b0])\n" //a4,a5 = b
    "sw a5, (%[b1])\n"
    :
    : [a0] "r" (a), [a1] "r" (a+a_stride), [b0] "r" (b), [b1] "r" (b+b_stride)
    : "t4", "t3", "t2", "t1","a7", "a6", "a5", "a4", "a3", "a2", "a1", "a0", "memory");
#elif defined USE_MASKED_EXT && defined USE_MASKED_ISA
    asm volatile (
    "li a0, 1\n"
    "mask.b.mask (a1,a0), a0\n" //a1,a0 = 1
    "lw a2, (%[a0])\n" //a3,a2 = a
    "lw a3, (%[a1])\n"
    "lw a4, (%[b0])\n" //a4,a5 = b
    "lw a5, (%[b1])\n"
    "mask.b.xor (a7,a6), (a3,a2), (a5,a4)\n" //a7,a6 = ab
    "mask.b.srli (t6,t5),(a3,a2),1\n" //t6,t5 shifted a
    "mask.b.srli (t4,t3),(a5,a4),1\n" //t4,t3 shifted b
    "mask.b.cmpgt (t2,t1),(t4,t3),(t6,t5)\n" //t2,t1 = res
    "mask.b.cmov (a1,a0),(a7,a6),(t2,t1)\n" //cmov tmp
    "mask.b.xor (a3,a2),(a3,a2),(a1,a0)\n"
    "mask.b.xor (a5,a4),(a5,a4),(a1,a0)\n"

    "sw a2, (%[a0])\n" //a3,a2 = a
    "sw a3, (%[a1])\n"
    "sw a4, (%[b0])\n" //a4,a5 = b
    "sw a5, (%[b1])\n"
    :
    : [a0] "r" (a), [a1] "r" (a+a_stride), [b0] "r" (b), [b1] "r" (b+b_stride)
    : "t6", "t5", "t4", "t3", "t2", "t1","a7", "a6", "a5", "a4", "a3", "a2", "a1", "a0", "memory");
#else
    b_mask(tmp, 1, 0);
    b_xor(ab, 1, a, a_stride, b, b_stride);
    // the lower two bits in a and b are payload, by shifting to the right
    // we reach a 31-bit operand and can use the faster cmpg implementation with sub
    b_srl1(a_shifted, 1, a, a_stride);
    b_srl1(b_shifted, 1, b, b_stride);
    b_cmpg(res, 1, b_shifted, 1,  a_shifted, 1);
    b_cmov(tmp, 1, ab, 1, res, 1);
    // if a > b : tmp = a^b, else tmp = 0
    b_xor(a, a_stride, a, a_stride, tmp, 1);
    b_xor(b, b_stride, b, b_stride, tmp, 1);
#endif
}

static void masked_crypto_sort(uint32_t x[NSHARES][N]) {
    long long top, p, q, r, i;

    if (N < 2)
        return;

    top = 1;
    while (top < N - top)
        top += top;

    for (p = top; p > 0; p >>= 1) {
        for (i = 0; i < N - p; ++i) {
            if (!(i & p))
                masked_minmax(&x[0][i], N,  &x[0][i + p], N);
        }
        i = 0;
        for (q = top; q > p; q >>= 1) {
            for (; i < N - q; ++i) {
                if (!(i & p)) {
                    for (r = q; r > p; r >>= 1)
                        masked_minmax(&x[0][i + p], N,  &x[0][i + r], N);                  
                }
            }
        }
    }
}

// expects randomness in the upper 30bits for all entries in x
void sample_sort(uint32_t x[NSHARES][N]) {
    int i;
    uint32_t one[NSHARES];
    uint32_t bitmask[NSHARES];

    b_mask(one, 1, 1);
    b_mask(bitmask, 1, 0xFFFFFFFC);

    for (i = 0; i < W; i++) {
        // set first W coeffs to 1:
        b_or(&x[0][i], N, &x[0][i], N, one, 1);
    }
    for (; i < N; i++) {
        // set remaining coeffs to 0:
        b_and(&x[0][i], N, &x[0][i], N, bitmask, 1);
    }

    masked_crypto_sort(x);
}