#ifndef CORE_OPS_H
#define CORE_OPS_H

#include <stddef.h>
#include <stdint.h>

#include "gadgets.h"
#include "masked.h"

#ifdef USE_MASKED_ISA

static inline void b_mask(uint32_t *a, size_t a_stride, uint32_t b) {

    register uint32_t *loca0 asm("t6") = 0;
    register uint32_t *loca1 asm("t5") = 0;
    
    uint32_t result;
    asm volatile (
        "mask.b.mask (t6,t5), %[input]\n"
        : "=r" (loca0), "=r" (loca1)
        : [input] "r" (b)
        :
    );

    *a = loca0;
    *(a+a_stride) = loca1;
}
static inline void b_unmask(uint32_t *z, const uint32_t *a, size_t a_stride) {

    register uint32_t *loca0 asm("t1") = *a;
    register uint32_t *loca1 asm("t2") = *(a+a_stride);
    
    uint32_t result;
    asm volatile (
        "mask.b.unmask %[result],(t2,t1)\n"
        : [result] "=r" (result)
        : "r" (loca0), "r" (loca1)
        :
    );

    *z = result;
}
static inline void b_not(uint32_t *z, size_t z_stride, const uint32_t *a,
                  size_t a_stride) {

    register uint32_t *loca0 asm("t1") = *a;
    register uint32_t *loca1 asm("t2") = *(a+a_stride);
    
    volatile register uint32_t *locz0 asm("t4") = 0;
    volatile register uint32_t *locz1 asm("t3") = 0;
    
    asm volatile (
        "mask.b.not (t4,t3),(t2,t1)\n"
        : "=r" (locz0), "=r" (locz1)
        : "r" (loca0), "r" (loca1)
        :
    );
    *z = locz0;
    *(z+z_stride) = locz1;
}
static inline void b_and(uint32_t *z, size_t z_stride, const uint32_t *a,
                  size_t a_stride, const uint32_t *b, size_t b_stride) {
    
    register uint32_t *loca0 asm("t1") = *a;
    register uint32_t *loca1 asm("t2") = *(a+a_stride);

    register uint32_t *locb0 asm("t3") = *b;
    register uint32_t *locb1 asm("t4") = *(b+b_stride);
    
    volatile register uint32_t *locz0 asm("t5") = 0;
    volatile register uint32_t *locz1 asm("t6") = 0;
    
    asm volatile (
        "mask.b.and (t6,t5),(t4,t3),(t2,t1)\n"
        : "=r" (locz0), "=r" (locz1)
        : "r" (loca0), "r" (loca1),"r" (locb0),"r" (locb1)
        :
    );
    *z = locz0;
    *(z+z_stride) = locz1;
}
static inline void b_or(uint32_t *z, size_t z_stride, const uint32_t *a,
                 size_t a_stride, const uint32_t *b, size_t b_stride) {
                
    register uint32_t *loca0 asm("t1") = *a;
    register uint32_t *loca1 asm("t2") = *(a+a_stride);

    register uint32_t *locb0 asm("t3") = *b;
    register uint32_t *locb1 asm("t4") = *(b+b_stride);
    
    volatile register uint32_t *locz0 asm("t5") = 0;
    volatile register uint32_t *locz1 asm("t6") = 0;
    
    asm volatile (
        "mask.b.ior (t6,t5),(t4,t3),(t2,t1)\n"
        : "=r" (locz0), "=r" (locz1)
        : "r" (loca0), "r" (loca1),"r" (locb0),"r" (locb1)
        :
    );
    *z = locz0;
    *(z+z_stride) = locz1;
}
static inline void b_xor(uint32_t *z, size_t z_stride, const uint32_t *a,
                  size_t a_stride, const uint32_t *b, size_t b_stride) {
    
    register uint32_t *loca0 asm("t1") = *a;
    register uint32_t *loca1 asm("t2") = *(a+a_stride);

    register uint32_t *locb0 asm("t3") = *b;
    register uint32_t *locb1 asm("t4") = *(b+b_stride);
    
    volatile register uint32_t *locz0 asm("t5") = 0;
    volatile register uint32_t *locz1 asm("t6") = 0;
    
    asm volatile (
        "mask.b.xor (t6,t5),(t4,t3),(t2,t1)\n"
        : "=r" (locz0), "=r" (locz1)
        : "r" (loca0), "r" (loca1),"r" (locb0),"r" (locb1)
        :
    );
    *z = locz0;
    *(z+z_stride) = locz1;       
}
static inline void b_slli(uint32_t *z, size_t z_stride, const uint32_t *a,
                   size_t a_stride, size_t shamt) {
    size_t i;
    for (i = 0; i < NSHARES; i++) {
        z[i * z_stride] = a[i * a_stride] << shamt;
    }              
}
static inline void b_srli(uint32_t *z, size_t z_stride, const uint32_t *a,
                   size_t a_stride, size_t shamt) {
    size_t i;
    for (i = 0; i < NSHARES; i++) {
        z[i * z_stride] = a[i * a_stride] >> shamt;
    }
}

static inline void b_sll1(uint32_t *z, size_t z_stride, const uint32_t *a,
                   size_t a_stride) {

    register uint32_t *loca0 asm("t1") = *a;
    register uint32_t *loca1 asm("t2") = *(a+a_stride);
    
    volatile register uint32_t *locz0 asm("t5") = 0;
    volatile register uint32_t *locz1 asm("t6") = 0;
    
    asm volatile (
        "mask.b.slli (t6,t5),(t2,t1),1\n"
        : "=r" (locz0), "=r" (locz1)
        : "r" (loca0), "r" (loca1)
        :
    );
    *z = locz0;
    *(z+z_stride) = locz1;             
}

static inline void b_srl1(uint32_t *z, size_t z_stride, const uint32_t *a,
                   size_t a_stride) {

    register uint32_t *loca0 asm("t1") = *a;
    register uint32_t *loca1 asm("t2") = *(a+a_stride);
    
    volatile register uint32_t *locz0 asm("t5") = 0;
    volatile register uint32_t *locz1 asm("t6") = 0;
    
    asm volatile (
        "mask.b.srli (t6,t5),(t2,t1),1\n"
        : "=r" (locz0), "=r" (locz1)
        : "r" (loca0), "r" (loca1)
        :
    );
    *z = locz0;
    *(z+z_stride) = locz1;    
}

static inline void b_add(uint32_t *z, size_t stride_z, const uint32_t *a,
                   size_t stride_a, const uint32_t *b, size_t stride_b) {

    register uint32_t *loca0 asm("t1") = *a;
    register uint32_t *loca1 asm("t2") = *(a+stride_a);

    register uint32_t *locb0 asm("t3") = *b;
    register uint32_t *locb1 asm("t4") = *(b+stride_b);
    
    volatile register uint32_t *locz0 asm("t5") = 0;
    volatile register uint32_t *locz1 asm("t6") = 0;
    
    asm volatile (
        "mask.b.add (t6,t5),(t2,t1),(t4,t3)\n"
        : "=r" (locz0), "=r" (locz1)
        : "r" (loca0), "r" (loca1),"r" (locb0),"r" (locb1)
        :
    );
    *z = locz0;
    *(z+stride_z) = locz1;
}
static inline void b_sub(uint32_t *z, size_t stride_z, const uint32_t *a,
                   size_t stride_a, const uint32_t *b, size_t stride_b) {

    register uint32_t *loca0 asm("t1") = *a;
    register uint32_t *loca1 asm("t2") = *(a+stride_a);

    register uint32_t *locb0 asm("t3") = *b;
    register uint32_t *locb1 asm("t4") = *(b+stride_b);
    
    volatile register uint32_t *locz0 asm("t5") = 0;
    volatile register uint32_t *locz1 asm("t6") = 0;
    
    asm volatile (
        "mask.b.sub (t6,t5),(t2,t1),(t4,t3)\n"
        : "=r" (locz0), "=r" (locz1)
        : "r" (loca0), "r" (loca1),"r" (locb0),"r" (locb1)
        :
    );
    *z = locz0;
    *(z+stride_z) = locz1;
    
}

#else  // DONT USE_MASKED_ISA

static inline void b_mask(uint32_t *a, size_t a_stride, uint32_t b) {
    size_t i;
    a[0] = b;
    for (i = 1; i < NSHARES; i++) {
        uint32_t r = randomint();
        a[i * a_stride] = r;
        a[0] ^= r;
    }
}

static inline void b_unmask(uint32_t *z, const uint32_t *a, size_t a_stride) {
    size_t i;

    *z = a[0 * a_stride];
    for (i = 1; i < NSHARES; i++) {
        *z ^= a[i * a_stride];
    }
}

static inline void b_not(uint32_t *z, size_t z_stride, const uint32_t *a,
                  size_t a_stride) {
    z[0] = ~a[0];
    if (z != a) // this will be left out due to inline
    {
        size_t i;
        for (i = 1; i < NSHARES; i++) {
            z[i * z_stride] = a[i * a_stride];
        }
    }
}

static inline void b_and(uint32_t *z, size_t z_stride, const uint32_t *a,
                  size_t a_stride, const uint32_t *b, size_t b_stride) {
    masked_and(NSHARES, z, z_stride, a, a_stride, b, b_stride);
}

static inline void b_or(uint32_t *z, size_t z_stride, const uint32_t *a,
                 size_t a_stride, const uint32_t *b, size_t b_stride) {
    masked_or(NSHARES, z, z_stride, a, a_stride, b, b_stride);
}

static inline void b_xor(uint32_t *z, size_t z_stride, const uint32_t *a,
                  size_t a_stride, const uint32_t *b, size_t b_stride) {
    masked_xor(NSHARES, z, z_stride, a, a_stride, b, b_stride);
}

static inline void b_slli(uint32_t *z, size_t z_stride, const uint32_t *a,
                   size_t a_stride, size_t shamt) {
    size_t i;
    for (i = 0; i < NSHARES; i++) {
        z[i * z_stride] = a[i * a_stride] << shamt;
    }
}

static inline void b_srli(uint32_t *z, size_t z_stride, const uint32_t *a,
                   size_t a_stride, size_t shamt) {
    size_t i;
    for (i = 0; i < NSHARES; i++) {
        z[i * z_stride] = a[i * a_stride] >> shamt;
    }
}

static inline void b_sll1(uint32_t *z, size_t z_stride, const uint32_t *a,
                   size_t a_stride) {
    size_t i;
    for (i = 0; i < NSHARES; i++) {
        z[i * z_stride] = a[i * a_stride] << 1;
    }
}

static inline void b_srl1(uint32_t *z, size_t z_stride, const uint32_t *a,
                   size_t a_stride) {
    size_t i;
    for (i = 0; i < NSHARES; i++) {
        z[i * z_stride] = a[i * a_stride] >> 1;
    }
}

static inline void b_add(uint32_t *z, size_t stride_z, const uint32_t *a,
                   size_t stride_a, const uint32_t *b, size_t stride_b) {
  add(z, stride_z, a, stride_a, b, stride_b, 32);
}

// UNSAFE, only for functional testing!
static inline void b_sub(uint32_t *z, size_t stride_z, const uint32_t *a,
                   size_t stride_a, const uint32_t *b, size_t stride_b) {
    uint32_t A,B;
    b_unmask(&A, a, stride_a);
    b_unmask(&B, b, stride_b);
    b_mask(z, stride_z, A-B);
}
#endif // USE_MASKED_ISA


// Our instruction extension:

static inline void b_cmov(uint32_t *z, size_t z_stride, const uint32_t *a,
                   size_t a_stride, const uint32_t *cond, size_t cond_stride) {
    // expects cond to be one bit
#ifdef USE_MASKED_EXT
    register uint32_t *loccond0 asm("t1") = *cond;
    register uint32_t *loccond1 asm("t2") = *(cond+cond_stride);

    register uint32_t *locn0 asm("t3") = *a;
    register uint32_t *locn1 asm("t4") = *(a+a_stride);
    
    volatile register uint32_t *locz0 asm("t5") = *z;
    volatile register uint32_t *locz1 asm("t6") = *(z+z_stride);
    
    asm volatile (
        "mask.b.cmov (t6,t5),(t4,t3),(t2,t1)\n"
        : "=r" (locz0), "=r" (locz1)
        : "r" (loccond0), "r" (loccond1),"r" (locn0),"r" (locn1)
        :
    );
    *z = locz0;
    *(z+z_stride) = locz1;

#elif defined USE_MASKED_ISA
    uint32_t cond_mask[NSHARES];
    uint32_t xor [NSHARES];

    for (int s = 0; s < NSHARES; s++) {
        // expand from 1 bit to registerwidth:
        cond_mask[s] = ~((cond[s * cond_stride] & 1) - 1);
    }
    // xor = old ^ new -> if cond : old = old ^ old ^ new = new, else : old =
    // old ^ 0 = old

    b_xor(xor, 1, z, z_stride, a, a_stride);
    b_and(xor, 1, xor, 1, cond_mask, 1);
    b_xor(z, z_stride, z, z_stride, xor, 1);    
#else
    masked_sel(cond, cond_stride, z, z_stride, a, a_stride);
#endif
}

// returns 1 if b > a
static inline void b_cmpg(uint32_t *z, size_t stride_z, const uint32_t *a,
                   size_t stride_a, const uint32_t *b, size_t stride_b) {
#ifdef USE_MASKED_EXT

    register uint32_t *loca0 asm("t1") = *a;
    register uint32_t *loca1 asm("t2") = *(a+stride_a);

    register uint32_t *locb0 asm("t3") = *b;
    register uint32_t *locb1 asm("t4") = *(b+stride_b);
    
    volatile register uint32_t *locz0 asm("t5") = 0;
    volatile register uint32_t *locz1 asm("t6") = 0;
    
    asm volatile (
        "mask.b.cmpgt (t6,t5),(t2,t1),(t4,t3)\n"
        : "=r" (locz0), "=r" (locz1)
        : "r" (loca0), "r" (loca1),"r" (locb0),"r" (locb1)
        :
    );
    *z = locz0;
    *(z+stride_z) = locz1;

#elif defined USE_MASKED_ISA
    b_sub(z, stride_z, a, stride_a, b, stride_b);
    b_srli(z, stride_z, z, stride_z, 31);
#else
    cmpg(z, stride_z, a, stride_a, b, stride_b, 32);
#endif
}

static inline void b_cmpeq(uint32_t *z, size_t stride_z, const uint32_t *a,
                    size_t stride_a, const uint32_t *b, size_t stride_b) {
#ifdef USE_MASKED_EXT

    register uint32_t *loca0 asm("t1") = *a;
    register uint32_t *loca1 asm("t2") = *(a+stride_a);

    register uint32_t *locb0 asm("t3") = *b;
    register uint32_t *locb1 asm("t4") = *(b+stride_b);
    
    volatile register uint32_t *locz0 asm("t5") = 0;
    volatile register uint32_t *locz1 asm("t6") = 0;
    
    asm volatile (
        "mask.b.cmpeq (t6,t5),(t2,t1),(t4,t3)\n"
        : "=r" (locz0), "=r" (locz1)
        : "r" (loca0), "r" (loca1),"r" (locb0),"r" (locb1)
        :
    );
    *z = locz0;
    *(z+stride_z) = locz1;

#elif defined USE_MASKED_ISA
    // for the not extended masked ISA, we can compute the cmpeq
    // more efficiently with subtractions:
    uint32_t tmp1[NSHARES];
    uint32_t tmp2[NSHARES];
    b_sub(tmp1, 1, a, stride_a, b, stride_b);
    b_srli(tmp1, 1, tmp1, 1, 31);
    b_sub(tmp2, 1, b, stride_b, a, stride_a);
    b_srli(tmp2, 1, tmp2, 1, 31);
    b_or(tmp1, 1, tmp1, 1, tmp2, 1);
    b_not(tmp1, 1, tmp1, 1);
    b_mask(tmp2, 1, 1);
    b_and(z, stride_z, tmp1, 1, tmp2, 1);
#else
    masked_eq(z, stride_z, a, stride_a, b, stride_b, 32);
#endif 
}

#endif // CORE_OPS_H

