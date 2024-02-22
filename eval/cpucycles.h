#ifndef CPUCYCLES_H
#define CPUCYCLES_H

#include <stdint.h>

static inline uint64_t cpucycles(void) {
    uint32_t low;
    uint32_t high;

    asm volatile("rdcycleh %[x] \n"
    "rdcycle %[y]\n" 
    : [x] "=r" (high), [y] "=r" (low) 
    : 
    : "memory");

    uint64_t result = high;
    result = result << 32;
    result = result | low;

    /*
    __asm__ volatile("rdtsc; shlq $32,%%rdx; orq %%rdx,%%rax"
                     : "=a"(result)
                     :
                     : "%rdx");
    */
    return result;
}

static uint64_t cpucycles_overhead(void) {
    uint64_t t0, t1, overhead = -1LL;
    unsigned int i;
    
    for (i = 0; i < 100000; i++) {
        t0 = cpucycles();
        __asm__ volatile("");
        t1 = cpucycles();
        if (t1 - t0 < overhead)
            overhead = t1 - t0;
    }

    return overhead;
}

#endif // CPUCYCLES_H