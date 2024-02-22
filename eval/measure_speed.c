#include "bike_rotr.h"
#include "core_ops.h"
#include "cpucycles.h"
#include "frodo_sampler.h"
#include "haetae_sampler.h"
#include "hawk_sampler.h"
#include "masked.h"
#include "masked_fwsampling.h"
#include "masked_keccak.h"
#include "masked_poly_frommsg.h"
#include <math.h>
#include <stdint.h>
//#include <stdio.h>
#include <string.h>

#define NTESTS 1

void send(char byte) {
    volatile uint32_t* uart_tx = 0x40001004;
    *uart_tx = byte & 0xFF;
}

void send_str(char* str) {
    uint32_t i = 0;
    while(str[i] != 0) {
        send(str[i]);
        i++;
        for (int i = 0; i < 100; i++) { 
                __asm__("nop");
        }
    }
}


uint64_t cycles[NTESTS];

static uint64_t average(uint64_t *t, size_t tlen) {
  size_t i;
  uint64_t acc=0;

  for(i=0;i<tlen;i++)
    acc += t[i];

  return acc/tlen;
}

static void print_results(const char *s, uint64_t *t, size_t tlen) {
  static uint64_t overhead = -1;
  if(overhead  == (uint64_t)-1)
    overhead = cpucycles_overhead();
   
   uint64_t result = (unsigned long long)average(t, tlen) - overhead;
   result /= 1000; // we only report kcycle numbers
   send_str(s);
   send('\n');
   for (int i = 0; i < 100; i++) { 
                __asm__("nop");
    }
    for(int k = 0; k < 8; k++) {
        send((result >> k*8) & 0xFF);
        for (int i = 0; i < 100; i++) { 
                __asm__("nop");
        }
    }
    send('\n');
    for (int i = 0; i < 100; i++) { 
                __asm__("nop");
    }
  //printf("%s ", s);
  //printf("avg: %llu cycles\n", (unsigned long long)average(t, tlen) - overhead);
}

static void measure_bs_sample_reject() {
    uint32_t p_bs[2][W_PADDED];
    uint32_t rand[2][B_PADDED];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {
        randombytes((uint8_t *)&rand[0][0], sizeof(uint32_t) * 2 * B_PADDED);
        cycle_start = cpucycles();
        bs_sample_reject(p_bs, rand);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_sample_reject() {
    uint32_t p[2][W];
    uint32_t rand[2][B_PADDED];
    uint64_t cycle_start;
    //uint32_t ret;

    for(int i = 0; i < NTESTS; ++i) {
        randombytes((uint8_t *)&rand[0][0], sizeof(uint32_t) * 2 * B_PADDED);
        cycle_start = cpucycles();
        sample_reject(p, rand);
        //ret = sample_reject(p, rand);
        cycles[i] = cpucycles() - cycle_start;
        /*
        for(int k = 0; k < 4; k++) {
            send((ret >> k*8) & 0xFF);
            for (int i = 0; i < 100; i++) { 
                    __asm__("nop");
            }
        }
        send('\n');
        for (int i = 0; i < 100; i++) { 
                    __asm__("nop");
        }
        */
    }
}

static void measure_bs_sample_fisher_yates() {
    uint32_t p[2][W_PADDED];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {    
        randombytes((uint8_t *)&p[0][0], sizeof(uint32_t) * 2 * W_PADDED);
        cycle_start = cpucycles();
        bs_sample_fisher_yates(p);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_bs_i2c() {
    uint32_t p[2][W_PADDED], pp[2][N_PADDED/32];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {    
        randombytes((uint8_t *)&p[0][0], sizeof(uint32_t) * 2 * W_PADDED);
        bs_sample_fisher_yates(p);
        cycle_start = cpucycles();
        bs_coeff_to_regular(pp, p);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_sample_fisher_yates() {
    uint32_t p[2][W_PADDED];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {  
        randombytes((uint8_t *)&p[0][0], sizeof(uint32_t) * 2 * W_PADDED);
        cycle_start = cpucycles();
        sample_fisher_yates(p);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_i2c() {
    uint32_t p[2][W_PADDED], pp[2][N_PADDED/32];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {    
        randombytes((uint8_t *)&p[0][0], sizeof(uint32_t) * 2 * W_PADDED);
        bs_sample_fisher_yates(p);
        cycle_start = cpucycles();
        coeff_to_regular(pp, p);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_bs_sample_sort() {
    uint32_t p_bs[2][PAD64(N)];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {     
        randombytes((uint8_t *)&p_bs[0][0], sizeof(uint32_t) * 2 * N);
        cycle_start = cpucycles();
        bs_sample_sort(p_bs);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_sample_sort() {
    uint32_t p[2][N];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {    
        randombytes((uint8_t *)&p[0][0], sizeof(uint32_t) * 2 * N);
        cycle_start = cpucycles();
        sample_sort(p);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_bs_sample_repand() {
    uint32_t p[2][(N + 31) / 32];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {       
        memset(&p[0][0], 0, sizeof(uint32_t) * 2 * ((N + 31) / 32));
        cycle_start = cpucycles();
        bs_sample_repand(p);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_sample_cmp() {
    uint32_t randomness[NSHARES][N];
    uint32_t p[NSHARES][(N + 31) / 32];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {   
        cycle_start = cpucycles(); 
        // runtime of one iteration is indep of randomness
        // avg number of iterations is precomputed and must be factored in
        sample_cmp(p, randomness, CMP_THR);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_bs_sample_cmp() {
    uint32_t randomness[PAD32(N)/32][NSHARES][L];
    uint32_t p[NSHARES][(N + 31) / 32];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {   
        cycle_start = cpucycles(); 
        // runtime of one iteration is indep of randomness
        // avg number of iterations is precomputed and must be factored in
        bs_sample_cmp(p, randomness);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_bike_rotr() {
    uint32_t in[NSHARES][3 * R_QWORDS];
    uint32_t out[NSHARES][3 * R_QWORDS];
    uint32_t qw_num[NSHARES];
    uint32_t rot;
    uint64_t cycle_start;

    randombytes((uint8_t *)&in[0][0], sizeof(uint32_t) * R_QWORDS);
    randombytes((uint8_t *)&in[1][0], sizeof(uint32_t) * R_QWORDS);
    memcpy(&in[0][1 * R_QWORDS], &in[0][0], sizeof(uint32_t) * R_QWORDS);
    memcpy(&in[0][2 * R_QWORDS], &in[0][0], sizeof(uint32_t) * R_QWORDS);
    memcpy(&in[1][1 * R_QWORDS], &in[1][0], sizeof(uint32_t) * R_QWORDS);
    memcpy(&in[1][2 * R_QWORDS], &in[1][0], sizeof(uint32_t) * R_QWORDS);

    memcpy(&out[0][0], &in[0][0], sizeof(uint32_t) * 6 * R_QWORDS);

    rot = randomint();
    rot &= (UPTOPOW2(R_QWORDS) - 1); // log2(rot) <= log2(R_QWORDS)
    rot = 1;
    b_mask(qw_num, 1, rot);

    for(int i = 0; i < NTESTS; ++i) { 
        cycle_start = cpucycles();
        rotr_big(out, qw_num);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_frodo_sampler() {
    uint16_t samples[NSHARES * FRODO_N];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {    
        randombytes((uint8_t *)samples, sizeof(uint16_t) * NSHARES * FRODO_N);
        cycle_start = cpucycles();
        frodo_sampler(samples);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_bs_frodo_sampler() {
    uint32_t samples[NSHARES][(FRODO_N+31)/32][16];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {     
        randombytes((uint8_t *)samples, sizeof(uint32_t) * NSHARES * ((FRODO_N+31)/32) * 16);
        cycle_start = cpucycles();
        frodo_sampler_bs(samples);
        cycles[i] = cpucycles() - cycle_start;
    }
}

int __attribute__((naked)) entry() {
    __asm__("j _start");
}
static void measure_haetae_sampler() {
    uint16_t samples[NSHARES * HAETAE_N];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {     
        randombytes((uint8_t *)samples, sizeof(uint16_t) * NSHARES * HAETAE_N);
        cycle_start = cpucycles();
        haetae_sampler(samples);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_bs_haetae_sampler() {
    uint32_t samples[NSHARES][(HAETAE_N+31)/32][16];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {     
        randombytes((uint8_t *)samples, sizeof(uint32_t) * NSHARES * ((HAETAE_N+31)/32) * 16);
        cycle_start = cpucycles();
        haetae_sampler_bs(samples);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_bs_hawk_sampler() {
    uint32_t samples[NSHARES][(HAWK_N+31)/32][80];
    uint64_t cycle_start;

    for(int i = 0; i < NTESTS; ++i) {     
        randombytes((uint8_t *)samples, sizeof(uint32_t) * NSHARES * ((HAWK_N+31)/32) * 80);
        cycle_start = cpucycles();
        hawk_sampler_bs(samples);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_keccak() {
    uint32_t in[NSHARES][(KECCAK_INLEN+3)/4], out[NSHARES][16];
    uint64_t cycle_start;

    randombytes((uint8_t *)in, sizeof(uint32_t) * NSHARES * ((KECCAK_INLEN+3)/4));
    for(int i = 0; i < NTESTS; ++i) {     
        cycle_start = cpucycles();
        masked_sha3_512((uint8_t*)out, 64, 1, (uint8_t*)in, KECCAK_INLEN, ((KECCAK_INLEN+3)/4)*4, 1);
        cycles[i] = cpucycles() - cycle_start;
    }
}

static void measure_poly_frommsg() {
    uint32_t r[NSHARES][KYBER_N];
    uint8_t msg[NSHARES][KYBER_SYMBYTES];
    uint64_t cycle_start;

    randombytes(msg[0], NSHARES * KYBER_SYMBYTES);

    for (size_t i = 0; i < NTESTS; i++)
    {  
      cycle_start = cpucycles();
      masked_poly_frommsg(r, msg);
      cycles[i] = cpucycles() - cycle_start;
    }
}

int main() {
    measure_keccak();
    print_results("keccak", cycles, NTESTS);
    measure_sample_fisher_yates();
    print_results("fy", cycles, NTESTS);
    measure_bs_sample_fisher_yates();
    print_results("bs fy", cycles, NTESTS);
    measure_i2c();
    print_results("i2c", cycles, NTESTS);
    measure_bs_i2c();
    print_results("bs i2c", cycles, NTESTS);
    measure_sample_sort();
    print_results("sort", cycles, NTESTS);
    measure_bs_sample_sort();
    print_results("bs sort", cycles, NTESTS);
    // Only relevant for SET1 - SET4:
#if defined(SET1) || defined(SET2) || defined(SET3) || defined(SET4)
    measure_sample_reject();
    print_results("reject", cycles, NTESTS);
    measure_bs_sample_reject();
    print_results("bs reject", cycles, NTESTS);
    measure_bs_sample_repand();
    print_results("bs repAND", cycles, NTESTS);
    measure_sample_cmp();
    print_results("cmp", cycles, NTESTS);
    measure_bs_sample_cmp();
    print_results("bs cmp", cycles, NTESTS);
#endif
    // Only relevant for SET1 - SET3:
#if defined(SET1) || defined(SET2) || defined(SET3)
    measure_bike_rotr();   
    print_results("bike rot", cycles, NTESTS); 
    measure_frodo_sampler();
    print_results("frodo sample", cycles, NTESTS);
    measure_bs_frodo_sampler();
    print_results("bs frodo sample", cycles, NTESTS);
    measure_bs_hawk_sampler();
    print_results("bs hawk sample", cycles, NTESTS);

#endif

#if defined(SET1)
    measure_haetae_sampler();
    print_results("haetae sample", cycles, NTESTS);
    measure_bs_haetae_sampler();
    print_results("bs haetae sample", cycles, NTESTS);
    measure_poly_frommsg();
    print_results("kyber poly_frommsg", cycles, NTESTS);
    
#endif
    while(1) {
        __asm__ volatile ("nop");
    }
}
