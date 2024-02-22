#ifndef MASKED_CMP__H
#define MASKED_CMP__H

#if N == 12323   // BIKE-1
#define L 9
#elif N == 24659 // BIKE-3
#define L 8
#elif N == 40973 // BIKE-5
#define L 11
#elif N == 17669 // HQC-128
#define L 8
#elif N == 35851 // HQC-192
#define L 10
#elif N == 57637 // HQC-256
#define L 12
#elif N == 3488  // McE 348864
#define L 8
#elif N == 4608  // McE 460896
#define L 9
#elif N == 6688  // McE 6688128
#define L 8
#elif N == 6960  // McE 6960119
#define L 9
#elif N == 8192  // McE 8192128
#define L 6
#elif N == 509
#define L 1
#elif N == 677
#define L 3
#elif N == 821
#define L 3
#elif N == 653
#define L 4
#elif N == 761
#define L 3
#elif N == 857
#define L 3
#elif N == 953
#define L 5
#elif N == 1013
#define L 4
#elif N == 1277
#define L 3
#else
#define L 1 //default
#endif

#ifdef NTRULPR
#if N == 653
#define L 3
#elif N == 761
#define L 6
#elif N == 857
#define L 3
#elif N == 953
#define L 3
#elif N == 1013
#define L 3
#elif N == 1277
#define L 5
#else
#error
#endif
#endif


#endif
