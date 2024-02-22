/* Copyright 2022 UCLouvain, Belgium and PQM4 contributors
 *
 * This file is part of pqm4_masked.
 *
 * pqm4_masked is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, version 3.
 *
 * pqm4_masked is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * pqm4_masked. If not, see <https://www.gnu.org/licenses/>.
 *
 * Edited by (anonymized)
 */

#ifndef MASKED_H
#define MASKED_H
#include "randombytes.h"
#include <stdint.h>

#define NSHARES 2
#define BSSIZE 32

#if defined(SET2)
#define N 953
#define W 396
#define LOG_N 10
#define LOG_W 9
#define CMP_REP 45.20
#define CMP_THR 13
#define FRODO_N 976
#define BIKE_N 24659
#define KECCAK_INLEN 73 // 2 keccak permutations in masked_sha3_512
#define HAWK_N 512
#elif defined(SET3)
#define N 6688
#define W 128
#define LOG_N 13
#define LOG_W 7
#define CMP_REP 28.88
#define CMP_THR 5
#define FRODO_N 1344
#define BIKE_N 40973
#define KECCAK_INLEN 145 // 3 keccak permutations in masked_sha3_512
#define HAWK_N 1024
#elif defined(SET4)
#define N 12323
#define W 71
#define LOG_N 14
#define LOG_W 7
#define CMP_REP 21.30
#define CMP_THR 3
#define FRODO_N 640
#define BIKE_N 12323
#define KECCAK_INLEN 217 // 4 keccak permutations in masked_sha3_512
#define HAWK_N 256
#elif defined(SET5)
#define N 49318
#define W 199
#define LOG_N 16
#define LOG_W 8
#define CMP_REP 1
#define CMP_THR 1
#define FRODO_N 640
#define BIKE_N 12323
#define KECCAK_INLEN 289 // 5 keccak permutations in masked_sha3_512
#define HAWK_N 256
#elif defined(SET6)
#define N 35851
#define W 114
#define LOG_N 16
#define LOG_W 7
#define CMP_REP 1
#define CMP_THR 1
#define FRODO_N 640
#define BIKE_N 12323
#define KECCAK_INLEN 361 // 5 keccak permutations in masked_sha3_512
#define HAWK_N 256
#else // SET 1 is the default
#define N 677
#define W 254
#define LOG_N 10
#define LOG_W 8
#define CMP_REP 31.59
#define CMP_THR 3
#define FRODO_N 640
#define BIKE_N 12323
#define KECCAK_INLEN 16 // 1 keccak permutation in masked_sha3_512
#define HAWK_N 256
#endif // SET


#define PAD32(X) ((((X) + 31) / 32) * 32)
#define N_PADDED PAD32(N)
#define W_PADDED PAD32(W)
#define B_PADDED (W_PADDED * 10) 

#define PAD64(X) ((((X) + 63) / 64) * 64)

union U64 {
    uint64_t v64;
    uint32_t v32[2];
};

#endif // MASKED_H
