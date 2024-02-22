#ifndef BIKE_ROTR_H
#define BIKE_ROTR_H

#include "masked.h"
#include <stddef.h>
#include <stdint.h>

#define R_QWORDS ((BIKE_N + 31) / 32)

// Copied from (Kaz answer)
// https://stackoverflow.com/questions/466204/rounding-up-to-next-power-of-2
#define UPTOPOW2_0(v) ((v)-1)
#define UPTOPOW2_1(v) (UPTOPOW2_0(v) | (UPTOPOW2_0(v) >> 1))
#define UPTOPOW2_2(v) (UPTOPOW2_1(v) | (UPTOPOW2_1(v) >> 2))
#define UPTOPOW2_3(v) (UPTOPOW2_2(v) | (UPTOPOW2_2(v) >> 4))
#define UPTOPOW2_4(v) (UPTOPOW2_3(v) | (UPTOPOW2_3(v) >> 8))
#define UPTOPOW2_5(v) (UPTOPOW2_4(v) | (UPTOPOW2_4(v) >> 16))

#define UPTOPOW2(v) (UPTOPOW2_5(v) + 1)

#define R_QWORDS_HALF_LOG2 UPTOPOW2(R_QWORDS / 2)

void rotr_big(uint32_t in[NSHARES][3 * R_QWORDS], uint32_t qw_num[NSHARES]);

#endif // BIKE_ROTR_H