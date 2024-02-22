#ifndef MASKED_POLY_FROMMSG__H
#define MASKED_POLY_FROMMSG__H

#include <stdint.h>
#include "masked.h"
#include "core_ops.h"

#define KYBER_Q 3329
#define KYBER_N 256
#define KYBER_SYMBYTES 32

void masked_poly_frommsg(uint32_t r[NSHARES][KYBER_N], /*const*/ uint8_t msg[NSHARES][KYBER_SYMBYTES]);

#endif