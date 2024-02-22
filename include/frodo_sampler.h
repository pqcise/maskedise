#ifndef FRODO_SAMPLER_H
#define FRODO_SAMPLER_H

#include <stdint.h>
#include <stdio.h>
#include "masked.h"

void frodo_sampler(uint16_t samples[NSHARES * FRODO_N]);
void frodo_sampler_bs(uint32_t samples[NSHARES][(FRODO_N+31)/32][16]);

#endif
