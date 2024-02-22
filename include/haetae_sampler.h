#ifndef HAETAE_SAMPLER_H
#define HAETAE_SAMPLER_H

#include <stdint.h>
#include "masked.h"

#define HAETAE_N 256
void haetae_sampler(uint16_t samples[NSHARES * HAETAE_N]);
void haetae_sampler_bs(uint32_t samples[NSHARES][(HAETAE_N+31)/32][16]);

#endif
