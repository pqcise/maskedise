#ifndef HAWK_SAMPLER_H
#define HAWK_SAMPLER_H

#include <stdint.h>
#include "masked.h"

void hawk_sampler_bs(uint32_t samples[NSHARES][(HAWK_N+31)/32][80]);

#endif