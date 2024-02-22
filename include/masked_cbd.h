#ifndef MASKED_CBD__H
#define MASKED_CBD__H

void masked_cbd2(uint32_t poly[NSHARES][N], uint32_t randomness[NSHARES][(N*2+31)/32]);
void masked_cbd3(uint32_t poly[NSHARES][N], uint32_t randomness[NSHARES][(N*3+31)/32]);

#endif
