#include "masked_poly_frommsg.h"

// this function outputs a boolean-shared polynomial
void masked_poly_frommsg(uint32_t r[NSHARES][KYBER_N], /*const*/ uint8_t msg[NSHARES][KYBER_SYMBYTES]) 
{
	size_t i;
  uint32_t qhalf[NSHARES], tmp[NSHARES];

  // set output polynomial to zero
  for (i = 0; i < KYBER_N; i++)
  {
    b_mask(&r[0][i], KYBER_N, 0);
  }

  // store masked qhalf
  b_mask(qhalf, 1, (KYBER_Q + 1) / 2);

  // iterate over polynomial
  for (i = 0; i < KYBER_N; i++)
  {
    if ((i%32) == 1)
    {
      b_srl1(tmp, 1, (uint32_t*)&msg[0][i/8], KYBER_SYMBYTES/4);
    } else if ((i%32) > 1) {
      b_srl1(tmp, 1, tmp, 1);
    }

    if ((i%32) == 0)
    {
      b_cmov(&r[0][i], KYBER_N, qhalf, 1, (uint32_t*)&msg[0][i/8], KYBER_SYMBYTES/4);
    } else {
      b_cmov(&r[0][i], KYBER_N, qhalf, 1, tmp, 1);
    }
  }
}