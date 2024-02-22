#include "bike_rotr.h"
#include "core_ops.h"
#include "masked.h"
#include "masked_cmp.h"
#include "masked_fwsampling.h"
#include "masked_keccak.h"
#include "masked_poly_frommsg.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

static void b_unmask_poly(uint32_t *p, size_t len) {
    for (size_t i = 0; i < len; i++) {
        p[0 + i] ^= p[len + i];
    }
}


/*static void print_array(uint32_t* p, size_t len) {
    for(size_t i=0; i<len; i++) {
        printf("%u, ", p[i]);
    }
    printf("\n");
}*/

// verifies if all indices are between 0 and N-1 and unique
// returns 0 for success, and 1 otherwise
static uint32_t verify_index_poly(uint32_t p[W]) {
    for (int i = 0; i < W; i++) {
        if (p[i] > N - 1) {
            printf("ERROR: Invalid Index: %u\n", p[i]);
            return 1;
        }

        for (int j = i + 1; j < W; j++) {
            if (p[i] == p[j]) {
                printf("ERROR: Index Collision\n");
                return 1;
            }
        }
    }
    return 0;
}

// verifies if all coefficients are 0 or 1 and sum up to W
// returns 0 for success, and 1 otherwise
static uint32_t verify_coeff_poly(uint32_t *p) {
    uint32_t weight = 0;
    uint32_t coeff;

    for (int i = 0; i < N; i++) {
        coeff = p[i] & 1;
        if (coeff > 1) {
            printf("ERROR: Invalid Coefficient: %u\n", coeff);
            return 1;
        }
        weight += coeff;
    }

    if (weight != W) {
        printf("ERROR: Incorrect Weight: %u\n", weight);
        return 1;
    }
    return 0;
}

// verifies if the weight is W
// returns 0 for success, and 1 otherwise
static uint32_t verify_dense_coeff_poly(uint32_t *p) {
    uint32_t weight = 0;

    uint32_t mask = 0xFFFFFFFF >> (PAD32(N) - N);
    p[((N + 31) / 32) - 1] &= mask;

    for (int i = 0; i < (N + 31) / 32; i++) {
        for (int j = 0; j < 32; j++) {
            weight += (p[i] >> j) & 1;
        }
    }

    if (weight != W) {
        printf("ERROR: Incorrect Weight: %u (dense)\n", weight);
        return 1;
    }
    return 0;
}

static uint32_t verify_sort_order(uint32_t *p, int len) {
    for (int i = 1; i < len; i++) {

        uint32_t val0 = p[i - 1] & 0xFFFFFFFC;
        uint32_t val1 = p[i - 0] & 0xFFFFFFFC;
        if (val1 < val0) {
            printf("ERROR: Incorrect sorting\n");
            return 1;
        }
    }
    return 0;
}

static int test_reject() {
    uint32_t p_bs[2][W_PADDED];
    uint32_t p[2][W];
    uint32_t rand[2][B_PADDED];
    int ret;

    randombytes((uint8_t *)&rand[0][0], sizeof(uint32_t) * 2 * B_PADDED);
    ret = bs_sample_reject(p_bs, rand);
    if (ret)
        printf("ERROR: BS Rejection reached end of randomness\n");
    b_unmask_poly(&p_bs[0][0], W_PADDED);
    ret |= verify_index_poly(&p_bs[0][0]);
    if (ret)
        return ret;

    randombytes((uint8_t *)&rand[0][0], sizeof(uint32_t) * 2 * B_PADDED);
    ret = sample_reject(p, rand);
    if (ret)
        printf("ERROR: Rejection reached end of randomness\n");
    b_unmask_poly(&p[0][0], W);
    ret |= verify_index_poly(&p[0][0]);

    return ret;
}

static int test_fisher_yates() {
    uint32_t p[2][W_PADDED];
    int ret;

    randombytes((uint8_t *)&p[0][0], sizeof(uint32_t) * 2 * W_PADDED);
    bs_sample_fisher_yates(p);
    b_unmask_poly(&p[0][0], W_PADDED);
    ret = verify_index_poly(&p[0][0]);
    if (ret)
        return ret;

    randombytes((uint8_t *)&p[0][0], sizeof(uint32_t) * 2 * W_PADDED);
    sample_fisher_yates(p);
    b_unmask_poly(&p[0][0], W_PADDED);
    ret = verify_index_poly(&p[0][0]);

    return ret;
}

static int test_sort() {
    uint32_t p_bs[2][PAD64(N)];
    uint32_t p[2][N];
    int ret;

    randombytes((uint8_t *)&p_bs[0][0], sizeof(uint32_t) * 2 * N);
    bs_sample_sort(p_bs);
    b_unmask_poly(&p_bs[0][0], PAD64(N));
    ret = verify_coeff_poly(&p_bs[0][0]);
    ret |= verify_sort_order(&p_bs[0][0], PAD64(N));
    if (ret)
        return ret;

    randombytes((uint8_t *)&p[0][0], sizeof(uint32_t) * 2 * N);
    sample_sort(p);
    b_unmask_poly(&p[0][0], N);
    ret = verify_coeff_poly(&p[0][0]);
    ret |= verify_sort_order(&p[0][0], N);

    return ret;
}

static int test_repand() {
    uint32_t p[2][(N + 31) / 32];
    int ret;

    memset(&p[0][0], 0, sizeof(uint32_t) * 2 * ((N + 31) / 32));
    bs_sample_repand(p);
    b_unmask_poly(&p[0][0], (N + 31) / 32);
    ret = verify_dense_coeff_poly(&p[0][0]);

    return ret;
}

static int test_bike_rotr() {
    uint32_t in[NSHARES][3 * R_QWORDS];
    uint32_t out[NSHARES][3 * R_QWORDS];
    uint32_t qw_num[NSHARES];
    uint32_t rot;

    randombytes((uint8_t *)&in[0][0], sizeof(uint32_t) * R_QWORDS);
    randombytes((uint8_t *)&in[1][0], sizeof(uint32_t) * R_QWORDS);
    memcpy(&in[0][1 * R_QWORDS], &in[0][0], sizeof(uint32_t) * R_QWORDS);
    memcpy(&in[0][2 * R_QWORDS], &in[0][0], sizeof(uint32_t) * R_QWORDS);
    memcpy(&in[1][1 * R_QWORDS], &in[1][0], sizeof(uint32_t) * R_QWORDS);
    memcpy(&in[1][2 * R_QWORDS], &in[1][0], sizeof(uint32_t) * R_QWORDS);

    memcpy(&out[0][0], &in[0][0], sizeof(uint32_t) * 6 * R_QWORDS);

    rot = randomint();
    rot &= (UPTOPOW2(R_QWORDS) - 1); // log2(rot) <= log2(R_QWORDS)
    b_mask(qw_num, 1, rot);

    rotr_big(out, qw_num);

    b_unmask_poly(&in[0][0], 3 * R_QWORDS);
    b_unmask_poly(&out[0][0], 3 * R_QWORDS);

    for(int i=0; i < R_QWORDS; i++) {
        if(out[0][i] != in[0][i + rot]) {
            printf("ERROR: BIKE rotation failed\n");
            return 1;
        }
    }
    return 0;
}

static int test_cmp() {
    uint32_t randomness[NSHARES][N], p[NSHARES][(N + 31) / 32];
    int ret;

    do {
      randombytes((uint8_t*)randomness[0], NSHARES * N * sizeof(uint32_t));
      // ensure that random input is < 2^L
      for (size_t i = 0; i < N; i++)
      {
        randomness[0][i] &= (1<<L)-1;
        for (size_t j = 1; j < NSHARES; j++)
        {
          randomness[0][i] ^= (~((1<<L)-1)) & randomness[j][i];
        }
      }
    } while (sample_cmp(p, randomness, round((((uint64_t)W)<<L)*1.0f/(N*1.0f))) != W);

    b_unmask_poly(p[0], (N+31)/32);
    ret = verify_dense_coeff_poly(p[0]);
    return ret;
}

static int test_bs_cmp() {
    uint32_t randomness[PAD32(N)/32][NSHARES][L], p[NSHARES][(N + 31) / 32];
    int ret;

    do {
      randombytes((uint8_t*)randomness[0][0], NSHARES * PAD32(N)/32 * L * sizeof(uint32_t));
    } while (bs_sample_cmp(p, randomness) != (W %256));

    b_unmask_poly(p[0], (N+31)/32);
    ret = verify_dense_coeff_poly(p[0]);
    return ret;
}

static int test_i2c() {
    uint32_t p[2][W_PADDED], pp[2][N_PADDED/32];
    int ret;

    randombytes((uint8_t *)&p[0][0], sizeof(uint32_t) * 2 * W_PADDED);
    bs_sample_fisher_yates(p);
    coeff_to_regular(pp, p);
    b_unmask_poly(&pp[0][0], N_PADDED/32);
    ret = verify_dense_coeff_poly(&pp[0][0]);
    if (ret)
        return ret;

    randombytes((uint8_t *)&p[0][0], sizeof(uint32_t) * 2 * W_PADDED);
    sample_fisher_yates(p);
    coeff_to_regular(pp, p);
    b_unmask_poly(&pp[0][0], N_PADDED/32);
    ret = verify_dense_coeff_poly(&pp[0][0]);

    return ret;
}

static int test_bs_i2c() {
    uint32_t p[2][W_PADDED], pp[2][N_PADDED/32];
    int ret;

    randombytes((uint8_t *)&p[0][0], sizeof(uint32_t) * 2 * W_PADDED);
    bs_sample_fisher_yates(p);
    bs_coeff_to_regular(pp, p);
    b_unmask_poly(&pp[0][0], N_PADDED/32);
    ret = verify_dense_coeff_poly(&pp[0][0]);
    if (ret)
        return ret;

    randombytes((uint8_t *)&p[0][0], sizeof(uint32_t) * 2 * W_PADDED);
    sample_fisher_yates(p);
    bs_coeff_to_regular(pp, p);
    b_unmask_poly(&pp[0][0], N_PADDED/32);
    ret = verify_dense_coeff_poly(&pp[0][0]);

    return ret;
}

static int test_keccak() {
  const uint8_t in[2][4] = {{0xab, 0xad, 0x1d, 0xea}, {0,0,0,0}};
  const uint8_t inlong[2][140] = {{
    0x4a, 0xe0, 0x7f, 0x34, 0x16, 0x8c, 0x09, 0x64, 
    0x61, 0xa4, 0x85, 0x04, 0x8f, 0xa3, 0x32, 0xfd, 
    0x3e, 0x94, 0xdb, 0xe3, 0x17, 0x17, 0xc7, 0xd2, 
    0x84, 0x44, 0xa7, 0xfd, 0xce, 0xc8, 0x5b, 0xc5, 
    0x3d, 0x61, 0x70, 0xc0, 0xcd, 0xbc, 0xbd, 0xd9, 
    0xac, 0x3b, 0x13, 0x2c, 0x1a, 0x50, 0x7a, 0x34, 
    0xca, 0x50, 0xe6, 0x6c, 0xa2, 0xb5, 0x69, 0x16, 
    0x3e, 0xfa, 0x94, 0x01, 0xf7, 0x2f, 0xe0, 0xb0, 
    0x6a, 0x4d, 0xee, 0xd3, 0x90, 0x50, 0xca, 0x2d, 
    0x02, 0xf9, 0xcf, 0xe0, 0xe7, 0x31, 0x7d, 0x9d, 
    0x2f, 0x4b, 0x79, 0xf0, 0x5a, 0x51, 0x05, 0x0a, 
    0x2b, 0x4f, 0x1e, 0x4c, 0x39, 0xab, 0x3b, 0xfe, 
    0x6e, 0x8f, 0xff, 0x59, 0xa0, 0x98, 0xf3, 0x5d, 
    0x89, 0x11, 0xb2, 0xa7, 0x9e, 0x77, 0x06, 0x83, 
    0xe2, 0x5e, 0xd9, 0x4a, 0x5a, 0xa5, 0x28, 0x9e, 
    0xae, 0xd8, 0x02, 0x92, 0x74, 0x17, 0xef, 0x8a, 
    0x02, 0x14, 0xaa, 0xab, 0x1e, 0xbf, 0x29, 0x97, 
    0xa2, 0x2b, 0x0f, 0}, {0}};
  uint8_t out[2][64];
  const uint8_t out_expected_256[2][32] = {{
    0xe7, 0x0e, 0x9f, 0x95, 0x8c, 0x5a, 0xab, 0x53, 
    0x8d, 0x7b, 0x31, 0x83, 0x9b, 0x53, 0xb8, 0x73, 
    0x35, 0xdf, 0x0e, 0x8c, 0xbd, 0xe7, 0xfe, 0x55, 
    0x9c, 0x68, 0xa9, 0x51, 0xe3, 0x8d, 0x7f, 0x79
    }, {
    0xc6, 0x3e, 0xee, 0x29, 0xfe, 0x34, 0xf7, 0x0b, 
    0x2c, 0xf6, 0x72, 0x01, 0xfd, 0x0d, 0x7d, 0x14, 
    0x15, 0x70, 0x85, 0x01, 0x80, 0x0f, 0x0e, 0x1b, 
    0x99, 0x38, 0x4d, 0x32, 0xb2, 0x8f, 0x52, 0x71
    }};
  
  const uint8_t out_expected_512[2][64] = {{
    0x01, 0x6b, 0x9a, 0x09, 0x1e, 0xa6, 0x31, 0xfc, 
    0x7b, 0x3f, 0x14, 0x18, 0x4a, 0x55, 0x4c, 0xb0, 
    0x67, 0x10, 0xc1, 0xc1, 0x7f, 0x5f, 0x3d, 0xb7, 
    0x56, 0x2b, 0x0a, 0xfc, 0x72, 0xdf, 0x23, 0xc9, 
    0xf8, 0xba, 0xca, 0x82, 0x0a, 0xeb, 0xae, 0xe3, 
    0x46, 0x2b, 0xf6, 0x6b, 0x19, 0xb2, 0x11, 0x94, 
    0x06, 0x07, 0x59, 0x82, 0xbd, 0x80, 0x7f, 0x77, 
    0x5b, 0x6c, 0x26, 0x28, 0xcd, 0x12, 0xcf, 0xbb
  },{
    0xfe, 0x2c, 0x1c, 0x5a, 0xe2, 0x43, 0x81, 0x13, 
    0x87, 0x27, 0x0a, 0xb6, 0xe5, 0x99, 0x46, 0xb2, 
    0x45, 0xe5, 0x69, 0xcb, 0xa9, 0xcf, 0x30, 0xe6, 
    0xc1, 0x3d, 0xf1, 0xda, 0x9b, 0x8c, 0xf3, 0x31, 
    0x89, 0x6a, 0xee, 0x56, 0xf4, 0x98, 0xbb, 0x7d, 
    0x11, 0x16, 0x8f, 0x90, 0x42, 0x42, 0xe9, 0xbd, 
    0x24, 0xea, 0xdd, 0x5b, 0x9d, 0xcb, 0x46, 0x5d, 
    0xe6, 0x50, 0xc9, 0x29, 0xbc, 0x5c, 0xb9, 0xeb
  }};

  const uint8_t out_long_expected_256[32] = {
    0xca, 0xe8, 0x17, 0x94, 0x88, 0xfe, 0xfd, 0x70, 
    0x64, 0x97, 0x20, 0xf5, 0x41, 0x3f, 0x0f, 0xf1, 
    0xc1, 0xac, 0x62, 0x15, 0xea, 0xa7, 0xe1, 0x0c, 
    0x19, 0x5e, 0x4e, 0x5c, 0x95, 0x53, 0x32, 0xe4
  };

  const uint8_t out_long_expected_512[64] = {
    0xee, 0xfa, 0x55, 0x82, 0x4d, 0xad, 0xfc, 0xb9, 
    0x4d, 0x2e, 0x2f, 0x7d, 0x1c, 0xe0, 0x91, 0xa4, 
    0x08, 0x8b, 0x0d, 0xb1, 0x3e, 0x65, 0xf7, 0x20, 
    0x21, 0x66, 0x9b, 0x03, 0xc8, 0xe2, 0xe7, 0x52, 
    0x0e, 0x15, 0x46, 0xde, 0x8f, 0x0a, 0x4c, 0x4b, 
    0x1e, 0x0b, 0x9a, 0x59, 0x7e, 0x45, 0x6c, 0xe1, 
    0xe4, 0x27, 0x08, 0x51, 0x55, 0x19, 0xaa, 0xc7, 
    0xfa, 0xde, 0x2c, 0x13, 0x13, 0xff, 0xb3, 0x61,
  };

  const uint8_t out_long_expected_shake256[32] = {
    0x54, 0xdb, 0x65, 0xbd, 0x1b, 0x7b, 0x0e, 0x43, 
    0x8a, 0xae, 0x52, 0x1a, 0x0f, 0x2c, 0xa4, 0x17, 
    0x00, 0x4a, 0xd5, 0x41, 0x42, 0xc2, 0x82, 0x08, 
    0x3c, 0xd7, 0x73, 0x4e, 0x1c, 0x81, 0x24, 0x3d
  };

  //////////////////////////////////////////////// SHA3-256
  masked_sha3_256(out[0], 64, 1, in[0], 4, 4, 1);

  for (size_t i = 0; i < 32; i++)
  {
    if ((out[0][i]^out[1][i]) != out_expected_256[0][i])
    {
      printf("bad sha3 256 (4 byte in)\n");
      return 1;
    }
  }

  masked_sha3_256(out[0], 64, 1, in[0], 3, 4, 1);

  for (size_t i = 0; i < 32; i++)
  {
    if ((out[0][i]^out[1][i]) != out_expected_256[1][i])
    {
      printf("bad sha3 256 (3 byte in)\n");
      return 1;
    }
  }

  masked_sha3_256(out[0], 64, 1, inlong[0], 139, 140, 1);

  for (size_t i = 0; i < 32; i++)
  {
    if ((out[0][i]^out[1][i]) != out_long_expected_256[i])
    {
      printf("bad sha3 256 (139 byte in)\n");
      return 1;
    }
  }

  //////////////////////////////////////////////// SHA3-512
  masked_sha3_512(out[0], 64, 1, in[0], 4, 4, 1);

  for (size_t i = 0; i < 64; i++)
  {
    if ((out[0][i]^out[1][i]) != out_expected_512[0][i])
    {
      printf("bad sha3 512 (4 byte in)\n");
      return 1;
    }
  }

  masked_sha3_512(out[0], 64, 1, in[0], 3, 4, 1);

  for (size_t i = 0; i < 64; i++)
  {
    if ((out[0][i]^out[1][i]) != out_expected_512[1][i])
    {
      printf("bad sha3 512 (3 byte in)\n");
      return 1;
    }
  }

  masked_sha3_512(out[0], 64, 1, inlong[0], 139, 140, 1);

  for (size_t i = 0; i < 64; i++)
  {
    if ((out[0][i]^out[1][i]) != out_long_expected_512[i])
    {
      printf("bad sha3 512 (139 byte in)\n");
      return 1;
    }
  }


  //////////////////////////////////////////////// SHAKE

  for (size_t olen = 1; olen < 32; olen++)
  {
    masked_shake256(out[0], olen, 64, 1, inlong[0], 139, 140, 1);

    for (size_t i = 0; i < olen; i++)
    {
      if ((out[0][i]^out[1][i]) != out_long_expected_shake256[i])
      {
        printf("bad shake256 (139 byte in, %lu byte out)\n", olen);
        return 1;
      }
    }
    //printf("+");
  }
  //printf("\n");
  return 0;
}


static int test_poly_frommsg(void)
{
  uint32_t r[NSHARES][KYBER_N];
  uint8_t msg[NSHARES][KYBER_SYMBYTES];

  randombytes(msg[0], NSHARES * KYBER_SYMBYTES);
  masked_poly_frommsg(r, msg);

  for (size_t i = 0; i < KYBER_N; i++)
  {
    r[0][i] ^= r[1][i];
    if (r[0][i] != ((KYBER_Q + 1) / 2) * (((msg[0][i/8]^msg[1][i/8])>>(i%8))&1))
    {
      printf("%lu: bit is %d, coefficient is %d\n", i, (((msg[0][i/8]^msg[1][i/8])>>(i%8))&1), r[0][i]);
      return 1;
    }
  }
  return 0;
}

int main() {
    int ret = 0;

    ret |= test_keccak();
    for (int i = 0; i < 10; i++) {
        ret |= test_bike_rotr();
        ret |= test_reject();
        ret |= test_fisher_yates();
        ret |= test_sort();
        ret |= test_repand(); 
        ret |= test_cmp();
        ret |= test_bs_cmp();
        ret |= test_i2c();
        ret |= test_bs_i2c();
        ret |= test_poly_frommsg();

        if (ret)
            break;
        printf(".");
        fflush(stdout);
    }
    printf("\n%i\n", ret);
    return ret;
}
