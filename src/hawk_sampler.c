#include "hawk_sampler.h"
#include "core_ops.h"
#include "gadgets.h"



#if HAWK_N == 256
#define CDF_TABLE_LEN 20
static const uint16_t sig_gauss_hi_Hawk[CDF_TABLE_LEN] = {
	0x4D70, 0x268B,
	0x0F80, 0x04FA,
	0x0144, 0x0041,
	0x000A, 0x0001,
  0,0,0,0,
  0,0,0,0,
  0,0,0,0
};
static const uint64_t sig_gauss_lo_Hawk[CDF_TABLE_LEN] = {
	0x71FBD58485D45050, 0x1408A4B181C718B1,
	0x54114F1DC2FA7AC9, 0x614569CC54722DC9,
	0x42F74ADDA0B5AE61, 0x151C5CDCBAFF49A3,
	0x252E2152AB5D758B, 0x23460C30AC398322,
	0x0FDE62196C1718FC, 0x01355A8330C44097,
	0x00127325DDF8CEBA, 0x0000DC8DE401FD12,
	0x000008100822C548, 0x0000003B0FFB28F0,
	0x0000000152A6E9AE, 0x0000000005EFCD99,
	0x000000000014DA4A, 0x0000000000003953,
	0x000000000000007B, 0x0000000000000000
};



#elif HAWK_N == 512
#define CDF_TABLE_LEN 26
static const uint16_t sig_gauss_hi_Hawk[CDF_TABLE_LEN] = {
	0x580B, 0x35F9,
	0x1D34, 0x0DD7,
	0x05B7, 0x020C,
	0x00A2, 0x002B,
	0x000A, 0x0001,
  0,0,0,0,
  0,0,0,0,
  0,0,0,0,
  0,0,0,0
};
static const uint64_t sig_gauss_lo_Hawk[CDF_TABLE_LEN] = {
	0x0C27920A04F8F267, 0x3C689D9213449DC9,
	0x1C4FF17C204AA058, 0x7B908C81FCE3524F,
	0x5E63263BE0098FFD, 0x4EBEFD8FF4F07378,
	0x56AEDFB0876A3BD8, 0x4628BC6B23887196,
	0x061E21D588CC61CC, 0x7F769211F07B326F,
	0x2BA568D92EEC18E7, 0x0668F461693DFF8F,
	0x00CF0F8687D3B009, 0x001670DB65964485,
	0x000216A0C344EB45, 0x00002AB6E11C2552,
	0x000002EDF0B98A84, 0x0000002C253C7E81,
	0x000000023AF3B2E7, 0x0000000018C14ABF,
	0x0000000000EBCC6A, 0x000000000007876E,
	0x00000000000034CF, 0x000000000000013D,
	0x0000000000000006, 0x0000000000000000
};




#elif HAWK_N == 1024
#define CDF_TABLE_LEN 26
static const uint16_t sig_gauss_hi_Hawk[CDF_TABLE_LEN] = {
	0x58B0, 0x36FE,
	0x1E3A, 0x0EA0,
	0x0632, 0x024A,
	0x00BC, 0x0034,
	0x000C, 0x0002,
  0,0,0,0,
  0,0,0,0,
  0,0,0,0,
  0,0,0,0
};
static const uint64_t sig_gauss_lo_Hawk[CDF_TABLE_LEN] = {
	0x3AAA2EB76504E560, 0x01AE2B17728DF2DE,
	0x70E1C03E49BB683E, 0x6A00B82C69624C93,
	0x55CDA662EF2D1C48, 0x2685DB30348656A4,
	0x31E874B355421BB7, 0x430192770E205503,
	0x57C0676C029895A7, 0x5353BD4091AA96DB,
	0x3D4D67696E51F820, 0x09915A53D8667BEE,
	0x014A1A8A93F20738, 0x0026670030160D5F,
	0x0003DAF47E8DFB21, 0x0000557CD1C5F797,
	0x000006634617B3FF, 0x0000006965E15B13,
	0x00000005DBEFB646, 0x0000000047E9AB38,
	0x0000000002F93038, 0x00000000001B2445,
	0x000000000000D5A7, 0x00000000000005AA,
	0x0000000000000021, 0x0000000000000000
};
#endif

static const uint8_t LOGTABLE[] = {
  1, 2, 2, 3, 3, 3, 3, 4, 
  4, 4, 4, 4, 4, 4, 4, 5, 
  5, 5, 5, 5, 5, 5, 5, 5, 
  5, 5, 5, 5, 5, 5, 5, 6
};

/**
 * 80 bit slices represent 32 input values
 *
 * Usage: put the input randomness bitsliced into the 3d-array samples, then call the function.
 * Output will be written to the same array.
 */
void hawk_sampler_bs(uint32_t samples[NSHARES][(HAWK_N + 31)/32][80]) {
    size_t i, j, k;
    uint32_t tmp[NSHARES] = {0};

    for (i = 0; i < (HAWK_N + 31) / 32; i++) 
    {
        uint32_t result[NSHARES][5] = {0};
        for (j = 0; j < CDF_TABLE_LEN; j++) {
            uint32_t cmpres[NSHARES] = {0};
            for (k = 0; k < 80; k++)
            {
                b_not(tmp, 1,
                      &samples[0][i][k],
                      ((HAWK_N+31)/32)*80);
                if (k < 64)
                {
                  if ((sig_gauss_lo_Hawk[j]>>k)&1) {
                      b_and(cmpres, 1, cmpres, 1, tmp, 1);
                  } else {
                      b_or(cmpres, 1, cmpres, 1, tmp, 1);
                  }
                } else {
                  if ((sig_gauss_hi_Hawk[j]>>(k-64))&1) {
                      b_and(cmpres, 1, cmpres, 1, tmp, 1);
                  } else {
                      b_or(cmpres, 1, cmpres, 1, tmp, 1);
                  }
                }
            }
            if (j == 0) {
                b_not(result[0], 5, cmpres, 1);
            } else {
                b_not(cmpres, 1, cmpres, 1);
                // ripple-carry adder LOGTABLE[j] bit + 1 bit with tmp as carry
                b_and(tmp, 1, result[0], 5, cmpres, 1);
                b_xor(result[0], 5, result[0], 5, cmpres, 1);
                for (k = 1; k < (size_t) LOGTABLE[j]; k++) {
                    uint32_t tmp2[2] = {0};
                    b_xor(tmp2, 1, &result[0][k], 5, tmp, 1);
                    b_and(tmp, 1, &result[0][k], 5, tmp, 1);
                    copy_sharing(NSHARES, &result[0][k], 5, tmp2, 1);
                }
            }
        }
        for (j = 0; j < 5; j++) 
        {
          copy_sharing(NSHARES, &samples[0][i][j], ((HAWK_N+31)/32)*80, &result[0][j], 5);
        }

        // zeroize remaining random input
        for (k = 0; k < NSHARES; k++)
        {
          for (j = 5; j < 80; j++)
          {
            samples[k][i][j] = 0;
          }
        }
    }
}
