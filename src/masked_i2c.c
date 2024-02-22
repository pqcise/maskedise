#include <stdint.h>
#include <stddef.h>
#include "core_ops.h"
#include "gadgets.h"
#include "masked_fwsampling.h"

// expects non bitsliced inputs
void coeff_to_regular(uint32_t poly[NSHARES][N_PADDED/32], uint32_t index[NSHARES][W_PADDED]) {
  size_t i,j;
  uint32_t cntmsk[NSHARES], cmpres[NSHARES], one[NSHARES], shifted_one[NSHARES], current_or_one[NSHARES];
  b_mask(one, 1, 1);
  for (i = 0; i < N_PADDED/32; i++) // init to zero
  {
    b_mask(&poly[0][i], N_PADDED/32, 0);
  }

  for (i = 0; i < N; i++)
  {
    if ((i%32) == 0) {
      copy_sharing(NSHARES, shifted_one, 1, one, 1);
    }
    b_or(current_or_one, 1, &poly[0][i/32], N_PADDED/32, shifted_one, 1);
    b_mask(cntmsk, 1, i);
    for (j = 0; j < W; j++)
    {
#if defined(USE_MASKED_ISA) && !defined(USE_MASKED_EXT)
      asm volatile (
      "lw a0, (%[in0])\n"
      "lw a1, (%[in1])\n"
      "lw a2, (%[in2])\n"
      "lw a3, (%[in3])\n"

      "mask.b.sub (a5,a4),(a3,a2),(a1,a0)\n"
      "mask.b.srli (a5,a4),(a5,a4),31\n"
      "mask.b.sub (a7,a6), (a1,a0),(a3,a2)\n"
      "mask.b.srli (a7,a6), (a7,a6), 31\n"
      "mask.b.ior (a5,a4),(a5,a4),(a7,a6)\n"
      "mask.b.not (a5,a4),(a5,a4)\n"
      "li t0, 1\n"
      "mask.b.mask (a7,a6),t0\n"
      "mask.b.and (a5,a4),(a5,a4),(a7,a6)\n"

      "mask.b.and (a5,a4), (a5,a4), (a7,a6)\n"
      "mask.b.sub (a5,a4), (a5,a4), (a7,a6)\n"
      "mask.b.not (a5,a4), (a5,a4)\n"

      "lw a0, (%[poly0])\n"
      "lw a1, (%[poly1])\n"
      "lw a2, (%[curr0])\n"
      "lw a3, (%[curr1])\n"
      
      "mask.b.xor (a7,a6),(a3,a2),(a1,a0)\n"
      "mask.b.and (t6,t5),(a7,a6),(a5,a4)\n"
      "mask.b.xor (a1,a0),(a1,a0),(t6,t5)\n"

      "sw a0, (%[poly0])\n"
      "sw a1, (%[poly1])\n"
      :
      : [in0] "r" (cntmsk), 
      [in1] "r" (cntmsk+1), 
      [in2] "r" (&index[0][j]),
      [in3] "r" (&index[0][j]+W_PADDED),
      [poly0] "r" (&poly[0][i/32]),
      [poly1] "r" (&poly[0][i/32]+N_PADDED/32),
      [curr0] "r" (current_or_one),
      [curr1] "r" (current_or_one+1)
      : "t6","t5","t0","a7","a6", "a5", "a4", "a3", "a2", "a1", "a0", "memory");
#elif defined USE_MASKED_EXT && defined USE_MASKED_ISA
      asm volatile (
      "lw a0, (%[in0])\n"
      "lw a1, (%[in1])\n"
      "lw a2, (%[in2])\n"
      "lw a3, (%[in3])\n"
      "mask.b.cmpeq (a5,a4),(a3,a2),(a1,a0)\n"
      "lw a0, (%[poly0])\n"
      "lw a1, (%[poly1])\n"
      "lw a2, (%[curr0])\n"
      "lw a3, (%[curr1])\n"
      "mask.b.cmov (a1,a0),(a3,a2),(a5,a4)\n"
      "sw a0, (%[poly0])\n"
      "sw a1, (%[poly1])\n"
      :
      : [in0] "r" (cntmsk), 
      [in1] "r" (cntmsk+1), 
      [in2] "r" (&index[0][j]),
      [in3] "r" (&index[0][j]+W_PADDED),
      [poly0] "r" (&poly[0][i/32]),
      [poly1] "r" (&poly[0][i/32]+N_PADDED/32),
      [curr0] "r" (current_or_one),
      [curr1] "r" (current_or_one+1)
      : "a5", "a4", "a3", "a2", "a1", "a0", "memory");
#else
      b_cmpeq(cmpres, 1, cntmsk, 1, &index[0][j], W_PADDED);
      b_cmov(&poly[0][i/32], N_PADDED/32, current_or_one, 1, cmpres, 1);
#endif
    }
    if ((i%32) != 31 && i < N-1) {
      b_sll1(shifted_one, 1, shifted_one, 1);
    }
  }
}

// expects non bitsliced inputs
void bs_coeff_to_regular(uint32_t poly[NSHARES][N_PADDED/32], uint32_t index[NSHARES][W_PADDED]) 
{
  uint32_t i_index[NSHARES][W*32], p_index[NSHARES][32] = {0};
  uint32_t res[NSHARES]; // result of comparison
  size_t i,j,s;

  // transpose masked indices 32 times repeatedly
  for (s = 0; s < NSHARES; s++) 
  {
    for (i = 0; i < W; i++) 
    {
      for (j = 0; j < 32; j++)
      {
        i_index[s][i*32 + j] = index[s][i];
      }
      transpose32(&i_index[s][i*32]);
    }
  }


  // init poly to zero
#if defined(USE_MASKED_ISA) || defined(USE_MASKED_ISE)
  for (i = 0; i < N_PADDED/32; i++)
  {
    b_mask(&poly[0][i], N_PADDED/32, 0);
  }
#else
  for (s = 0; s < NSHARES; s++)
  {
    for (i = 0; i < N_PADDED/32; i++)
    {
      poly[s][i] = 0;
    }
  }
#endif


  for (i = 0; i < N; i += 32) // iterate over whole polynomial
  {
    for (j = i; j < i+32; j++) // generate bitsliced indices
    { 
      p_index[0][j-i] = j;
    }
    transpose32(p_index[0]);
#if defined(USE_MASKED_ISA) || defined(USE_MASKED_ISE)
    for (j = 0; j < 32; j++) // masked instructions require already masked input
    {
      b_mask(&p_index[0][j], 32, p_index[0][j]);
    }
#endif
    
    for(j = 0; j < W; j++) // iterate over all nonzero indices
    {
      // index of coefficient equal to the current nonzero coefficient?
      masked_bs_eq(res, 1, p_index[0], 32, &i_index[0][j*32], W*32, LOG_N);
      // if yes, set coefficient in poly to 1
      masked_bs_sel(res, 1, &poly[0][i/32], N_PADDED/32, res, 1, 1);
    }
  }
}

