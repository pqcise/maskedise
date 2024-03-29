/* Copyright 2022 UCLouvain, Belgium and PQM4 contributors
 *
 * This file is part of pqm4_masked.
 *
 * pqm4_masked is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, version 3.
 *
 * pqm4_masked is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * pqm4_masked. If not, see <https://www.gnu.org/licenses/>.
 *
 * Edited by (anonymized)
 */
#ifndef GADGETS_H
#define GADGETS_H

#include "masked.h"
#include <stddef.h>
#include <stdint.h>

// Atomic gadgets
void masked_xor(size_t nshares, uint32_t *out, size_t out_stride,
                const uint32_t *ina, size_t ina_stride, const uint32_t *inb,
                size_t inb_stride);
void masked_and(size_t nshares, uint32_t *z, size_t z_stride, const uint32_t *a,
                size_t a_stride, const uint32_t *b, size_t b_stride);
void masked_or(size_t nshares, uint32_t *z, size_t z_stride, const uint32_t *a,
               size_t a_stride, const uint32_t *b, size_t b_stride);
void masked_or_1bit(size_t nshares, uint32_t *z, size_t z_stride,
                    const uint32_t *a, size_t a_stride, const uint32_t *b,
                    size_t b_stride);
void masked_and_1bit(size_t nshares, uint32_t *z, size_t z_stride,
                     const uint32_t *a, size_t a_stride, const uint32_t *b,
                     size_t b_stride);
void copy_sharing(size_t nshares, uint32_t *out, size_t out_stride,
                  const uint32_t *in, size_t in_stride);
void RefreshIOS_rec(size_t nshares, size_t d, uint32_t *x, size_t x_msk_stride);
void seccompress(size_t nshares, uint32_t q, uint32_t c, uint32_t *out,
                 size_t out_msk_stride, size_t out_data_stride,
                 const int16_t *in, size_t in_msk_stride,
                 size_t in_data_stride);

uint32_t unmask_boolean(size_t nshares, const uint32_t *in, size_t in_stride);

// old = cond ? new : old
// expects cond to be one bit
void masked_sel(const uint32_t *cond, size_t cond_stride, uint32_t *old,
                size_t old_stride, const uint32_t *new, size_t new_stride);
// old = cond ? new : old
// bits <= 32
// expects datastride of 1
void masked_bs_sel(uint32_t *cond, size_t cond_stride, uint32_t *old,
                   size_t old_stride, uint32_t *new, size_t new_stride,
                   size_t bits);

// returns 1 if equal
void masked_eq(uint32_t *res, size_t res_stride, const uint32_t *a,
               size_t a_stride, const uint32_t *b, size_t b_stride,
               size_t bits);
// bits <= 32
// expects datastride of 1
// returns 1 per element if equal
void masked_bs_eq(uint32_t *res, size_t res_stride, uint32_t *a,
                  size_t a_stride, uint32_t *b, size_t b_stride, size_t bits);
// returns 1 if b > a
void cmpg(uint32_t *res, uint32_t stride_res, const uint32_t *a,
          uint32_t stride_a, const uint32_t *b, uint32_t stride_b, size_t bits);
// returns 1 if b > a
void add(uint32_t *res, uint32_t stride_res, const uint32_t *a,
          uint32_t stride_a, const uint32_t *b, uint32_t stride_b, size_t bits);
// returns 1 if b > a
void bscmpg(uint32_t *res, const uint32_t *a, const uint32_t stride_a,
            const uint32_t *b, const uint32_t stride_b);

void masked_shiftr(uint32_t *z, size_t z_stride, const uint32_t *a,
                   size_t a_stride, size_t shiftamt);

// Adders
void secfulladd(size_t nshares, uint32_t *co, size_t co_msk_stride, uint32_t *w,
                size_t w_msk_stide, const uint32_t *ci, size_t ci_msk_stride,
                const uint32_t *x, size_t x_msk_stride, const uint32_t *y,
                size_t y_msk_stride);
void sechalfadd(size_t nshares, uint32_t *c, size_t c_msk_stride, uint32_t *w,
                size_t w_msk_stide, const uint32_t *x, size_t x_msk_stride,
                const uint32_t *y, size_t y_msk_stride);
void sechalfsub(size_t nshares, uint32_t *c, size_t c_msk_stride, uint32_t *w,
                size_t w_msk_stide, const uint32_t *x, size_t x_msk_stride,
                const uint32_t *y, size_t y_msk_stride);
void secadd(size_t nshares, size_t kbits, size_t kbits_out, uint32_t *out,
            size_t out_msk_stride, size_t out_data_stride, const uint32_t *in1,
            size_t in1_msk_stride, size_t in1_data_stride, const uint32_t *in2,
            size_t in2_msk_stride, size_t in2_data_stride);
void secadd_modp(size_t nshares, size_t kbits, uint32_t p, uint32_t *out,
                 size_t out_msk_stride, size_t out_data_stride,
                 const uint32_t *in1, size_t in1_msk_stride,
                 size_t in1_data_stride, const uint32_t *in2,
                 size_t in2_msk_stride, size_t in2_data_stride);
void secadd_constant_bmsk(size_t nshares, size_t kbits, size_t kbits_out,
                          uint32_t *out, size_t out_msk_stride,
                          size_t out_data_stride, const uint32_t *in1,
                          size_t in1_msk_stride, size_t in1_data_stride,
                          uint32_t constant, const uint32_t *bmsk,
                          size_t bmsk_msk_stride);
void secadd_constant(size_t nshares, size_t kbits, size_t kbits_out,
                     uint32_t *out, size_t out_msk_stride,
                     size_t out_data_stride, const uint32_t *in1,
                     size_t in1_msk_stride, size_t in1_data_stride,
                     uint32_t constant);

// Conversions
void seca2b(size_t nshares, size_t kbits, uint32_t *in, size_t in_msk_stride,
            size_t in_data_stride);
void seca2b_modp(size_t nshares, size_t kbits, uint32_t p, uint32_t *in,
                 size_t in_msk_stride, size_t in_data_stride);
void secb2a_1bit(size_t nshares, int16_t *a, size_t a_msk_stride, uint32_t *x,
                 size_t x_msk_stride);
void secb2a_modp(size_t nshares, uint32_t p, uint32_t *in, size_t in_msk_stride,
                 size_t in_data_stride);
void secb2a(size_t nshares, size_t kbits, uint32_t *in, size_t in_msk_stride,
            size_t in_data_stride);

// FW Sampling specific
void secb2a_32to48(size_t nshares, uint32_t *in, size_t in_msk_stride,
                   size_t in_data_stride);

#endif // GADGETS_H
