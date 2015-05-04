/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include <math.h>
#include "arb_poly.h"

void
_arb_poly_binomial_transform_convolution(arb_ptr b, arb_srcptr a, long alen, long len, long prec)
{
    long i;
    arb_ptr c, d;

    alen = FLINT_MIN(alen, len);

    c = _arb_vec_init(alen);
    d = _arb_vec_init(len);

    _arb_poly_borel_transform(c, a, alen, prec);
    for (i = 1; i < alen; i += 2)
        arb_neg(c + i, c + i);

    arb_one(d);
    for (i = 1; i < len; i++)
        arb_div_ui(d + i, d + i - 1, i, prec);

    _arb_poly_mullow(b, d, len, c, alen, len, prec);

    _arb_poly_inv_borel_transform(b, b, len, prec);

    _arb_vec_clear(c, alen);
    _arb_vec_clear(d, len);
}

void
arb_poly_binomial_transform_convolution(arb_poly_t b, const arb_poly_t a, long len, long prec)
{
    if (len == 0 || a->length == 0)
    {
        arb_poly_zero(b);
        return;
    }

    if (b == a)
    {
        arb_poly_t c;
        arb_poly_init2(c, len);
        _arb_poly_binomial_transform_convolution(c->coeffs, a->coeffs, a->length, len, prec);
        arb_poly_swap(b, c);
        arb_poly_clear(c);
    }
    else
    {
        arb_poly_fit_length(b, len);
        _arb_poly_binomial_transform_convolution(b->coeffs, a->coeffs, a->length, len, prec);
    }

    _arb_poly_set_length(b, len);
    _arb_poly_normalise(b);
}

