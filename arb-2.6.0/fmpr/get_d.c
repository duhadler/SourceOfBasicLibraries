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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmpr.h"

double
fmpr_get_d(const fmpr_t x, fmpr_rnd_t rnd)
{
    double r;
    mpfr_rnd_t mrnd;
    mpfr_t t;
    mp_limb_t tmp[2];

    if (fmpr_is_zero(x))
        return 0.0;

    mrnd = rnd_to_mpfr(rnd);

    t->_mpfr_prec = 53;
    t->_mpfr_sign = 1;
    t->_mpfr_exp = 0;
    t->_mpfr_d = tmp;

    fmpr_get_mpfr(t, x, mrnd);
    r = mpfr_get_d(t, mrnd);

    return r;
}

