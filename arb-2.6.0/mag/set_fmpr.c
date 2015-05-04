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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "mag.h"
#include "arf.h"

void
mag_set_fmpr(mag_t x, const fmpr_t y)
{
    if (fmpr_is_special(y))
    {
        if (fmpr_is_zero(y))
            mag_zero(x);
        else
            mag_inf(x);
    }
    else
    {
        mag_set_fmpz_2exp_fmpz(x, fmpr_manref(y), fmpr_expref(y));
    }
}

