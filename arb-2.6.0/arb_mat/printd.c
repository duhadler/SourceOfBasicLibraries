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

#include "arb_mat.h"

void
arb_mat_printd(const arb_mat_t mat, long digits)
{
    long i, j;

    for (i = 0; i < arb_mat_nrows(mat); i++)
    {
        printf("[");

        for (j = 0; j < arb_mat_ncols(mat); j++)
        {
            arb_printd(arb_mat_entry(mat, i, j), digits);

            if (j < arb_mat_ncols(mat) - 1)
                printf(", ");
        }

        printf("]\n");
    }
}
