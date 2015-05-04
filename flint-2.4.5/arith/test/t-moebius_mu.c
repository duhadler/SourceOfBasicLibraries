/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"
#include "profiler.h"


void check(fmpz_t n, int expected)
{
    int mu;

    mu = arith_moebius_mu(n);
    if (mu != expected)
    {
        flint_printf("FAIL:");
        fmpz_print(n);
        flint_printf("\n");
    }
}

int main(void)
{
    fmpz_t x;
    ulong p;
    slong i, j, k, l;
    FLINT_TEST_INIT(state);

    flint_printf("moebius_mu....");
    fflush(stdout);

    fmpz_init(x);

    for (i = -1000; i < 1000; i++)
    {
        fmpz_set_si(x, i);
        check(x, n_moebius_mu(FLINT_ABS(i)));
    }

    for (i = 0; i < 1000; i++)
    {
        fmpz_set_ui(x, 1);
        /* Product of some primes */
        k = n_randtest(state) % 10;
        l = n_randtest(state) % 10;
        for (j = 0; j < k; j++)
        {
            l += (n_randtest(state) % 10) + 1;
            fmpz_mul_ui(x, x, n_nth_prime(l+1));
        }

        check(x, (k % 2 ? -1 : 1));
        fmpz_neg(x, x);

        check(x, (k % 2 ? -1 : 1));
        fmpz_abs(x, x);

        /* No longer square-free */
        p = n_nth_prime(n_randtest(state) % 100 + 1);
        fmpz_mul_ui(x, x, p*p);
        check(x, 0);
    }

    fmpz_clear(x);

    

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
