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

#include "arf.h"


int main()
{
    long iter;
    flint_rand_t state;

    printf("set_fmpq....");
    fflush(stdout);

    flint_randinit(state);

    /* test exact roundtrip R -> Q -> R */
    for (iter = 0; iter < 100000; iter++)
    {
        long bits, res;
        arf_t x, z;
        fmpq_t y;

        bits = 2 + n_randint(state, 200);

        arf_init(x);
        arf_init(z);
        fmpq_init(y);

        arf_randtest(x, state, bits, 10);
        arf_randtest(z, state, bits, 10);

        arf_get_fmpq(y, x);
        res = arf_set_fmpq(z, y, bits, ARF_RND_DOWN);

        if (!arf_equal(x, z) || res != 0)
        {
            printf("FAIL\n\n");
            printf("bits: %ld\n", bits);
            printf("x = "); arf_print(x); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            printf("z = "); arf_print(z); printf("\n\n");
            abort();
        }

        arf_clear(x);
        arf_clear(z);
        fmpq_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
