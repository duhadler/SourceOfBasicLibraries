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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#ifdef T

#include "flint.h"
#include "templates.h"
#include "profiler.h"

#define nalgs 2
#define ncases 20
#define cpumin 2

int
main(int argc, char** argv)
{
    double s[nalgs];

    int c, n, lenf, ext, reps = 0;
    fmpz_t p, temp;
    TEMPLATE(T, poly_t) f;
    TEMPLATE(T, ctx_t) ctx;
    
    FLINT_TEST_INIT(state);
    
    fmpz_init(p);
    fmpz_set_str(p, argv[1], 10);

    fmpz_init(temp);
       
    fmpz_set_str(temp, argv[2], 10);
    ext = fmpz_get_si(temp);

    lenf = atol(argv[3]);

    TEMPLATE(T, ctx_init)(ctx, p, ext, "a");

    TEMPLATE(T, poly_init)(f, ctx);

    for (c = 0; c < nalgs; c++)
    {
        s[c] = 0.0;
    }
       
    for (n = 0; n < ncases; n++)
    {
        double t[nalgs];
        int l, loops = 1;

        /*
           Construct random elements of fq
        */
        {
            TEMPLATE(T, poly_randtest_monic)(f, state, lenf, ctx);
        }
        
    loop:
        t[0] = 0.0;
        init_clock(0);
        prof_start();
        for (l = 0; l < loops; l++)
        {
            TEMPLATE(T, poly_is_irreducible_ben_or)(f, ctx);
        }
        prof_stop();
        t[0] += get_clock(0);

        t[1] = 0.0;
        init_clock(0);
        prof_start();
        for (l = 0; l < loops; l++)
        {
            TEMPLATE(T, poly_is_irreducible_ddf)(f, ctx);
        }
        prof_stop();
        t[1] += get_clock(0);

        for (c = 0; c < nalgs; c++)
            if (t[c] * FLINT_CLOCK_SCALE_FACTOR <= cpumin)
            {
                loops *= 10;
                goto loop;
            }
        
        for (c = 0; c < nalgs; c++)
            s[c] += t[c];
        reps += loops;
    }
        
    for (c = 0; c < nalgs; c++)
    {
        flint_printf("%20f ", s[c] / (double) reps);
        fflush(stdout);
    }
    printf("\n");
        
    TEMPLATE(T, poly_clear)(f, ctx);
    TEMPLATE(T, ctx_clear)(ctx);
    fmpz_clear(p);
    fmpz_clear(temp);

    FLINT_TEST_CLEANUP(state);
    
    return 0;
}

#endif
