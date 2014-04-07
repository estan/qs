#include "rho.h"

// Interface to the GMP random number functions.
gmp_randclass rng(gmp_randinit_default);

/*
 * A (very) basic Pollard's rho method implementation.
 *
 * Takes an integer N as input and returns a factor of N.
 */
mpz_class rho(const mpz_class &N)
{
    if (N % 2 == 0)
        return 2;

    mpz_class c = rng.get_z_range(N);
    mpz_class x = rng.get_z_range(N);
    mpz_class y = x;
    mpz_class d = 1;
    mpz_class z;

    while (d == 1) {
        x = (x*x + c) % N;
        y = (y*y + c) % N;
        y = (y*y + c) % N;
        z = x - y;
        mpz_gcd(d.get_mpz_t(), z.get_mpz_t(), N.get_mpz_t());
    }
    
    return d;
}

