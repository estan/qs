#include <gmpxx.h>
#include <gmp.h>

#include "fermat.h"

/*
 * The basic Fermat method of factorization.
 *
 * Takes an odd integer as input N and returns a factor of N.
 */
mpz_class fermat(const mpz_class& N) {
    mpz_class a = sqrt(N);
    mpz_class b2 = a * a - N;

    while (!mpz_perfect_square_p(b2.get_mpz_t())) {
        ++a;
        b2 = a * a - N;
    }
    return a - sqrt(b2);
}
