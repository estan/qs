#include <iostream>
#include <vector>
#include <stack>
#include <ctime>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmpxx.h>

#include "fermat.h"
#include "rho.h"
#include "qs.h"

// Input size threshold below which we resort to trial division.
const static uint32_t TRIAL_THRESHOLD = 1000000000;

// Maximum input size we can handle (bits).
const static uint32_t MAX_DIGITS = 100;

// GMP random number generator.
extern gmp_randclass rng;

enum Algorithm {
    Fermat,
    PollardRho,
    QuadraticSieve,
    Heuristic
};

int main(int argc, char *argv[]) {

    // Determine the algorithm to use.
    Algorithm algorithm = Heuristic;
    bool badArg = false;
    if (argc == 2) {
        std::string arg(argv[1]);
        if (arg.compare("-f") == 0) {
            algorithm = Fermat;
        } else if (arg.compare("-p") == 0) {
            algorithm = PollardRho;
        } else if (arg.compare("-q") == 0) {
            algorithm = QuadraticSieve;
        } else {
            badArg = true;
        }
    }
    if (badArg) {
        std::cerr << "Usage: factor [-f | -p | -q ]" << std::endl;
        std::cerr << "    -f    Use only Fermat algorithm" << std::endl;
        std::cerr << "    -p    Use only Pollard's rho algorithm" << std::endl;
        std::cerr << "    -q    Use only Quadratic Sieve algorithm" << std::endl;
        std::cerr << "By default, Quadratic Sieve is used, but with" << std::endl;
        std::cerr << "fallback to trial division for numbers < 10^9." << std::endl;
        return 1;
    }

    // Seed the GMP random number generator.
    rng.seed(time(0));

    // Seed the standard library random number generator.
    srand(time(0));

    // Find some primes for trial division.
    std::vector<uint32_t> primes;
    uint32_t max = ceil(sqrt(TRIAL_THRESHOLD)) + 1;
    std::vector<bool> sieve(max, false);
    for (uint32_t p = 2; p < max; ++p) {
        if (sieve[p])
            continue;
        primes.push_back(p);
        for (uint32_t i = p; i < max; i += p)
            sieve[i] = true;
    }

    // Factor each of the 100 integers.
    mpz_class N;
    while (std::cin >> N) {
        if (mpz_sizeinbase(N.get_mpz_t(), 2) > MAX_DIGITS) {
            std::cout << "fail" << std::endl << std::endl; // Too many digits.
            continue;
        }

        if (mpz_probab_prime_p(N.get_mpz_t(), 10)) {
            // N is prime.
            std::cout << N << std::endl << std::endl;
            continue;
        }

        std::stack<mpz_class> factors;
        factors.push(N);

        while (!factors.empty()) {
            mpz_class factor = factors.top();
            factors.pop();

            if (mpz_probab_prime_p(factor.get_mpz_t(), 10)) {
                // N is prime.
                std::cout << factor << std::endl;
                continue;
            }

            if (algorithm == Fermat) {
                mpz_class result = fermat(factor);
                factors.push(result);
                factors.push(factor/result);
            } else if (algorithm == PollardRho) {
                mpz_class result = rho(factor);
                factors.push(result);
                factors.push(factor/result);
            } else if (algorithm == QuadraticSieve) {
                mpz_class result = quadraticSieve(factor);
                factors.push(result);
                factors.push(factor/result);
            } else {
                // Use a combination of techniques.

                if (factor < TRIAL_THRESHOLD) {
                    // Run trial division if below threshold.
                    uint32_t smallFactor = factor.get_ui();
                    for (uint32_t i = 0; i < primes.size(); ++i) {
                        if (smallFactor % primes[i] == 0) {
                            factors.push(primes[i]);
                            factors.push(factor / primes[i]);
                            break;
                        }
                    }
                } else {
                    // Run some trial division before starting the sieve.
                    bool foundFactor = false;
                    for (uint32_t i = 0; i < primes.size(); ++i) {
                        if (mpz_divisible_ui_p(factor.get_mpz_t(), primes[i])) {
                            factors.push(primes[i]);
                            factors.push(factor / primes[i]);
                            foundFactor = true;
                            break;
                        }
                    }
                    if (foundFactor)
                        continue; // Trial division was successful.

                    // Handle perfect powers separately (QS doesn't like them).
                    if (mpz_perfect_power_p(factor.get_mpz_t())) {
                        mpz_class root, r;
                        uint32_t max = mpz_sizeinbase(factor.get_mpz_t(), 2) / 2;
                        for (uint32_t n = 2; n < max; ++n) {
                            mpz_rootrem(root.get_mpz_t(), r.get_mpz_t(), factor.get_mpz_t(), n);
                            if (r == 0) {
                                for (uint32_t i = 0; i < n; ++i)
                                    factors.push(root);
                            }
                        }
                    } else {
                        // Run the QS algorithm.
                        mpz_class result = quadraticSieve(factor);
                        factors.push(result);
                        factors.push(factor / result);
                    }
                }
            }
        }
        std::cout << std::endl;
    }
    return 0;
}

