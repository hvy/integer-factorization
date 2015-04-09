#ifndef QUADRATIC_SIEVE_H
#define QUADRATIC_SIEVE_H

#include <stdio.h>
#include <gmp.h>

void tonelli_shanks(mpz_t x, const mpz_t p, const mpz_t n);
void quadratic_sieve(mpz_t factor, mpz_t n, const int primec, mpz_t *primev);

#endif
