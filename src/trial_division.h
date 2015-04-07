#ifndef TRIAL_DIVISION_H
#define TRIAL_DIVISION_H

#include <stdio.h>
#include <gmp.h>

void get_fst_primes(const int primec, mpz_t *primev);
void trial_division(mpz_t factor, const mpz_t n, const int primec, mpz_t *primev);

#endif
