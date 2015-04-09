#ifndef QUADRATIC_SIEVE_H
#define QUADRATIC_SIEVE_H

#include <stdio.h>
#include <gmp.h>

void tonelli_shanks(mpz_t x, const mpz_t p, const mpz_t n);

/* Lifts a solution f(r) = 0 (mod p^k) to p^(k+1), ie. finds an x such that
   f(s) = 0 (mod p^(k+1)) and r = x mod p^k 
   @param x such that f(s) = 0 (mod p^(k+1)) (output)
   @param r the current solution
   @param p a prime number
   @param k the current prime 
   @param n the prime number to factorize
   @param n_floor_sqrt sqrt(n) floored */
void hensel_lift(mpz_t x, const mpz_t r, const mpz_t p, const int k, const mpz_t n, const mpz_t n_floor_sqrt);

void quadratic_sieve(mpz_t factor, mpz_t n, const int primec, mpz_t *primev);

void get_q(mpz_t q, const mpz_t x, const mpz_t n, const mpz_t sqrt_n);

#endif
