#include "quadratic_sieve.h"

#include <math.h>

/* Computes the smoothness bound B, given n. The smoothness bound B will be 
   the size of the factor base. The algorithm is inspired by Torbjorn 
   Granlund */
long smoothness_bound(mpz_t n) {
  double nd = mpz_get_d(n);
  double b = 2 * exp(0.5 * sqrt(log(nd) * log(log(nd))));
  return (long) b;
}

/* Computes the sieving interval M, given b. The algorithm is inspired 
   by Eric Landquist. The quadratic sieve factoring algorithm, 2001. */
long sieving_interval(long b) {
  return 4 * b * b;  
}

void factor_base(long factorc, mpz_t *factorv, long b, long n) {
  // TODO Continue working here  
}


void quadratic_sieve(mpz_t factor, mpz_t n) {

  mpz_init(factor);

  if (0 != mpz_perfect_square_p(n) /* n is a perfect square */) {
    mpz_sqrt(factor, n);
    return;
  }

  long b = smoothness_bound(n);
  long m = sieving_interval(b);
  
   

  
}
