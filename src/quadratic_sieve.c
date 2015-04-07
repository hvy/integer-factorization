#include "quadratic_sieve.h"

#include <stdlib.h>
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

int factor_base(long *c, mpz_t *v, long b, mpz_t n, const int primec, mpz_t *primev) {

  mpz_t b_mpz;
  mpz_init_set_si(b_mpz, b);
 
  /* If b is larger than the largest precomputed prime, then we won't have
     enough primes up until b. If that is the case, return an error code. */ 
  if(0 < mpz_cmp(b_mpz, primev[primec - 1])) {
    return -1; // TODO Change this to a typedef enum error code
  }

  *c = 0;
  
  for(int i = 0; i < primec && 0 < mpz_cmp(b_mpz, primev[i]); ++i) {
    
    /* Compute the Legendre Symbol. Note that primev[i] needs to be an odd
       positive prime. */
    if(0 != mpz_cmp_ui(primev[i], 2) /* primev[i] is not 2, hence an odd prime */) {
      
      int legenre_symbol = mpz_legendre(n, primev[i]);
      if (1 == legenre_symbol) {
        mpz_init_set(v[*c], primev[i]);

        printf("\tAdded to factor base: ");
        mpz_out_str(stdout, 10, v[*c]);
        printf("\n");

        ++(*c);
      }
    }  
  }
  
  return 0;
}

void quadratic_sieve(mpz_t factor, mpz_t n, const int primec, mpz_t *primev) {

  mpz_init(factor);

  if (0 != mpz_perfect_square_p(n) /* n is a perfect square */) {
    mpz_sqrt(factor, n);
    return;
  }

  long b = smoothness_bound(n);
  long m = sieving_interval(b);
 
  printf("\tSmoothness bound, B: %ld\n", b);
  printf("\tSieving interval, M: %ld\n", m );
  
  long factor_basec;
  mpz_t *factor_basev = (mpz_t *) malloc(b * sizeof(mpz_t));
  factor_base(&factor_basec, factor_basev, b, n, primec, primev); 

  printf("\tSize of factor base: %ld\n", factor_basec);

  free(factor_basev);
}
