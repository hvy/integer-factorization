#include "trial_division.h"

void get_fst_primes(const int primec, mpz_t *primev) {

  int i; /* number of primes found */
  mpz_t n;
  
  /* Add first prime number, 2 */
  mpz_init_set_ui(primev[0], 2);
  i = 1;

  /* Initial prime candidate is set to 3 */
  mpz_init_set_ui(n, 3);
  
  while (i < primec) {
    if(mpz_probab_prime_p(n, 25)) {
      mpz_init_set(primev[i], n);
      ++i;
    }
    mpz_add_ui(n, n, 2);
  }

  mpz_clear(n);
}

void trial_division(mpz_t factor, const mpz_t n, const int primec, mpz_t *primev) {

  mpz_t tmp_mod;
  mpz_init(factor);
  mpz_init(tmp_mod);
  
  for(int i = 0; i < primec; ++i) {
    if(0 > mpz_cmp(n, primev[i]) /* n > primev[i] */) {
      break;
    }  

    mpz_mod(tmp_mod, n, primev[i]);
    
    if(0 == mpz_cmp_ui(tmp_mod, 0) /* primev[i] divides n */) {
      mpz_set(factor, primev[i]);
      break;
    }
  }

  mpz_clear(tmp_mod);
}
