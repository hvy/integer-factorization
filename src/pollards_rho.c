#include "pollards_rho.h"

#include <assert.h>

void pollard_rho(mpz_t factor, const mpz_t n, gmp_randstate_t rnd_state) {

  /* Assert n != 1 */
  assert(0 != mpz_cmp_ui(n, 1));

  mpz_t n_minus_1, a, xi, xm, s, tmp_abs;

  mpz_init(factor);
  mpz_init(xi);
  mpz_init(s);
  mpz_init(tmp_abs);
  mpz_init_set(n_minus_1, n);
  mpz_sub_ui(n_minus_1, n_minus_1, 1);
  mpz_init_set_ui(a, 1);
  mpz_urandomm(xi, rnd_state, n_minus_1);
  
  if(0 == mpz_cmp_ui(xi, 0)) {
    mpz_add_ui(xi, xi, 3);
  }

  mpz_init_set(xm, xi);

  for(int i = 0; i < 1000000; ++i) {
    mpz_mul(xi, xi, xi); 
    mpz_add(xi, xi, a); 
    mpz_mod(xi, xi, n);
    mpz_sub(tmp_abs, xi, xm);
    mpz_abs(tmp_abs, tmp_abs);
    mpz_gcd(s, tmp_abs, n);

    if(0 != mpz_cmp_ui(s, 1) && 0 != mpz_cmp(s, n) /* found a factor */) {
      mpz_set(factor, s);
      break;
    }

    if(0 == (i & (i - 1)) /* i is a power of 2 */) {
      mpz_set(xm, xi);  
    }
  }
  
  mpz_clear(n_minus_1);
  mpz_clear(a);
  mpz_clear(xi);
  mpz_clear(xm);
  mpz_clear(s);
  mpz_clear(tmp_abs);
}
