#include "trial_division.h"
#include "pollards_rho.h"
#include "quadratic_sieve.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>

void test_tonelli_shanks();
void test_hensel_lift();

int main() {
  
  printf("[INFO] Starting test... \n");

  test_tonelli_shanks();
  test_hensel_lift();

  printf("[INFO] Finished all tests without errors.\n");

  return 0;
}

void test_tonelli_shanks() {

  printf("[INFO] Starting Tonelli-Shanks...");
  
  mpz_t x, p, n;
  
  mpz_init(x);
  mpz_init(p);
  mpz_init(n);

  /* Test 1 */
  mpz_set_ui(p, 13);
  mpz_set_ui(n, 10);
  tonelli_shanks(x, p, n);
  assert(7 == mpz_get_si(x));
  
  /* Test 2 */
  mpz_set_ui(p, 7);
  mpz_set_ui(n, 90283);
  tonelli_shanks(x, p, n);
  assert(2 == mpz_get_si(x));

  /* Test 3 */
  mpz_set_ui(p, 23);
  mpz_set_ui(n, 15347);
  tonelli_shanks(x, p, n);
  assert(12 == mpz_get_si(x));

  printf(" done.\n");

  mpz_clear(x);
  mpz_clear(p);
  mpz_clear(n);
}


void test_hensel_lift() {

  printf("[INFO] Starting Hensel Lift...");
  
  int k = 1;
  mpz_t bi0, bi1, p, n, sqrt_n, i0, i1, q, p_pow;
  
  mpz_init(bi0);
  mpz_init(bi1);
  mpz_init(p);
  mpz_init(n);
  mpz_init(sqrt_n);
  mpz_init(i0);
  mpz_init(i1);
  mpz_init(q);
  mpz_init(p_pow);

  /* Test 1 */
  mpz_set_ui(p, 23);
  mpz_set_ui(n, 15347);
  mpz_set_ui(sqrt_n, 124);
  mpz_set_ui(i0, 2);
  mpz_set_ui(i1, 3);

  for(int i = 0; i < 5; ++i) {
    hensel_lift(bi0, i0, p, k, n, sqrt_n);
    hensel_lift(bi1, i1, p, k, n, sqrt_n);

    mpz_pow_ui(p_pow, p, (k + 1));
    
    get_q(q, bi0, n, sqrt_n);
    mpz_mod(q, q, p_pow);
    assert(0 == mpz_cmp_ui(q, 0));
    
    get_q(q, bi1, n, sqrt_n);
    mpz_mod(q, q, p_pow);
    assert(0 == mpz_cmp_ui(q, 0));

    mpz_set(i0, bi0);
    mpz_set(i1, bi1); 
    ++k;
  }
   
  printf(" done.\n");

  mpz_clear(bi0);
  mpz_clear(bi1);
  mpz_clear(p);
  mpz_clear(n);
  mpz_clear(sqrt_n);
  mpz_clear(i0);
  mpz_clear(i1);
  mpz_clear(q);
  mpz_clear(p_pow);
}
