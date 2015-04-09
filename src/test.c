#include "trial_division.h"
#include "pollards_rho.h"
#include "quadratic_sieve.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>

void test_tonelli_shanks();

int main() {
  
  printf("[INFO] Starting test... \n");

  test_tonelli_shanks();

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
