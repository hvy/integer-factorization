#include "trial_division.h"
#include "pollards_rho.h"
#include "quadratic_sieve.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>

void test_tonelli_shanks();
void test_hensel_lift();
void test_left_null_space();
void test_smooth_numbers();

int main() {
  
  printf("[INFO] Starting test... \n");

  test_tonelli_shanks();
  test_hensel_lift();
  test_left_null_space();
  //test_smooth_numbers();
  
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


void test_left_null_space() {

  printf("[INFO] Starting Left Null Space...\n");

  int max = 5;
  int left_null_spacec;
  mpz_t left_null_spacev[5];
  mpz_t roots[5];
  mpz_t roots_unmod[5];
  mpz_t tmp;

  mpz_init(tmp);

  for(int i = 0; i < 5; ++i) {
    mpz_init2(roots[i], 5);
    mpz_init2(roots_unmod[i], 5);
  }
 
  mpz_setbit(roots[0], 0); 
  mpz_setbit(roots[0], 1); 
  mpz_setbit(roots[1], 2); 
  mpz_setbit(roots[1], 4); 
  mpz_setbit(roots[2], 1); 
  mpz_setbit(roots[2], 3); 
  mpz_setbit(roots[3], 0); 
  mpz_setbit(roots[3], 3); 
  mpz_setbit(roots[4], 0); 
  mpz_setbit(roots[4], 1); 
  
  for(int i = 0; i < 5; ++i) {
    mpz_set(roots_unmod[i], roots[i]);
  }

  left_null_spacec = get_left_null_space(left_null_spacev, max, roots);

  printf("Number of left null spaces: %d\n", left_null_spacec);

  for(int i = 0; i < left_null_spacec; ++i) {
    printf("{");
    mpz_out_str(stdout, 2, left_null_spacev[i]);
    printf("}\n");
  
    mpz_set_ui(tmp, 0);

    for(int j = 0; j < 5; ++j) {
      if(1 == mpz_tstbit(left_null_spacev[i], 4 - j)) {
        mpz_xor(tmp, tmp, roots_unmod[j]);  
      }
    }
    
    assert(0 == mpz_popcount(tmp));  
  }

  mpz_clear(tmp);

  printf(" done.\n");
}

void test_smooth_numbers() {
  
  printf("[INFO] Starting Smooth Numbers...\n");

  int trial_division_primec;
  mpz_t n, *trial_division_primev, *factor_basev, *smooth_numberv;
  long b, m, factor_basec, smooth_numberc;
  char *n_str = "34567876234234234242342343242342342342345231432345";
 
  mpz_init_set_str(n, n_str, 10);
  b = smoothness_bound(n);
  m = sieving_interval(b);

  mpz_t *factorv = (mpz_t *) malloc(b * sizeof(mpz_t));
  
  printf("[INFO] Smoothness bound, B: %ld\n", b);                               
  printf("[INFO] Sieving interval, M: %ld\n", m);        

  /* Precompute the 1000000 first primes */
  trial_division_primec = 1000000;
  trial_division_primev = 
    (mpz_t *) malloc(trial_division_primec * sizeof(mpz_t));
  get_fst_primes(trial_division_primec, trial_division_primev);        
  
  /* Generate the factor base */
  factor_basec = 0;
  factor_basev = (mpz_t *) malloc(b * sizeof(mpz_t));
  int factor_base_result = factor_base(&factor_basec, factor_basev, b, n, 
    trial_division_primec, trial_division_primev);
  if(0 == factor_base_result) {                                                 
    printf(" done.\n");                                                         
    printf("[INFO] Size of factor base: %ld\n", factor_basec);                  
  } else {                                                                      
    printf(" failed.\n");                                                       
  }    
 
  /* Find the smooth numbers using the precomputed primes and 
   the factor base */ 
  smooth_numberc = 0;
  smooth_numberv = (mpz_t *) malloc(b * sizeof(mpz_t));
  smooth_numbers(&smooth_numberc, smooth_numberv, factorv, factor_basec, factor_basev,
    m, n);
  
  printf(" done.\n");
   
  for(int i = 0; i < smooth_numberc; ++i) {
    printf("Smooth number %d:\t%ld\n", (i + 1), mpz_get_si(smooth_numberv[i])); 
  }
  
  free(trial_division_primev);
  free(factor_basev);
  free(smooth_numberv);

  mpz_clear(n);
}
