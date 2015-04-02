/* Pollard's rho algorithm implementation in C using the GMP library to find a 
   prime factor given a number. Brent's cycle finding method is used to speed 
   up the algorithm. 
   
   Note: When compiling, specify the path to gmp.h using the -I flag if it is 
   not automatically found by the compiler and make sure to compile with the 
   -lgmp flag as well. E.g. gcc -I/usr/local/include/ main.cpp -lgmp */

#include <stdio.h>
#include <assert.h>
#include <gmp.h>

#define BASE 10

mpz_t one, two;
gmp_randstate_t rnd_state;

void pollard_rho(mpz_t factor, mpz_t n);

int main(int argc, char *argv[]) {
    
  mpz_t n, factor;
  mpz_init(factor);
  mpz_init_set_ui(one, 1);
  mpz_init_set_ui(two, 2);
  gmp_randinit_default(rnd_state);

  if(argc < 2) {
    printf("No argument provided, will factorize build in number\n");  
    mpz_init_set_str(n, "502560280658509", BASE);
  } else {
    mpz_init_set_str(n, argv[1], BASE);
  }

  /* Print the number that will be factorized */
  printf("\tFactorizing: ");
  mpz_out_str(stdout, BASE, n);
  printf("\n");

  /* Set the factor to n if n == 1 or if n is most likely a prime, else use 
     Pollard's rho algorithm to find a prime factor */
  if(0 == mpz_cmp(n, one) || mpz_probab_prime_p(n, 25)) {
    mpz_set(factor, n);
  } else {
    pollard_rho(factor, n);
  }  
  
  printf("\tFactor: ");
  mpz_out_str(stdout, BASE, factor);
  printf("\n");

  return 0;
}


void pollard_rho(mpz_t factor, mpz_t n) {

  /* Assert n != 1 */
  assert(0 != mpz_cmp(n, one));

  mpz_t n_minus_1, a, xi, xm, s, tmp_abs, tmp_gcd;

  mpz_init(xi);
  mpz_init(s);
  mpz_init(tmp_abs);
  mpz_init(tmp_gcd);
  mpz_init_set(n_minus_1, n);
  mpz_sub_ui(n_minus_1, n_minus_1, 1);
  mpz_init_set_ui(a, 1);
  mpz_urandomm(xi, rnd_state, n_minus_1);
  mpz_init_set(xm, xi);

  for(int i = 0; i < 10000; ++i) {
    mpz_mul(xi, xi, xi); 
    mpz_add(xi, xi, a); 
    mpz_mod(xi, xi, n);
    mpz_sub(tmp_abs, xi, xm);
    mpz_abs(tmp_abs, tmp_abs);
    mpz_gcd(tmp_gcd, tmp_abs, n);
    mpz_set(s, tmp_gcd);

    if(0 != mpz_cmp(s, one) && 0 != mpz_cmp(s, n) /* found a factor */) {
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
  mpz_clear(tmp_gcd); 
}
