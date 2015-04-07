/* Pollard's rho algorithm implementation in C using the GMP library to find a 
   prime factor given a number. Brent's cycle finding method is used to speed 
   up the algorithm. 
   
   Note: When compiling, specify the path to gmp.h using the -I flag if it is 
   not automatically found by the compiler and make sure to compile with the 
   -lgmp flag as well. E.g. gcc -I/usr/local/include/ main.cpp -lgmp */

#include "trial_division.h"
#include "pollards_rho.h"

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#define BASE                  10
#define TRIAL_DIVISION_PRIMES 1000000

gmp_randstate_t rnd_state;

int main(int argc, char *argv[]) {
  
  int trial_division_primec;
  mpz_t *trial_division_primev;  
  mpz_t n, factor;
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
  printf("\n\n");

  /* Set the factor to n if n == 1 or if n is most likely a prime, else use 
     Pollard's rho algorithm to find a prime factor */
  if(0 == mpz_cmp_ui(n, 1) || mpz_probab_prime_p(n, 25)) {
    mpz_set(factor, n);
  } else {
    /* Pollard's Rho algorithm */
    //pollard_rho(factor, n, rnd_state);
    
    /* Trial division */
    trial_division_primec = TRIAL_DIVISION_PRIMES;
    trial_division_primev = (mpz_t *) malloc(trial_division_primec * sizeof(mpz_t));

    printf("\tComputing the %d first primes...", trial_division_primec);
    fflush(stdout);
    get_fst_primes(trial_division_primec, trial_division_primev);
    printf(" done (%luMB).\n\n", trial_division_primec * sizeof(mpz_t) / 1000000);

    trial_division(factor, n, trial_division_primec, trial_division_primev);

    /* Quadratic Sieve */
  }  
  
  printf("\tFactor: ");
  mpz_out_str(stdout, BASE, factor);
  printf("\n");

  free(trial_division_primev);

  return 0;
}
