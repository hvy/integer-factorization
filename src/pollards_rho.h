#ifndef POLLARDS_RHO_H
#define POLLARDS_RHO_H

#include <stdio.h>
#include <gmp.h>

void pollard_rho(mpz_t factor, mpz_t n, gmp_randstate_t rnd_state);

#endif
