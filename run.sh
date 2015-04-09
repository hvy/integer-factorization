#!/bin/bash
gcc -o intfact src/main.c src/trial_division.c src/pollards_rho.c src/quadratic_sieve.c -lgmp
./intfact
rm intfact
