#!/bin/bash
gcc -o test src/test.c src/trial_division.c src/pollards_rho.c src/quadratic_sieve.c -lgmp
./test
rm test
