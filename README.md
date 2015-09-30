# Integer Factorization

An integer factorizer implmented in C using the GNU Multi-Precision Arithmetic Library (GMP, https://gmplib.org). The factorizer combines trial division, Pollard's Rho algorithm and the Quadratic Sieve for optimal results.

## How to Run
1. Install the GMP library from https://gmplib.org.
2. Run the Bash script *run.sh*.
3. Enter the number to factorize and press *Enter*.
4. Wait and hope for the best...

## List of algorithms

- Trial Division - A set of smaller prime numbers are precomputed when the factorizer is started. If any of the primes divides the given number, we have found a prime factor.

- The Pollard's Rho Algorithm - The Pollard's Rho algorithm with Brent's cycle detection is used to quickly find smaller prime factors that were not found using the trial division.

- Quadratc Sieve - For larger integers that is not divided by any precomputed prime or that the Pollard's Rho algorithm can't handle, the Quadratic Sieve is used. This implementation of the Quadratic Sieve includes the Tonelli-Shanks algorithm to find square roots module primes and the Hensel's lifting lemma to speed up the factorization by generating the exponent vector for higher power larger than 1.
