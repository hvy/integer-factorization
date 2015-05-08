#include "quadratic_sieve.h"

#include <stdlib.h>
#include <limits.h>
#include <math.h>

#define BATCH_SIZE  100000000

struct mpz_t_pair {
  mpz_t x, y;  
};

/* Computes the smoothness bound B, given n. The smoothness bound B will be 
   the size of the factor base. The algorithm is inspired by Torbjorn 
   Granlund */
long smoothness_bound(const mpz_t n) {
  double nd = mpz_get_d(n);
  double b = 2 * exp(0.5 * sqrt(log(nd) * log(log(nd))));
  return (long) b;
}

/* Computes the sieving interval M, given b. The algorithm is inspired 
   by Eric Landquist. The quadratic sieve factoring algorithm, 2001. */
long sieving_interval(const long b) {
  return 4 * b * b;  
}

void get_q(mpz_t q, const mpz_t x, const mpz_t n, const mpz_t sqrt_n) {
  mpz_add(q, x, sqrt_n);
  mpz_pow_ui(q, q, 2);
  mpz_sub(q, q, n);
}

void q_prim(mpz_t q, const mpz_t x, const mpz_t sqrt_n) {
  mpz_add(q, x, sqrt_n);
  mpz_mul_ui(q, q, 2);
}

void q_prim_inverse(mpz_t q, const mpz_t x, const mpz_t sqrt_n, const mpz_t p) {
  q_prim(q, x, sqrt_n);
  int successfull = mpz_invert(q, q, p);
  if(0 == successfull /* Error, modular inverse didn't exist */) {
    printf("[ERROR] Failed in q_prim_inverse");
    mpz_clear(q);
  }
}

int factor_base(long *c, mpz_t *v, long b, const mpz_t n, const int primec, mpz_t *primev) {

  mpz_t b_mpz;
  mpz_init_set_si(b_mpz, b);
 
  /* If b is larger than the largest precomputed prime, then we won't have
     enough primes up until b. If that is the case, return an error code. */ 
  if(0 < mpz_cmp(b_mpz, primev[primec - 1])) {
    printf("[ERROR] Could not generate a factor base. B was too large.");
    return -1; // TODO Change this to a typedef enum error code
  }

  *c = 0;
  
  for(int i = 0; i < primec && 0 < mpz_cmp(b_mpz, primev[i]); ++i) {
    
    /* Compute the Legendre Symbol. Note that primev[i] needs to be an odd
       positive prime. */
    if(0 != mpz_cmp_ui(primev[i], 2) /* primev[i] is not 2, hence an odd prime */) {
      
      int legenre_symbol = mpz_legendre(n, primev[i]);
      if (1 == legenre_symbol) {
        mpz_init_set(v[*c], primev[i]);
        ++(*c);
      }
    }  
  }

  mpz_clear(b_mpz);
    
  return 0;
}

/* Tonelli-Shanks is used for solving a congruence on form x^2 = n (mod p). 
   It is used to find the start indices when finding the smooth numbers.
   @param x a start index (output)
   @param p an odd prime
   @prime n an integer which is a quadratic residue (mod p), meaning that
            the Legendre symbol (n/p) = 1 */ 
void tonelli_shanks(mpz_t x, const mpz_t p, const mpz_t n) {

  mpz_t two, minus_one, p_minus_one, q, s, z, c, r, t, m, tmp, tmp2, t_sqrt, b;

  mpz_init_set_ui(two, 2);
  mpz_init_set_ui(minus_one, -1);
  mpz_init(p_minus_one);
  mpz_init(q);
  mpz_init_set_ui(s, 0);
  mpz_init_set_ui(z, 2);
  mpz_init(c);
  mpz_init(r);
  mpz_init(t);
  mpz_init(m);
  mpz_init(tmp);
  mpz_init(tmp2);
  mpz_init(t_sqrt);
  mpz_init(b);
  
  mpz_sub_ui(p_minus_one, p, 1);
  mpz_set(q, p_minus_one);
  while(0 != mpz_even_p(q) /* while q is even */) {
    mpz_add_ui(s, s, 1);
    mpz_divexact_ui(q, q, 2);
  } 

  /* Find the quadratic non-residue */ 
  while(-1 != mpz_legendre(z, p)) {
    mpz_add_ui(z, z, 1);
  }
  
  mpz_powm(c, z, q, p);
  mpz_add_ui(tmp, q, 1); /* q is odd, so tmp will be even */
  mpz_divexact_ui(tmp, tmp, 2);
  mpz_powm(r, n, tmp, p);
  mpz_powm(t, n, q, p);
  mpz_set(m, s);

  while(1) {
    if(0 == mpz_cmp_ui(t, 1) /* t == 1 */) {
      mpz_set(x, r);
      break;
    } else {
      long msl = mpz_get_si(m);
      for(long i = 0; i < msl; ++i) {
        mpz_pow_ui(tmp, two, i);
        mpz_powm(t_sqrt, t, tmp, p);
        if(0 == mpz_cmp_ui(t_sqrt, 1) /* t^(1/2) == 1 */) {
          mpz_pow_ui(tmp, two, msl - i - 1);
          mpz_powm(b, c, tmp, p);
          mpz_mul(r, r, b);
          mpz_mod(r, r, p);
          mpz_pow_ui(tmp, b, 2);
          mpz_mul(t, t, tmp);
          mpz_mod(t, t, p);
          mpz_mod(c, tmp, p);
          mpz_set_ui(m, i);
          break;
        }    
      }
    }
  }

  mpz_clear(two);
  mpz_clear(minus_one);
  mpz_clear(p_minus_one);
  mpz_clear(q);
  mpz_clear(s);
  mpz_clear(z);
  mpz_clear(c);
  mpz_clear(r);
  mpz_clear(t);
  mpz_clear(m);
  mpz_clear(tmp);
  mpz_clear(tmp2);
  mpz_clear(t_sqrt);
  mpz_clear(b);
}

void hensel_lift(mpz_t x, const mpz_t r, const mpz_t p, const int k, const mpz_t n, const mpz_t n_floor_sqrt) {

  mpz_t p_pow_k, f, f_prim_inverse, t, tmp;

  mpz_init(p_pow_k);
  mpz_init(f);
  mpz_init(f_prim_inverse);
  mpz_init(t);
  mpz_init(tmp);

  mpz_pow_ui(p_pow_k, p, k);
  get_q(f, r, n, n_floor_sqrt);
  q_prim_inverse(f_prim_inverse, r, n_floor_sqrt, p_pow_k);
  mpz_mul(t, f, f_prim_inverse);
  mpz_cdiv_q(t, t, p_pow_k); /* TODO use mpz_divexact if we know in advance that p_pow_k divides t */
  mpz_mul_si(t, t, -1);
  mpz_mul(tmp, t, p_pow_k);
  mpz_add(x, r, tmp);
  mpz_pow_ui(tmp, p, (k + 1));
  mpz_mod(x, x, tmp); 

  mpz_clear(p_pow_k);
  mpz_clear(f);
  mpz_clear(f_prim_inverse);
  mpz_clear(t);
  mpz_clear(tmp);
}

int smooth_numbers(long *smooth_numberscp, mpz_t *smooth_numbersv, mpz_t *factorv,
  const long factor_basec, const mpz_t *factor_basev, const long m, 
  const mpz_t n) {
  
  *smooth_numberscp = 0;

  int int_max = INT_MAX; 
  int batch_size = BATCH_SIZE;
  int max_factor_pow = 8;
  float *logsv = (float *) malloc(batch_size * sizeof(float));
  float current_log = 0;
  mpz_t smooth_num, n_floored_sqrt; // TODO optimize speed by implementing ceiling sqrt
  mpz_t mpz_offset, mpz_batch_size, offset_batch_size, q, q_offset_tmp, q_offset;
  mpz_t p_pow;
  mpz_t i0, i1, bi0, bi1, d_i;
  mpz_t tmp;
  mpz_t offset_diff, offset_multiples, i0_offset;
  mpz_t x;

  mpz_init(smooth_num);
  mpz_init(i0);
  mpz_init(i1);
  mpz_init(bi0);
  mpz_init(bi1);
  mpz_init(q);
  mpz_init(q_offset);
  mpz_init(q_offset_tmp);
  mpz_init(mpz_offset);
  mpz_init(mpz_batch_size);
  mpz_init(offset_batch_size);
  mpz_set_ui(mpz_batch_size, batch_size);
  mpz_init(n_floored_sqrt);
  mpz_sqrt(n_floored_sqrt, n);
  mpz_init(p_pow);
  mpz_init(d_i);
  mpz_init(tmp);
  mpz_init(offset_diff);
  mpz_init(offset_multiples);
  mpz_init(i0_offset);
  mpz_init(x);

  for(long offset = 0; offset * batch_size < m && *smooth_numberscp < factor_basec; ++offset) {
    mpz_set_si(mpz_offset, offset);
    mpz_mul(offset_batch_size, mpz_offset, mpz_batch_size);

    printf("[INFO] Offset: %ld Smooth numbers: %ld\n", offset, *smooth_numberscp);
    fflush(stdout);

    /* Initialize each element to the same value since the logs of the big 
       roots are very similar in size. */
    long q_offset_ui = (offset + 1) * batch_size;
    mpz_set_ui(q_offset_tmp, q_offset_ui);
    get_q(q_offset, q_offset_tmp, n, n_floored_sqrt);

    current_log = (float) log(fabs(mpz_get_d(q_offset)));

    for(int i = 0; i < batch_size; ++i) {
      logsv[i] = current_log;  
    }

    for(long i = 0; i < factor_basec; ++i) {
      double log_p = log(mpz_get_d(factor_basev[i]));

      /* Skip small number since they won't contribute to the sieving */
      if(0 == mpz_cmp_ui(factor_basev[i], 2) ||
          0 == mpz_cmp_ui(factor_basev[i], 3) ||
          0 == mpz_cmp_ui(factor_basev[i], 5)) {
        continue;
      }   
       
      mpz_set_ui(i0, 0);
      mpz_set_ui(i1, 0);
      mpz_set_ui(bi0, 0);
      mpz_set_ui(bi1, 0);    
    
      tonelli_shanks(bi0, factor_basev[i], n);
      mpz_sub(bi1, factor_basev[i], bi0);
      mpz_sub(i0, bi0, n_floored_sqrt);
      mpz_sub(i1, bi1, n_floored_sqrt);
      mpz_mod(i0, i0, factor_basev[i]);
      mpz_mod(i1, i1, factor_basev[i]);

      for(int t = 1; t <= max_factor_pow; ++t) {
        mpz_pow_ui(p_pow, factor_basev[i], t);
        int p_powi = -1;
        if(0 >= mpz_cmp_ui(p_pow, int_max) /* p_pow < INT_MAX */) {
          p_powi = mpz_get_ui(p_pow);
        }
        
        /* Higher powers than 1 requires Hensel Lift */
        if(t > 1) {
          hensel_lift(bi0, i0, factor_basev[i], t - 1, n, n_floored_sqrt);
          hensel_lift(bi1, i1, factor_basev[i], t - 1, n, n_floored_sqrt);
          mpz_set(i0, bi0);
          mpz_set(i1, bi1);
        }

        /* Compute the difference between the indices once, then we only need
           the first index. */
        mpz_sub(d_i, i1, i0);
        int d_ii = mpz_get_si(d_i);      
        if(d_ii < 0) {
          mpz_set(tmp, i0);
          mpz_set(i0, i1);
          mpz_set(i1, tmp);
          d_ii *= -1;
        }           
        
        /* Compute the new offset. */
        mpz_sub(offset_diff, offset_batch_size, i0);
        mpz_cdiv_q(offset_multiples, offset_diff, p_pow);
        mpz_mod(tmp, offset_diff, p_pow);
        if(0 != mpz_cmp_ui(tmp, 0)) {
          mpz_add_ui(offset_multiples, offset_multiples, 1);
        }

        /* Compute the new i0 with the offset */
        mpz_mul(i0_offset, offset_multiples, p_pow);
        mpz_add(i0_offset, i0_offset, i0);

        mpz_add_ui(tmp, offset_batch_size, batch_size);
        if(0 <= mpz_cmp(i0_offset, offset_batch_size) &&
            0 >= mpz_cmp(i0_offset, tmp) /* inside offset range */) {
          mpz_mod_ui(i0_offset, i0_offset, batch_size);
        
          if(p_powi > 0 /* p_pow <= INT_MAX */) {
            int j;
            for(j = mpz_get_si(i0_offset) + d_ii; j < batch_size &&
              j >= 0; j += p_powi) {
              logsv[j - d_ii] -= log_p;
              logsv[j] -= log_p;
            }

            /* Make sure that the lower index is inside the range */
            if((j - d_ii) < batch_size && (j - d_ii) >= 0) {
              logsv[j - d_ii] -= log_p;
            }
          } else /* p_pow is too big, we only add it once */ {
            int i0_offseti = mpz_get_si(i0_offset);
            logsv[i0_offseti] -= log_p;
            if((i0_offseti + d_ii) < batch_size && (i0_offseti + d_ii) >= 0) {
              logsv[i0_offseti + d_ii] -= log_p;
            }
          }
        }
      }       
    }
  
    long potential = 0;
    long actual = 0;
    double biggest_root_log = log(mpz_get_si(factor_basev[factor_basec - 1]));
    for(int i = 0; i < batch_size && *smooth_numberscp <= factor_basec; ++i) {
      if(fabs(logsv[i]) <= (biggest_root_log + 0.5) 
        /* potiential root, do trial division */) {
        ++potential;
        mpz_set_ui(x, offset * batch_size + i);
        get_q(q, x, n, n_floored_sqrt);
        for(long factor_base_i = 0; factor_base_i < factor_basec; 
          ++factor_base_i) {
          mpz_mod(tmp, q, factor_basev[factor_base_i]);
          if(0 == mpz_cmp_ui(tmp, 0) /* divisible by this root */) {
            do {
              mpz_cdiv_q(q, q, factor_basev[factor_base_i]);
              mpz_mod(tmp, q, factor_basev[factor_base_i]);
            } while (0 == mpz_cmp_ui(tmp, 0));
            if(0 == mpz_cmp_ui(q, 1)) {
              ++actual;
              //mpz_add(smooth_numbersv[*smooth_numberscp], x, n_floored_sqrt);
              // TODO Have another same length vector that corresponds to the pair.
              // TODO It needs to be assed as an argument and malloced.
              mpz_add(factorv[*smooth_numberscp], x, n_floored_sqrt);
              get_q(smooth_num, x, n, n_floored_sqrt);
              mpz_set(smooth_numbersv[*smooth_numberscp], smooth_num);
              ++(*smooth_numberscp);
              break;
            }
          } 
        }
      } 
    }

    printf("[INFO] Smooth numbers found in batch: %ld / %ld\n", actual, potential);
  }
  
  mpz_clear(d_i);

  free(logsv);
    
  return 0;
}


int get_left_null_space(mpz_t *left_null_space, const int smooth_numberc, mpz_t *factor_expv) {

  int left_null_spacec = 0;
  mpz_t id[smooth_numberc];

  /* Create the identity matrix first to create the left null space */
  for(int i = 0; i < smooth_numberc; ++i) {
    mpz_init2(id[i], smooth_numberc);
    mpz_setbit(id[i], smooth_numberc - 1 - i);
  }
 
  for(int ri = 0; ri < smooth_numberc; ++ri) {
    mp_bitcnt_t pivot = mpz_scan1(factor_expv[ri], 0);
    if(pivot < smooth_numberc /* at least one bit is set to 1 */) {
      for(int ri_inner = ri + 1; ri_inner < smooth_numberc; ++ri_inner) {
        if(1 == mpz_tstbit(factor_expv[ri_inner], pivot)) {
          mpz_xor(factor_expv[ri_inner], factor_expv[ri], factor_expv[ri_inner]);    
          mpz_xor(id[ri_inner], id[ri], id[ri_inner]);    
        }
      }    
    }
  } 

  for(int ri = 0; ri < smooth_numberc; ++ri) {
    if(0 == mpz_popcount(factor_expv[ri])) {
      mpz_set(left_null_space[left_null_spacec], id[ri]);
      ++left_null_spacec; 
    } 
  }

  return left_null_spacec;
}

void gaussian_elimination(mpz_t factor, const mpz_t n, int smooth_numberc, 
  mpz_t const *smooth_numberv, mpz_t const *factorv, const int factor_basec, mpz_t *factor_expv) {
 
  int left_null_spacec;
  mpz_t left_null_spacev[smooth_numberc];
  mpz_t prod_x_plus_sq, prod_smooth;
  mpz_t factor_candidate, tmp1, tmp2;

  mpz_init(prod_x_plus_sq);
  mpz_init(prod_smooth);
  mpz_init(factor_candidate);
  mpz_init(tmp1);
  mpz_init(tmp2);

  left_null_spacec = get_left_null_space(left_null_spacev, smooth_numberc, factor_expv);

  printf("Size of left null space: %d\n", left_null_spacec);

  for(int i = 0; i < left_null_spacec; ++i) {
 
 
    mpz_set_ui(prod_x_plus_sq, 0);
    mpz_set_ui(prod_smooth, 0);

    size_t num_bits = mpz_sizeinbase(left_null_spacev[i], 2);
    printf("Num of bits: %zu\n", num_bits);
    printf("Smooth number c %d\n", smooth_numberc);
    for(int j = 0; j < num_bits; ++j) {
      if(1 == mpz_tstbit(left_null_spacev[i], smooth_numberc - 1 - j)) { /* reversed bit index */
        if((0 == mpz_cmp_ui(prod_x_plus_sq, 0)) || (0 == mpz_cmp_ui(prod_smooth, 0))) {
          mpz_set(prod_x_plus_sq, factorv[j]);
          mpz_set(prod_smooth, smooth_numberv[j]);
        } else {
          mpz_mul(prod_x_plus_sq, prod_x_plus_sq, factorv[j]);
          mpz_mul(prod_smooth, prod_smooth, smooth_numberv[j]);
        } 
      }  
    }
    
    /* Now we have a candidate for the factor. Check if it is non trivial  */
    mpz_sqrt(tmp1, prod_smooth);
    mpz_sub(tmp2, prod_x_plus_sq, tmp1);
    mpz_gcd(factor_candidate, tmp2, n);
    if((0 != mpz_cmp_ui(factor_candidate, 0)) && (0 != mpz_cmp(factor_candidate, n))) {
      printf("[INFO] Found factor!");
      mpz_set(factor, factor_candidate);
      break;
    } else {
      mpz_add(tmp2, prod_x_plus_sq, tmp1);
      mpz_gcd(factor_candidate, tmp2, n);
      if((0 != mpz_cmp_ui(factor_candidate, 0)) && (0 != mpz_cmp(factor_candidate, n))) {
        printf("[INFO] Found factor!");
        mpz_set(factor, factor_candidate);
        break;
      }
    } 
  }

  mpz_clear(prod_x_plus_sq);
  mpz_clear(prod_smooth);
  mpz_clear(factor_candidate);
  mpz_clear(tmp1);
  mpz_clear(tmp2);
}

void quadratic_sieve(mpz_t factor, mpz_t n, const int primec, mpz_t *primev) {

  mpz_init(factor);

  mpz_t smooth_num, tmp;
  mpz_init(smooth_num); 
  mpz_init(tmp);
 
  if (0 != mpz_perfect_square_p(n) /* n is a perfect square */) {
    mpz_sqrt(factor, n);
    return;
  }

  long b = smoothness_bound(n);
  long m = sieving_interval(b);
 
  printf("[INFO] Smoothness bound, B: %ld\n", b);
  printf("[INFO] Sieving interval, M: %ld\n", m);
  
  long factor_basec = 0;
  mpz_t *factor_basev = (mpz_t *) malloc(b * sizeof(mpz_t));

  printf("[INFO] Computing factor base...");
  fflush(stdout);
  
  int factor_base_result = factor_base(&factor_basec, factor_basev, b, n, primec, primev); 

  if(0 == factor_base_result) {  
    printf(" done.\n");
    printf("[INFO] Size of factor base: %ld\n", factor_basec);
  } else {
    printf(" failed.\n");
  }

  long smooth_numberc = 0;
  mpz_t *smooth_numberv = (mpz_t *) malloc(b * sizeof(mpz_t));
  mpz_t *factorv = (mpz_t *) malloc(b * sizeof(mpz_t));

  printf("[INFO] Computing smooth numbers...");
  fflush(stdout);

  int smooth_number_result = smooth_numbers(&smooth_numberc, smooth_numberv, factorv, factor_basec, factor_basev, m, n);

  if(0 == smooth_number_result) {
    printf(" done.\n");
    printf("[INFO] Number of smooth numbers: %ld\n", smooth_numberc);
  } else {
    printf(" failed.\n");
  }

  /* make sure that the array allocation works */
  /*
  mpz_t factor_expv[smooth_numberc];
  int powc;
  for(int i = 0; i < smooth_numberc; ++i) {
    mpz_init2(factor_expv[i], factor_basec);
    mpz_set(smooth_num, smooth_numberv[i]);
    for(int j = 0; j < factor_basec; ++j) {
      powc = 0;
      mpz_mod(tmp, smooth_num, factor_basev[j]);
      while(0 == mpz_cmp_ui(tmp, 0)) {
        mpz_divexact(smooth_num, smooth_num, factor_basev[j]);
        ++powc; 
        mpz_mod(tmp, smooth_num, factor_basev[j]);
      }
      if(1 == (powc % 2)) {
        mpz_setbit(factor_expv[i], j);
      } 
    }
  } 
  */
  mpz_t factor_expv[smooth_numberc];
  int powc;
  for(int i = 0; i < smooth_numberc; ++i) {
    mpz_init2(factor_expv[i], factor_basec);
    mpz_set(smooth_num, smooth_numberv[i]);
    for(int j = 0; j < factor_basec; ++j) {
      powc = 0;
      mpz_mod(tmp, smooth_num, factor_basev[j]);
      while(0 == mpz_cmp_ui(tmp, 0)) {
        mpz_divexact(smooth_num, smooth_num, factor_basev[j]);
        ++powc; 
        mpz_mod(tmp, smooth_num, factor_basev[j]);
      }
      if(1 == (powc % 2)) {
        mpz_setbit(factor_expv[i], j);
      } 
    }
  } 
 
  gaussian_elimination(factor, n, smooth_numberc, smooth_numberv, factorv, factor_basec, factor_expv);

  free(factor_basev);
  free(smooth_numberv);
  free(factorv);
}
