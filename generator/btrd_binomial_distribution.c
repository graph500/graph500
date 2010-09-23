/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include "btrd_binomial_distribution.h"
#include <math.h>
#include <stddef.h>
#include "splittable_mrg.h"

static double f_c(unsigned int k) {
  static const double one_over_12 = 1./12;
  static const double one_over_360 = 1./360;
  static const double one_over_1260 = 1./1260;
  static const double values[] = { /* Constants for f_c(0) to f_c(9) */
    .08106146679532726,
    .04134069595540929,
    .02767792568499834,
    .02079067210376509,
    .01664469118982119,
    .01387612882307075,
    .01189670994589177,
    .01041126526197209,
    .009255462182712733,
    .008330563433362871
  };
  if (k <= 9) {
    return values[k];
  } else {
    double recip_k_plus_1 = 1. / (k + 1);
    return (one_over_12 - (one_over_360 - one_over_1260 * recip_k_plus_1 * recip_k_plus_1) * recip_k_plus_1 * recip_k_plus_1) * recip_k_plus_1;
  }
}

#ifdef __MTA__
#pragma mta expect parallel context
#pragma mta serial
#endif
size_t btrd_binomial_distribution(size_t n_orig, double p, mrg_state* state) {
  /* BTRD algorithm from pages 6--7 of "The Generation of Binomial Random
   * Variates" (Wolfgang Hoermann) --
   * http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.47.8407 */
  if (p == 0.) return 0;
  if (p > .5) return n_orig - btrd_binomial_distribution(n_orig, 1. - p, state);
  if (n_orig * p < 10) {
    /* First waiting time algorithm from page 525 of
     * http://cg.scs.carleton.ca/~luc/chapter_ten.pdf */
    int x = 0;
    int sum = 0;
    double recip_log_1_minus_p;
    /* Approximation to ln(1-p) from second sequence in exercise 1 of
     * http://cg.scs.carleton.ca/~luc/chapter_ten.pdf page 500; done to avoid
     * use of C99 log1p function. */
    {
      double r = 1. - 2. / p;
      double recip_r = 1. / r;
      double recip_r2 = recip_r * recip_r;
      double recip_r4 = recip_r2 * recip_r2;
      double recip_r6 = recip_r4 * recip_r2;
      double log_1_minus_p = 2. * recip_r * (1 + recip_r2 / 3. + recip_r4 / 5. + recip_r6 / 7.);
      recip_log_1_minus_p = 1. / log_1_minus_p;
    }
#if 0
    /* Old version using C99 log1p: */
    recip_log_1_minus_p = 1. / log1p(-p); /* 1 / ln(1 - p) */
#endif
    do {
      /* From page 500 of http://cg.scs.carleton.ca/~luc/chapter_ten.pdf
       * (geometric variate generator). */
      sum += (int)ceil(log(mrg_get_double_orig(state)) * recip_log_1_minus_p);
      ++x;
    } while (sum <= n_orig);
    return x - 1;
  }
  if (n_orig > 1000000000) {
    return btrd_binomial_distribution(1000000000, p, state) +
           btrd_binomial_distribution(n_orig - 1000000000, p, state);
  }
  {
    int n = (int)n_orig;
    /* Setup */
    int m = (int)floor((n + 1) * p);
    double r = p / (1. - p);
    double nr = (n + 1) * r;
    double npq = n * p * (1. - p);
    double sqrt_npq = sqrt(npq);
    double b = 1.15 + 2.53 * sqrt_npq;
    double a = -.0873 + .0248 * b + .01 * p;
    double c = n * p + .5;
    double alpha = (2.83 + 5.1 / b) * sqrt_npq;
    double v_r = .92 - 4.2 / b;
    double u_rv_r = .86 * v_r;
    while (1) {
      /* 1 */
      double v = mrg_get_double_orig(state);
      double u, us, temp_random_num;
      int k, km;
      if (v <= u_rv_r) {
        u = v / v_r - .43;
        return (int)floor((2 * a / (.5 + fabs(u)) + b) * u + c);
      }
      /* 2 */
      temp_random_num = mrg_get_double_orig(state);
      if (v >= v_r) {
        u = temp_random_num - .5;
      } else {
        u = v / v_r - .93;
        u = (u > 0. ? 1. : u < 0. ? -1. : 0.) * .5 - u;
        v = v_r * temp_random_num;
      }
      /* 3.0 */
      us = .5 - fabs(u);
      k = (int)floor((2 * a / us + b) * u + c);
      if (k < 0 || k > n) continue;
      v *= alpha / (a / (us * us) + b);
      km = (k >= m ? k - m : m - k);
      if (km > 15) {
        /* 3.2 */
        double rho = (km / npq) * (((km / 3 + .625) * km + 1. / 6) / npq + .5);
        double t = -km * km / (2 * npq);
        int nm, nk;
        double h, threshold;
        v = log(v);
        if (v < t - rho) return k;
        if (v > t + rho) continue;
        /* 3.3 */
        nm = n - m + 1;
        h = (m + .5) * log((m + 1) / (r * nm)) + f_c(m) + f_c(n - m);
        /* 3.4 */
        nk = n - k + 1;
        threshold = h +
                    (n + 1) * log((double)nm / nk) +
                    (k + .5) * log(nk * r / (k + 1)) -
                    f_c(k) -
                    f_c(n - k);
        if (v <= threshold) return k;
      } else {
        /* 3.1 */
        double f = 1.;
        int i;
        if (m < k) {
          for (i = m; i != k; ++i) f *= nr / i - r;
        } else if (m > k) {
          for (i = k; i != m; ++i) v *= nr / i + r;
        }
        if (v <= f) return k;
      }
    }
  }
}
