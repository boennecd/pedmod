#ifndef PNORM_H
#define PNORM_H

#include <Rmath.h>
// I ran into issues with likes those mentioned here
//   https://bugs.webkit.org/show_bug.cgi?id=59249
// Thus, we always use math.h
#ifdef beta
// we get an error if we do not undefine beta
#undef beta
#endif
#include <math.h> // may be included anyway?

/**
 * evaluates the standard normal CDF after avoiding some checks in the
 * R function.
 */
static inline double pnorm_std(double const x, int lower, int is_log) {
  int const is_inf = isinf(x),
            is_nan = isnan(x);

  if(!is_inf && !is_nan){
    double p, cp;
    p = x;
    Rf_pnorm_both(x, &p, &cp, lower ? 0 : 1, is_log);
    return lower ? p : cp;

  } else if(is_inf){
    if(is_log){
      if(lower)
        return x > 0 ? 0 : -INFINITY;
      else
        return x < 0 ? 0 :  INFINITY;
    } else {
      if(lower)
        return x > 0 ? 1 : 0;
      else
        return x < 0 ? 1 : 0;
    }
  }

  // isnan == TRUE
  return NAN;
}

/**
 * evaluates the normal CDF after avoiding some checks in the R function.
 */
static inline double pnorm_w(double const x, double const mu, double const sigma,
                             int lower, int is_log) {
  return pnorm_std((x - mu) / sigma, lower, is_log);
}

#endif
