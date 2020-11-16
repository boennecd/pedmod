/* $Id: C_FORTRAN_interface.c 313 2015-09-16 20:20:04Z mmaechler $
*
*  wrapper for calling R's random number generator from
*  the original FORTRAN code
*
*/

#include <R.h>
#include "pnorm.h"

double F77_SUB(mvphi)(double const *z){
  return pnorm_std(*z, 1L, 0L);
}
