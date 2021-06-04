#ifndef NEW_MVT_H
#define NEW_MVT_H
#include <memory>
#include <algorithm>
#include <limits>
#include "threat-safe-random.h"
#include "ped-mem.h"
#include "kahan.h"
#include "sobol.h"
#include <vector>

namespace pedmod {
struct rand_Korobov_output {
  size_t minvls;
  double abserr;
  int inform;
};

constexpr int n_qmc_seqs() {
  return 64;
}

template<class Func>
class rand_Korobov {
  static cache_mem<double> dmem;
  static cache_mem<int   > imem;

public:
  static void alloc_mem
    (int const max_ndim, int const max_nf, int const max_threads) {
    dmem.set_n_mem(
      (5 + n_qmc_seqs()) * max_nf + (n_qmc_seqs() + 2) * max_ndim, max_threads);
    imem.set_n_mem(max_ndim                 , max_threads);
  }

  static rand_Korobov_output comp
    (Func &f, int const ndim, size_t const minvls, size_t const maxvls,
     int const nf, double const abseps, double const releps,
     double * const __restrict__ finest,
     double * const __restrict__ sdest, parallelrng::unif_drawer &sampler,
     unsigned const n_sequences){
    if(n_sequences < 2)
      throw std::invalid_argument("n_sequences is less than two");

    /* constants */
    constexpr int const plim(28L),
                        klim(100L);
    constexpr int const p[plim] = { 31L, 47L, 73L, 113L, 173L, 263L, 397L, 593L, 907L, 1361L, 2053L, 3079L, 4621L, 6947L, 10427L, 15641L, 23473L, 35221L, 52837L, 79259L, 118891L, 178349L, 267523L, 401287L, 601943L, 902933L, 1354471L, 2031713L };
    constexpr int const C[plim][klim - 1L] = {
      { 12L, 9L, 9L, 13L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 3L, 3L, 3L, 12L, 7L, 7L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 3L, 3L, 3L, 12L, 7L, 7L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 3L, 3L, 3L, 12L, 7L, 7L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 3L, 3L, 3L, 12L, 7L, 7L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 7L, 3L, 3L, 3L, 7L, 7L, 7L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L },
      { 13L, 11L, 17L, 10L, 15L, 15L, 15L, 15L, 15L, 15L, 22L, 15L, 15L, 6L, 6L, 6L, 15L, 15L, 9L, 13L, 2L, 2L, 2L, 13L, 11L, 11L, 10L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 6L, 6L, 6L, 15L, 15L, 9L, 13L, 2L, 2L, 2L, 13L, 11L, 11L, 10L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 6L, 6L, 6L, 15L, 15L, 9L, 13L, 2L, 2L, 2L, 13L, 11L, 11L, 10L, 10L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 6L, 2L, 3L, 2L, 3L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L },
      { 27L, 28L, 10L, 11L, 11L, 20L, 11L, 11L, 28L, 13L, 13L, 28L, 13L, 13L, 13L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 31L, 31L, 5L, 5L, 5L, 31L, 13L, 11L, 11L, 11L, 11L, 11L, 11L, 13L, 13L, 13L, 13L, 13L, 13L, 13L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 31L, 31L, 5L, 5L, 5L, 11L, 13L, 11L, 11L, 11L, 11L, 11L, 11L, 11L, 13L, 13L, 11L, 13L, 5L, 5L, 5L, 5L, 14L, 13L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L },
      { 35L, 27L, 27L, 36L, 22L, 29L, 29L, 20L, 45L, 5L, 5L, 5L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 29L, 17L, 17L, 17L, 17L, 17L, 17L, 17L, 17L, 17L, 17L, 23L, 23L, 23L, 23L, 23L, 23L, 23L, 23L, 23L, 23L, 23L, 23L, 21L, 27L, 3L, 3L, 3L, 24L, 27L, 27L, 17L, 29L, 29L, 29L, 17L, 5L, 5L, 5L, 5L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 17L, 17L, 17L, 6L, 17L, 17L, 6L, 3L, 6L, 6L, 3L, 3L, 3L, 3L, 3L },
      { 64L, 66L, 28L, 28L, 44L, 44L, 55L, 67L, 10L, 10L, 10L, 10L, 10L, 10L, 38L, 38L, 10L, 10L, 10L, 10L, 10L, 49L, 49L, 49L, 49L, 49L, 49L, 49L, 49L, 49L, 49L, 49L, 49L, 38L, 38L, 31L, 4L, 4L, 31L, 64L, 4L, 4L, 4L, 64L, 45L, 45L, 45L, 45L, 45L, 45L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 11L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 66L, 45L, 11L, 7L, 3L, 2L, 2L, 2L, 27L, 5L, 3L, 3L, 5L, 5L, 2L, 2L, 2L, 2L, 2L, 2L, 2L },
      { 111L, 42L, 54L, 118L, 20L, 31L, 31L, 72L, 17L, 94L, 14L, 14L, 11L, 14L, 14L, 14L, 94L, 10L, 10L, 10L, 10L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 11L, 11L, 11L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 18L, 18L, 18L, 18L, 18L, 113L, 62L, 62L, 45L, 45L, 113L, 113L, 113L, 113L, 113L, 113L, 113L, 113L, 113L, 113L, 113L, 113L, 113L, 113L, 113L, 113L, 113L, 63L, 63L, 53L, 63L, 67L, 67L, 67L, 67L, 67L, 67L, 67L, 67L, 67L, 67L, 67L, 67L, 67L, 67L, 67L, 51L, 51L, 51L, 51L, 51L, 12L, 51L, 12L, 51L, 5L, 3L, 3L, 2L, 2L, 5L },
      { 163L, 154L, 83L, 43L, 82L, 92L, 150L, 59L, 76L, 76L, 47L, 11L, 11L, 100L, 131L, 116L, 116L, 116L, 116L, 116L, 116L, 138L, 138L, 138L, 138L, 138L, 138L, 138L, 138L, 138L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 116L, 116L, 116L, 116L, 116L, 116L, 100L, 100L, 100L, 100L, 100L, 138L, 138L, 138L, 138L, 138L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 38L, 38L, 38L, 38L, 38L, 38L, 38L, 38L, 3L, 3L, 3L, 3L, 3L },
      { 246L, 189L, 242L, 102L, 250L, 250L, 102L, 250L, 280L, 118L, 196L, 118L, 191L, 215L, 121L, 121L, 49L, 49L, 49L, 49L, 49L, 49L, 49L, 49L, 49L, 49L, 49L, 49L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 171L, 161L, 161L, 161L, 161L, 161L, 161L, 161L, 161L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 14L, 10L, 10L, 10L, 10L, 10L, 10L, 103L, 10L, 10L, 10L, 10L, 5L },
      { 347L, 402L, 322L, 418L, 215L, 220L, 339L, 339L, 339L, 337L, 218L, 315L, 315L, 315L, 315L, 167L, 167L, 167L, 167L, 361L, 201L, 124L, 124L, 124L, 124L, 124L, 124L, 124L, 124L, 124L, 124L, 124L, 231L, 231L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 48L, 48L, 48L, 48L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 90L, 243L, 243L, 243L, 243L, 243L, 243L, 243L, 243L, 243L, 243L, 283L, 283L, 283L, 283L, 283L, 283L, 283L, 283L, 283L, 16L, 283L, 16L, 283L, 283L },
      { 505L, 220L, 601L, 644L, 612L, 160L, 206L, 206L, 206L, 422L, 134L, 518L, 134L, 134L, 518L, 652L, 382L, 206L, 158L, 441L, 179L, 441L, 56L, 559L, 559L, 56L, 56L, 56L, 56L, 56L, 56L, 56L, 56L, 56L, 56L, 56L, 56L, 56L, 56L, 101L, 101L, 56L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 193L, 193L, 193L, 193L, 193L, 193L, 193L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 122L, 122L, 122L, 122L, 122L, 122L, 122L, 122L, 122L, 122L, 122L, 122L, 122L, 122L, 122L, 122L, 122L, 101L, 101L, 101L, 101L },
      { 794L, 325L, 960L, 528L, 247L, 247L, 338L, 366L, 847L, 753L, 753L, 236L, 334L, 334L, 461L, 711L, 652L, 381L, 381L, 381L, 652L, 381L, 381L, 381L, 381L, 381L, 381L, 381L, 226L, 326L, 326L, 326L, 326L, 326L, 326L, 326L, 126L, 326L, 326L, 326L, 326L, 326L, 326L, 326L, 326L, 326L, 326L, 195L, 195L, 55L, 55L, 55L, 55L, 55L, 55L, 55L, 55L, 55L, 55L, 55L, 55L, 55L, 55L, 55L, 55L, 55L, 55L, 55L, 195L, 195L, 195L, 195L, 195L, 195L, 195L, 132L, 132L, 132L, 132L, 132L, 132L, 132L, 132L, 132L, 132L, 132L, 387L, 387L, 387L, 387L, 387L, 387L, 387L, 387L, 387L, 387L, 387L, 387L, 387L },
      { 1189L, 888L, 259L, 1082L, 725L, 811L, 636L, 965L, 497L, 497L, 1490L, 1490L, 392L, 1291L, 508L, 508L, 1291L, 1291L, 508L, 1291L, 508L, 508L, 867L, 867L, 867L, 867L, 934L, 867L, 867L, 867L, 867L, 867L, 867L, 867L, 1284L, 1284L, 1284L, 1284L, 1284L, 1284L, 1284L, 1284L, 1284L, 563L, 563L, 563L, 563L, 1010L, 1010L, 1010L, 208L, 838L, 563L, 563L, 563L, 759L, 759L, 564L, 759L, 759L, 801L, 801L, 801L, 801L, 759L, 759L, 759L, 759L, 759L, 563L, 563L, 563L, 563L, 563L, 563L, 563L, 563L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L, 226L },
      { 1763L, 1018L, 1500L, 432L, 1332L, 2203L, 126L, 2240L, 1719L, 1284L, 878L, 1983L, 266L, 266L, 266L, 266L, 747L, 747L, 127L, 127L, 2074L, 127L, 2074L, 1400L, 1383L, 1383L, 1383L, 1383L, 1383L, 1383L, 1383L, 1383L, 1383L, 1383L, 1400L, 1383L, 1383L, 1383L, 1383L, 1383L, 1383L, 1383L, 507L, 1073L, 1073L, 1073L, 1073L, 1990L, 1990L, 1990L, 1990L, 1990L, 507L, 507L, 507L, 507L, 507L, 507L, 507L, 507L, 507L, 1073L, 1073L, 1073L, 1073L, 1073L, 1073L, 1073L, 1073L, 1073L, 1073L, 1073L, 1073L, 1073L, 1073L, 1073L, 1073L, 1073L, 22L, 22L, 22L, 22L, 22L, 22L, 1073L, 452L, 452L, 452L, 452L, 452L, 452L, 318L, 301L, 301L, 301L, 301L, 86L, 86L, 15L },
      { 2872L, 3233L, 1534L, 2941L, 2910L, 393L, 1796L, 919L, 446L, 919L, 919L, 1117L, 103L, 103L, 103L, 103L, 103L, 103L, 103L, 2311L, 3117L, 1101L, 3117L, 3117L, 1101L, 1101L, 1101L, 1101L, 1101L, 2503L, 2503L, 2503L, 2503L, 2503L, 2503L, 2503L, 2503L, 429L, 429L, 429L, 429L, 429L, 429L, 429L, 1702L, 1702L, 1702L, 184L, 184L, 184L, 184L, 184L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 105L, 784L, 784L, 784L, 784L, 784L, 784L, 784L, 784L, 784L, 784L, 784L, 784L, 784L },
      { 4309L, 3758L, 4034L, 1963L, 730L, 642L, 1502L, 2246L, 3834L, 1511L, 1102L, 1102L, 1522L, 1522L, 3427L, 3427L, 3928L, 915L, 915L, 3818L, 3818L, 3818L, 3818L, 4782L, 4782L, 4782L, 3818L, 4782L, 3818L, 3818L, 1327L, 1327L, 1327L, 1327L, 1327L, 1327L, 1327L, 1387L, 1387L, 1387L, 1387L, 1387L, 1387L, 1387L, 1387L, 1387L, 2339L, 2339L, 2339L, 2339L, 2339L, 2339L, 2339L, 2339L, 2339L, 2339L, 2339L, 2339L, 2339L, 3148L, 3148L, 3148L, 3148L, 3148L, 3148L, 3148L, 3148L, 3148L, 3148L, 3148L, 3148L, 3148L, 3148L, 3148L, 3148L, 3148L, 3148L, 1776L, 1776L, 1776L, 3354L, 3354L, 3354L, 925L, 3354L, 3354L, 925L, 925L, 925L, 925L, 925L, 2133L, 2133L, 2133L, 2133L, 2133L, 2133L, 2133L, 2133L },
      { 6610L, 6977L, 1686L, 3819L, 2314L, 5647L, 3953L, 3614L, 5115L, 423L, 423L, 5408L, 7426L, 423L, 423L, 487L, 6227L, 2660L, 6227L, 1221L, 3811L, 197L, 4367L, 351L, 1281L, 1221L, 351L, 351L, 351L, 7245L, 1984L, 2999L, 2999L, 2999L, 2999L, 2999L, 2999L, 3995L, 2063L, 2063L, 2063L, 2063L, 1644L, 2063L, 2077L, 2512L, 2512L, 2512L, 2077L, 2077L, 2077L, 2077L, 754L, 754L, 754L, 754L, 754L, 754L, 754L, 754L, 754L, 754L, 754L, 754L, 754L, 754L, 754L, 754L, 754L, 754L, 754L, 1097L, 1097L, 754L, 754L, 754L, 754L, 248L, 754L, 1097L, 1097L, 1097L, 1097L, 222L, 222L, 222L, 222L, 754L, 1982L, 1982L, 1982L, 1982L, 1982L, 1982L, 1982L, 1982L, 1982L, 1982L, 1982L },
      { 9861L, 3647L, 4073L, 2535L, 3430L, 9865L, 2830L, 9328L, 4320L, 5913L, 10365L, 8272L, 3706L, 6186L, 7806L, 7806L, 7806L, 8610L, 2563L, 11558L, 11558L, 9421L, 1181L, 9421L, 1181L, 1181L, 1181L, 9421L, 1181L, 1181L, 10574L, 10574L, 3534L, 3534L, 3534L, 3534L, 3534L, 2898L, 2898L, 2898L, 3450L, 2141L, 2141L, 2141L, 2141L, 2141L, 2141L, 2141L, 7055L, 7055L, 7055L, 7055L, 7055L, 7055L, 7055L, 7055L, 7055L, 7055L, 7055L, 7055L, 7055L, 7055L, 7055L, 2831L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 8204L, 4688L, 4688L, 4688L, 2831L, 2831L, 2831L, 2831L, 2831L, 2831L, 2831L, 2831L },
      { 10327L, 7582L, 7124L, 8214L, 9600L, 10271L, 10193L, 10800L, 9086L, 2365L, 4409L, 13812L, 5661L, 9344L, 9344L, 10362L, 9344L, 9344L, 8585L, 11114L, 13080L, 13080L, 13080L, 6949L, 3436L, 3436L, 3436L, 13213L, 6130L, 6130L, 8159L, 8159L, 11595L, 8159L, 3436L, 7096L, 7096L, 7096L, 7096L, 7096L, 7096L, 7096L, 7096L, 7096L, 7096L, 7096L, 7096L, 7096L, 7096L, 7096L, 7096L, 7096L, 7096L, 4377L, 7096L, 4377L, 4377L, 4377L, 4377L, 4377L, 5410L, 5410L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 4377L, 440L, 440L, 1199L, 1199L, 1199L },
      { 19540L, 19926L, 11582L, 11113L, 24585L, 8726L, 17218L, 419L, 4918L, 4918L, 4918L, 15701L, 17710L, 4037L, 4037L, 15808L, 11401L, 19398L, 25950L, 25950L, 4454L, 24987L, 11719L, 8697L, 1452L, 1452L, 1452L, 1452L, 1452L, 8697L, 8697L, 6436L, 21475L, 6436L, 22913L, 6434L, 18497L, 11089L, 11089L, 11089L, 11089L, 3036L, 3036L, 14208L, 14208L, 14208L, 14208L, 12906L, 12906L, 12906L, 12906L, 12906L, 12906L, 12906L, 12906L, 7614L, 7614L, 7614L, 7614L, 5021L, 5021L, 5021L, 5021L, 5021L, 5021L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 10145L, 4544L, 4544L, 4544L, 4544L, 4544L, 4544L, 8394L, 8394L, 8394L, 8394L },
      { 34566L, 9579L, 12654L, 26856L, 37873L, 38806L, 29501L, 17271L, 3663L, 10763L, 18955L, 1298L, 26560L, 17132L, 17132L, 4753L, 4753L, 8713L, 18624L, 13082L, 6791L, 1122L, 19363L, 34695L, 18770L, 18770L, 18770L, 18770L, 15628L, 18770L, 18770L, 18770L, 18770L, 33766L, 20837L, 20837L, 20837L, 20837L, 20837L, 20837L, 6545L, 6545L, 6545L, 6545L, 6545L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 30483L, 30483L, 30483L, 30483L, 30483L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 12138L, 9305L, 11107L, 11107L, 11107L, 11107L, 11107L, 11107L, 11107L, 11107L, 11107L, 11107L, 11107L, 11107L, 11107L, 9305L, 9305L },
      { 31929L, 49367L, 10982L, 3527L, 27066L, 13226L, 56010L, 18911L, 40574L, 20767L, 20767L, 9686L, 47603L, 47603L, 11736L, 11736L, 41601L, 12888L, 32948L, 30801L, 44243L, 53351L, 53351L, 16016L, 35086L, 35086L, 32581L, 2464L, 2464L, 49554L, 2464L, 2464L, 49554L, 49554L, 2464L, 81L, 27260L, 10681L, 2185L, 2185L, 2185L, 2185L, 2185L, 2185L, 2185L, 18086L, 18086L, 18086L, 18086L, 18086L, 17631L, 17631L, 18086L, 18086L, 18086L, 37335L, 37774L, 37774L, 37774L, 26401L, 26401L, 26401L, 26401L, 26401L, 26401L, 26401L, 26401L, 26401L, 26401L, 26401L, 26401L, 26401L, 12982L, 40398L, 40398L, 40398L, 40398L, 40398L, 40398L, 3518L, 3518L, 3518L, 37799L, 37799L, 37799L, 37799L, 37799L, 37799L, 37799L, 37799L, 37799L, 4721L, 4721L, 4721L, 4721L, 7067L, 7067L, 7067L, 7067L },
      { 40701L, 69087L, 77576L, 64590L, 39397L, 33179L, 10858L, 38935L, 43129L, 35468L, 35468L, 5279L, 61518L, 61518L, 27945L, 70975L, 70975L, 86478L, 86478L, 20514L, 20514L, 73178L, 73178L, 43098L, 43098L, 4701L, 59979L, 59979L, 58556L, 69916L, 15170L, 15170L, 4832L, 4832L, 43064L, 71685L, 4832L, 15170L, 15170L, 15170L, 27679L, 27679L, 27679L, 60826L, 60826L, 6187L, 6187L, 4264L, 4264L, 4264L, 4264L, 4264L, 45567L, 32269L, 32269L, 32269L, 32269L, 62060L, 62060L, 62060L, 62060L, 62060L, 62060L, 62060L, 62060L, 62060L, 1803L, 1803L, 1803L, 1803L, 1803L, 1803L, 1803L, 1803L, 1803L, 1803L, 1803L, 1803L, 1803L, 51108L, 51108L, 51108L, 51108L, 51108L, 51108L, 51108L, 51108L, 51108L, 51108L, 51108L, 51108L, 55315L, 55315L, 54140L, 54140L, 54140L, 54140L, 54140L, 13134L },
      { 103650L, 125480L, 59978L, 46875L, 77172L, 83021L, 126904L, 14541L, 56299L, 43636L, 11655L, 52680L, 88549L, 29804L, 101894L, 113675L, 48040L, 113675L, 34987L, 48308L, 97926L, 5475L, 49449L, 6850L, 62545L, 62545L, 9440L, 33242L, 9440L, 33242L, 9440L, 33242L, 9440L, 62850L, 9440L, 9440L, 9440L, 90308L, 90308L, 90308L, 47904L, 47904L, 47904L, 47904L, 47904L, 47904L, 47904L, 47904L, 47904L, 41143L, 41143L, 41143L, 41143L, 41143L, 41143L, 41143L, 36114L, 36114L, 36114L, 36114L, 36114L, 24997L, 65162L, 65162L, 65162L, 65162L, 65162L, 65162L, 65162L, 65162L, 65162L, 65162L, 65162L, 65162L, 65162L, 65162L, 47650L, 47650L, 47650L, 47650L, 47650L, 47650L, 47650L, 40586L, 40586L, 40586L, 40586L, 40586L, 40586L, 40586L, 38725L, 38725L, 38725L, 38725L, 88329L, 88329L, 88329L, 88329L, 88329L },
      { 165843L, 90647L, 59925L, 189541L, 67647L, 74795L, 68365L, 167485L, 143918L, 74912L, 167289L, 75517L, 8148L, 172106L, 126159L, 35867L, 35867L, 35867L, 121694L, 52171L, 95354L, 113969L, 113969L, 76304L, 123709L, 123709L, 144615L, 123709L, 64958L, 64958L, 32377L, 193002L, 193002L, 25023L, 40017L, 141605L, 189165L, 189165L, 141605L, 189165L, 189165L, 141605L, 141605L, 141605L, 189165L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127047L, 127785L, 127785L, 127785L, 127785L, 127785L, 127785L, 127785L, 127785L, 127785L, 127785L, 80822L, 80822L, 80822L, 80822L, 80822L, 80822L, 131661L, 131661L, 131661L, 131661L, 131661L, 131661L, 131661L, 131661L, 131661L, 131661L, 131661L, 131661L, 131661L, 131661L, 131661L, 131661L, 7114L, 131661L },
      { 130365L, 236711L, 110235L, 125699L, 56483L, 93735L, 234469L, 60549L, 1291L, 93937L, 245291L, 196061L, 258647L, 162489L, 176631L, 204895L, 73353L, 172319L, 28881L, 136787L, 122081L, 122081L, 275993L, 64673L, 211587L, 211587L, 211587L, 282859L, 282859L, 211587L, 242821L, 256865L, 256865L, 256865L, 122203L, 291915L, 122203L, 291915L, 291915L, 122203L, 25639L, 25639L, 291803L, 245397L, 284047L, 245397L, 245397L, 245397L, 245397L, 245397L, 245397L, 245397L, 94241L, 66575L, 66575L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 217673L, 210249L, 210249L, 210249L, 210249L, 210249L, 210249L, 210249L, 210249L, 210249L, 210249L, 94453L, 94453L, 94453L, 94453L, 94453L, 94453L, 94453L, 94453L, 94453L, 94453L, 94453L, 94453L, 94453L, 94453L, 94453L },
      { 333459L, 375354L, 102417L, 383544L, 292630L, 41147L, 374614L, 48032L, 435453L, 281493L, 358168L, 114121L, 346892L, 238990L, 317313L, 164158L, 35497L, 70530L, 70530L, 434839L, 24754L, 24754L, 24754L, 393656L, 118711L, 118711L, 148227L, 271087L, 355831L, 91034L, 417029L, 417029L, 91034L, 91034L, 417029L, 91034L, 299843L, 299843L, 413548L, 413548L, 308300L, 413548L, 413548L, 413548L, 308300L, 308300L, 308300L, 413548L, 308300L, 308300L, 308300L, 308300L, 308300L, 15311L, 15311L, 15311L, 15311L, 176255L, 176255L, 23613L, 23613L, 23613L, 23613L, 23613L, 23613L, 172210L, 204328L, 204328L, 204328L, 204328L, 121626L, 121626L, 121626L, 121626L, 121626L, 200187L, 200187L, 200187L, 200187L, 200187L, 121551L, 121551L, 248492L, 248492L, 248492L, 248492L, 248492L, 248492L, 248492L, 248492L, 248492L, 248492L, 248492L, 248492L, 13942L, 13942L, 13942L, 13942L, 13942L },
      { 500884L, 566009L, 399251L, 652979L, 355008L, 430235L, 328722L, 670680L, 405585L, 405585L, 424646L, 670180L, 670180L, 641587L, 215580L, 59048L, 633320L, 81010L, 20789L, 389250L, 389250L, 638764L, 638764L, 389250L, 389250L, 398094L, 80846L, 147776L, 147776L, 296177L, 398094L, 398094L, 147776L, 147776L, 396313L, 578233L, 578233L, 578233L, 19482L, 620706L, 187095L, 620706L, 187095L, 126467L, 241663L, 241663L, 241663L, 241663L, 241663L, 241663L, 241663L, 241663L, 241663L, 241663L, 241663L, 241663L, 321632L, 23210L, 23210L, 394484L, 394484L, 394484L, 78101L, 78101L, 78101L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 542095L, 277743L, 277743L, 277743L, 457259L, 457259L, 457259L, 457259L, 457259L, 457259L, 457259L, 457259L, 457259L, 457259L, 457259L, 457259L },
      { 858339L, 918142L, 501970L, 234813L, 460565L, 31996L, 753018L, 256150L, 199809L, 993599L, 245149L, 794183L, 121349L, 150619L, 376952L, 809123L, 809123L, 804319L, 67352L, 969594L, 434796L, 969594L, 804319L, 391368L, 761041L, 754049L, 466264L, 754049L, 754049L, 466264L, 754049L, 754049L, 282852L, 429907L, 390017L, 276645L, 994856L, 250142L, 144595L, 907454L, 689648L, 687580L, 687580L, 687580L, 687580L, 978368L, 687580L, 552742L, 105195L, 942843L, 768249L, 307142L, 307142L, 307142L, 307142L, 880619L, 880619L, 880619L, 880619L, 880619L, 880619L, 880619L, 117185L, 117185L, 117185L, 117185L, 117185L, 117185L, 117185L, 117185L, 117185L, 117185L, 117185L, 60731L, 60731L, 60731L, 60731L, 60731L, 60731L, 60731L, 60731L, 60731L, 60731L, 60731L, 178309L, 178309L, 178309L, 178309L, 74373L, 74373L, 74373L, 74373L, 74373L, 74373L, 74373L, 74373L, 214965L, 214965L, 214965L }
    };

    // working objects.
    int * const pr = imem.get_mem();

    double * const __restrict__ finval     = dmem.get_mem(),
           * const __restrict__ M          = finval + nf,
           * const __restrict__ finest_var = M + nf,
           * const __restrict__ kahan_comp = finest_var + nf,
           * const __restrict__ x          = kahan_comp + nf,
           * const __restrict__ r          = x + ndim * n_qmc_seqs(),
           * const __restrict__ vk         = r + ndim,
           * const __restrict__ values     = vk + ndim,
           * const __restrict__ fs         = values + nf;

    // initialize
    std::fill(finest    , finest     + nf, 0);
    std::fill(finest_var, finest_var + nf, 0);
    int sampls = n_sequences;
    int np(0L);
    {
      int i = 0L;
      for(; i < plim; ++i){
        np = i;
        if(minvls < static_cast<size_t>(2L * sampls * p[i]))
          break;
      }
      if(i >= plim)
        sampls = std::max<int>(
          n_sequences, static_cast<int>(minvls / (2L * p[np])));
    }

    constexpr size_t const maxit(1000L);
    size_t intvls(0L);
    int inform(1L);
    double abserr(std::numeric_limits<double>::infinity());

    for(size_t nit = 0; nit < maxit; nit++){
      *vk = 1. / static_cast<double>(p[np]);
      int k(1L);
      for(int i = 1; i < ndim; ++i){
        if(i < klim){
          k = fmod(
            C[np][std::min(ndim, klim) - 2L] *
              static_cast<double>(k), static_cast<double>(p[np]));
          vk[i] = k * vk[0];

        } else {
          vk[i] = static_cast<int>(
            p[np] * std::pow(2., static_cast<double>(i + 1 - klim) /
            static_cast<double>(ndim - klim + 1)));
          vk[i] = fmod(vk[i] / p[np], 1.);

        }
      }

      // reset variance estimator and estimator of the output
      std::fill(finval, finval + nf, 0.);
      std::fill(M     , M      + nf, 0.);

      auto mvkrsv =
        [&](double * __restrict__ const values, int const prime,
            double const * __restrict__ const vk,
            double * __restrict__ const x,
            double * __restrict__ const r,
            int * __restrict__  const pr,
            double * __restrict__  const fs){
          std::fill(values, values + nf, 0.);
          std::fill(kahan_comp, kahan_comp + nf, 0.);

          // random shift
          {
            double * rj = r;
            int * prj = pr;
            for(int j = 0; j < ndim; ++j, ++rj, ++prj){
              *rj = sampler();
              if(j < klim - 1L){
                int const jp = (j + 1) * *rj;
                if(jp < j)
                  *prj = pr[jp];
                pr[jp] = j;
              } else
                *prj = j;
            }
          }

          // apply lattice rule
          for(int k = 0; k < prime; ){
            // compute the points
            int i = 0;
            double * x_odd  = x,
                   * x_even = x_odd + ndim;
            for(; i < n_qmc_seqs() / 2 and k < prime;
                ++i, ++k, x_odd += 2 * ndim, x_even += 2 * ndim){
              for(int j = 0; j < ndim; ++j){
                r[j] += vk[pr[j]];
                if(r[j] > 1.)
                  r[j] -= 1.;
                x_odd [j] = std::abs(2 * r[j] - 1);
                x_even[j] = 1 - x_odd[j];
              }
            }

            // evaluate the integrand
            f(&ndim, x, &nf, fs, 2 * i);

            // update the sum (values)
            double *fsk = fs;
            for(int k = 0; k < 2 * i; ++k, fsk += nf)
              for(int j = 0; j < nf; ++j)
                kahan(values[j], kahan_comp[j], fsk[j]);
          }
          for(int j = 0; j < nf; ++j)
            values[j] /= static_cast<double>(2 * prime);
        };

      for(int i = 0; i < sampls; ++i){
        mvkrsv(values, p[np], vk, x, r, pr, fs);
        // stable version of Welford's online algorithm
        for(int k = 0; k < nf; ++k){
          double const term_diff = values[k] - finval[k];
          finval[k] +=  term_diff / (i + 1.);
          M     [k] += term_diff * (values[k] - finval[k]);
        }
      }

      intvls += 2 * sampls * p[np];
      bool passes_conv_check = true;
      for(int k = 0; k < nf; ++k){
        // update the mean estimator and variance estimator
        double const sig_new =
          M[k] / (sampls - 1.) / static_cast<double>(sampls);
        if(finest_var[k] <= 0){
          // no prior term or deterministic
          finest[k]     = finval[k];
          finest_var[k] = sig_new;

        } else {
          // there is a previous estimator. Update mean
          double const sig_old = finest_var[k],
                       w_denom = 1 / sig_old + 1 / sig_new;
          finest[k] = (finest[k] / sig_old + finval[k] / sig_new) / w_denom;

          // update the variance
          finest_var[k] = sig_old * sig_new / (sig_old + sig_new);
        }

        // passes criteria
        sdest[k] = std::sqrt(finest_var[k]);
        abserr = 7 / 2 * sdest[k];
        passes_conv_check &=
          abserr <= std::max(abseps, std::abs(finest[k]) * releps);
      }

      if(!passes_conv_check){
        if(np < plim - 1)
          ++np;
        else {
          sampls = std::min(3 * sampls / 2,
                            static_cast<int>((maxvls - intvls) / (2 * p[np])));
          sampls = std::max<int>(n_sequences, sampls);
        }
        if(intvls + 2 * sampls * p[np] > maxvls)
          break;
      } else {
        inform = 0L;
        break;
      }
    }

    return { intvls, abserr, inform };
  }
};

template<class Func>
cache_mem<double> rand_Korobov<Func>::dmem;
template<class Func>
cache_mem<int   > rand_Korobov<Func>::imem;

template<class Func>
class sobol_wrapper {
  static cache_mem<double> dmem;
  static unsigned max_n_sequences;

public:
  static void alloc_mem
  (int const max_ndim, int const max_nf, int const max_threads,
   unsigned const max_n_sequences_in) {
    max_n_sequences = std::max(max_n_sequences, max_n_sequences_in);
    dmem.set_n_mem(
      (n_qmc_seqs() + max_n_sequences) * max_nf +
        n_qmc_seqs() * max_ndim, max_threads);
  }

  static rand_Korobov_output comp
  (Func &f, int const ndim, size_t const minvls, size_t const maxvls,
   int const nf, double const abseps, double const releps,
   double * const __restrict__ finest,
   double * const __restrict__ sdest, parallelrng::unif_drawer &sampler,
   sobol::scrambling_type const method, unsigned const n_sequences){
    if(n_sequences < 2)
      throw std::invalid_argument("n_sequences is less than two");
    if(n_sequences > max_n_sequences)
      throw std::invalid_argument("n_sequences is larger then the set max_n_sequences");
    if(method == sobol::scrambling_type::none)
      throw std::invalid_argument("sobol::scrambling_type::none passed but it makes no sense");

    // the variables we will return
    size_t intvls(0L);
    int inform(1L);
    double abserr(std::numeric_limits<double>::infinity());

    // get the sequences we need
    std::vector<sobol> seqs;
    seqs.reserve(n_sequences);
    for(unsigned i = 0; i < n_sequences; ++i){
      int const i_seed =
        static_cast<int>(std::numeric_limits<int>::max() * sampler());
      seqs.emplace_back(ndim, method, i_seed < 0 ? 1 : i_seed);
    }

    // initialize the objects we need
    double * const __restrict__ integrand_vals = dmem.get_mem(),
           * const __restrict__ seq_means      = integrand_vals + nf * n_qmc_seqs(),
           * const __restrict__ draws          = seq_means + nf * max_n_sequences;
    std::fill(seq_means, seq_means + nf * n_sequences, 0);

    // main loop where we compute the result
    unsigned n_draw_next = minvls / n_sequences + 1L,
             n_drawn_per_seq = 0;
    for(; intvls < maxvls; ){
      // update the estimator for each of the Sobol sequences
      for(unsigned i = 0; i < n_sequences; ++i){
        double denom(n_drawn_per_seq);

        for(size_t k = 0; k < n_draw_next;){
          // get the next points
          double * __restrict__ d = draws;
          int const n_draw_j =  std::min<int>(n_qmc_seqs(), n_draw_next - k);
          for(int j = 0; j <n_draw_j; ++j, ++k, d += ndim){
            seqs[i].next();
            std::copy(seqs[i].get_qausi(), seqs[i].get_qausi() + ndim, d);
          }

          // evaluate the integrands
          f(&ndim, draws, &nf, integrand_vals, n_draw_j);

          // update the estimates
          double * int_val = integrand_vals;
          for(int l = 0; l < n_draw_j; ++l, int_val += nf){
            denom += 1;
            double *meas = seq_means + i * nf;
            for(int j = 0; j < nf; ++j)
              meas[j] += (int_val[j] - meas[j]) / denom;
          }
        }

        intvls += n_draw_next;
      }

      // update the denominator counter
      n_drawn_per_seq += n_draw_next;

      // compute the mean estimator and the metric to compute the standard error
      double * const __restrict__ M = integrand_vals;
      std::fill(finest, finest + nf, 0);
      std::fill(M     , M      + nf, 0);

      for(unsigned i = 0; i < n_sequences; ++i){
        // stable version of Welford's online algorithm
        double *values = seq_means + i * nf;
        for(int j = 0; j < nf; ++j){
          double const term_diff = values[j] - finest[j];
          finest[j] +=  term_diff / (i + 1.);
          M     [j] += term_diff * (values[j] - finest[j]);
        }
      }

      // check if we can exit early
      bool passes_conv_check = true;
      for(int j = 0; j < nf; ++j){
        double const sigma =
          M[j] / (n_sequences - 1.) / static_cast<double>(n_sequences);
        sdest[j] = std::sqrt(sigma);
        abserr = 7 / 2 * sdest[j];
        passes_conv_check &=
          abserr <= std::max(abseps, std::abs(finest[j]) * releps);
      }

      // exit or compute the number of samples to compute next time
      if(!passes_conv_check){
        n_draw_next = (n_draw_next * 3L) / 2L + 1L;
        if(n_draw_next * n_sequences + intvls > maxvls)
          n_draw_next = (maxvls - intvls) / n_sequences + 1L;
      } else {
        inform = 0L;
        break;
      }
    }

    return { intvls, abserr, inform };
  }
};

template<class Func>
cache_mem<double> sobol_wrapper<Func>::dmem;
template<class Func>
unsigned sobol_wrapper<Func>::max_n_sequences = 2;

} // namespace pedmod

#endif
