#include <x86intrin.h>
#include "score_functions.h"
// #include "libpopcnt.h"

void tanimoto_score(const uint64_t* fp1, const uint64_t* fp2, size_t size, double *t_score) {
  uint64_t fp_un = 0;
  uint64_t fp_int = 0;
  for (size_t i = 0; i < size; i++) {
    fp_un += _mm_popcnt_u64(fp1[i] | fp2[i]);
    fp_int += _mm_popcnt_u64(fp1[i] & fp2[i]);
  }
  *t_score = (double)fp_int / (double)fp_un;
}
