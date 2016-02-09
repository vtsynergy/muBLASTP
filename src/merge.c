#include "merge.h"
#include <algorithm>

void std_merge(const unsigned int *a, int aCount, const unsigned int *b,
               int bCount, unsigned int *dest) {
  std::merge(a, a + aCount, b, b + bCount, dest);
}

bool compInt(int i, int j) { return (i < j); }

bool compHit(int i, int j) {
  int t1 = ((i & 0xffff) << 16) + (i >> 16);
  int t2 = ((j & 0xffff) << 16) + (j >> 16);
  return (t1 < t2);
}

void std_sort(int *a, int aCount) { std::sort(a, a + aCount, compInt); }

void hit_sort(int *a, int aCount) { std::sort(a, a + aCount, compHit); }
