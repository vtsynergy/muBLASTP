#include <omp.h>
#include "blast.h"
#include <sched.h>
#include <stdlib.h>
#include <stdio.h>

int get_affinity(int *cpu_map) {
  char *kmp;
  kmp = getenv("KMP_AFFINITY");
#if defined(__ICC) || defined(__INTEL_COMPILER)
  if (kmp != NULL) {
    printf("Detect KMP_AFFINITY env: %s\n", kmp);
#pragma omp parallel
    { cpu_map[omp_get_thread_num()] = sched_getcpu(); }
  } else
#endif
  {
    printf("KMP_AFFINITY unavailable, use default affinity\n");
    int tt;
    for (tt = 0; tt < parameters_num_threads; tt++)
      cpu_map[tt] = tt;
  }
  return 0;
}

int comparePos_subjectOffset(const void *a, const void *b) {
    int pos_a = *((uint32_t *)a);
    int pos_b = *((uint32_t *)b);
    return pos_a - pos_b;
}

int compareWord(const void *a, const void *b) {
    struct initialWord_protein_db *word_a = (struct initialWord_protein_db *)a;
    struct initialWord_protein_db *word_b = (struct initialWord_protein_db *)b;
    return word_b->numSubPositions - word_a->numSubPositions;
}

int comparePos_seqId(const void *a, const void *b) {
    int pos_a = *((uint32_t *)a) & 0xffff;
    int pos_b = *((uint32_t *)b) & 0xffff;
    return pos_b - pos_a;
}

