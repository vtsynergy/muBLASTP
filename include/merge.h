#ifndef _MERGE_H_
#define _MERGE_H__

#include<stdint.h>

extern "C"{
void merge(const uint32_t *a, int aCount, const uint32_t *b, int bCount,
        int *dest);

void merge_2(const uint32_t *a, int aCount, const uint32_t *b, int bCount,
        int *dest);
void std_merge(const unsigned int *a, int aCount, const unsigned int *b, int bCount, unsigned int *dest);
void std_sort(int *a, int aCount);
}

#endif
