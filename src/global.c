// global.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license
// agreement,
// provided that this statement is retained.
//
// Global variables and functions used by BLAST

#include "blast.h"
#include <errno.h>

// Timing variables
int4 blast_prepTime = 0, blast_searchTime = 0;
int4 blast_gappedScoreTime = 0, blast_gappedExtendTime = 0,
     blast_finalizeTime = 0;
int4 blast_semiGappedScoreTime = 0, blast_copyTime = 0, blast_unpackTime = 0;
uint4 blast_numQuery = 0;
uint4 blast_numBlocks = 0;
uint2 blast_tgSize = 1;

int rank = 0, num_procs = 1;

uint4 longestQueryLength;
uint8 total_numberOfLetters = 0;

// BLAST statistics
uint4 blast_numHits = 0;
uint4 blast_numUngappedExtensions = 0, blast_numTriggerExtensions = 0,
      blast_numTriggerSequences = 0;
uint4 blast_numGapped = 0;
uint4 blast_numSemiGapped = 0;
uint4 blast_numExtensionsPruned = 0;
uint4 blast_numAttemptedJoin = 0, blast_numSuccessfullyJoined = 0;
uint4 blast_numGoodAlignments = 0;
uint4 blast_numGoodExtensions = 0;
uint4 blast_totalUnpacked = 0;
uint4 blast_totalCopied = 0;
uint4 blast_numExpandedSequences = 0;

// BLAST statistics
uint4 blast_numHits_multi[BATCH_SIZE];
uint4 blast_numUngappedExtensions_multi[BATCH_SIZE],
    blast_numTriggerExtensions_multi[BATCH_SIZE],
    blast_numTriggerSequences_multi[BATCH_SIZE];
uint4 blast_numGapped_multi[BATCH_SIZE];
uint4 blast_numGappedExtension_multi[BATCH_SIZE];
uint4 blast_numSemiGapped_multi[BATCH_SIZE];
uint4 blast_numExtensionsPruned_multi[BATCH_SIZE];
uint4 blast_numAttemptedJoin_multi[BATCH_SIZE],
    blast_numSuccessfullyJoined_multi[BATCH_SIZE];
uint4 blast_numGoodAlignments_multi[BATCH_SIZE];
uint4 blast_numGoodExtensions_multi[BATCH_SIZE];
uint4 blast_totalUnpacked_multi[BATCH_SIZE];
uint4 blast_totalCopied_multi[BATCH_SIZE];
uint4 blast_numExpandedSequences_multi[BATCH_SIZE];

// long long blast_hitDetectCycle = 0, blast_pesudoCycle = 0,
// blast_ungappedExtCycle = 0, blast_sortCycle = 0, blast_findGoodCycle = 0,
// blast_findFinalCycle = 0, blast_getFinalCycle = 0, blast_getTrackCycle = 0;

// BLAST global variables
int4 blast_ungappedNominalTrigger;
int4 blast_gappedNominalCutoff = 0, blast_nominalR1cutoff = 0,
     blast_nominalR2cutoff = 0;
int4 blast_dynamicGappedNominalCutoff = 0, blast_dynamicNominalR1cutoff = 0;
int4 blast_dloc;
char *blast_queryDescription;

// BLAST global variables
int4 blast_ungappedNominalTrigger_multi[BATCH_SIZE];
int4 blast_gappedNominalCutoff_multi[BATCH_SIZE] = { 0 },
     blast_nominalR1cutoff_multi[BATCH_SIZE] = { 0 },
     blast_nominalR2cutoff_multi[BATCH_SIZE] = { 0 };
int4 blast_dynamicGappedNominalCutoff_multi[BATCH_SIZE] = { 0 },
     blast_dynamicNominalR1cutoff_multi[BATCH_SIZE] = { 0 };
int4 blast_dloc_multi[BATCH_SIZE];
char *blast_queryDescription_multi[BATCH_SIZE] = {0};

uint4 dbIdx_block_size = 0;

char *global_string = NULL;

// Initialize global variables
void global_initialize_multi() {
  blast_prepTime = 0;
  blast_searchTime = 0;
  blast_gappedScoreTime = 0;
  blast_gappedExtendTime = 0;
  blast_finalizeTime = 0;
  blast_semiGappedScoreTime = 0;
  blast_copyTime = 0;
  blast_unpackTime = 0;

  // BLAST statistics
  memset(blast_numHits_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_numUngappedExtensions_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_numTriggerExtensions_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_numTriggerSequences_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_numGapped_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_numGappedExtension_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_numSemiGapped_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_numExtensionsPruned_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_numAttemptedJoin_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_numSuccessfullyJoined_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_numGoodAlignments_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_numGoodExtensions_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_totalUnpacked_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_totalCopied_multi, 0, BATCH_SIZE * sizeof(uint4));
  memset(blast_numExpandedSequences_multi, 0, BATCH_SIZE * sizeof(uint4));

  // BLAST global variables
  memset(blast_gappedNominalCutoff_multi, 0, BATCH_SIZE * sizeof(int4));
  memset(blast_nominalR1cutoff_multi, 0, BATCH_SIZE * sizeof(int4));
  memset(blast_nominalR2cutoff_multi, 0, BATCH_SIZE * sizeof(int4));
  memset(blast_dynamicGappedNominalCutoff_multi, 0, BATCH_SIZE * sizeof(int4));
  memset(blast_dynamicNominalR1cutoff_multi, 0, BATCH_SIZE * sizeof(int4));
}

// Initialize global variables
void global_initialize() {
  blast_prepTime = 0;
  blast_searchTime = 0;
  blast_gappedScoreTime = 0;
  blast_gappedExtendTime = 0;
  blast_finalizeTime = 0;
  blast_semiGappedScoreTime = 0;
  blast_copyTime = 0;
  blast_unpackTime = 0;

  // BLAST statistics
  blast_numHits = 0;
  blast_numUngappedExtensions = 0;
  blast_numTriggerExtensions = 0;
  blast_numTriggerSequences = 0;
  blast_numGapped = 0;
  blast_numSemiGapped = 0;
  blast_numExtensionsPruned = 0;
  blast_numAttemptedJoin = 0;
  blast_numSuccessfullyJoined = 0;
  blast_numGoodAlignments = 0;
  blast_numGoodExtensions = 0;
  blast_totalUnpacked = 0;
  blast_totalCopied = 0;
  blast_numExpandedSequences = 0;

  // BLAST global variables
  blast_gappedNominalCutoff = 0;
  blast_nominalR1cutoff = 0;
  blast_nominalR2cutoff = 0;
  blast_dynamicGappedNominalCutoff = 0;
  blast_dynamicNominalR1cutoff = 0;
}

// Convert a 32-bit integer into a string with commas
char *global_int4toString(uint4 number) {
  char string1[50];
  int4 length, count1, count2;

  // Convert integer to string
  sprintf(string1, "%u", number);
  length = strlen(string1);

  // Declare second string large enough to hold number with commas
  global_string =
      (char *)global_realloc(global_string, sizeof(char) * length * 4 / 3 + 1);
  global_string[0] = '\0';

  count1 = count2 = 0;

  // Convert number to version with commas
  while (count1 < length) {
    global_string[count2] = string1[count1];
    count1++;
    count2++;

    if (number >= 10000 && ((length - count1) % 3 == 0 && count1 < length)) {
      global_string[count2] = ',';
      count2++;
    }
  }

  // Null terminate string2
  global_string[count2] = '\0';

  return global_string;
}

// Convert a 64-bit integer into a string with commas
char *global_int8toString(uint8 number) {
  char string1[50];
  int4 length, count1, count2;

  // Convert integer to string
  sprintf(string1, "%lu", number);
  length = strlen(string1);

  // Declare second string large enough to hold number with commas
  global_string =
      (char *)global_realloc(global_string, sizeof(char) * length * 4 / 3 + 1);
  global_string[0] = '\0';

  count1 = count2 = 0;

  // Convert number to version with commas
  while (count1 < length) {
    global_string[count2] = string1[count1];
    count1++;
    count2++;

    if (number >= 10000 && ((length - count1) % 3 == 0 && count1 < length)) {
      global_string[count2] = ',';
      count2++;
    }
  }

  // Null terminate string2
  global_string[count2] = '\0';

  return global_string;
}

// Free the global convert string
void global_free() { free(global_string); }

uint8 global_totalMalloc = 0;

// Malloc new memory, check that malloc was successful
void *global_malloc(size_t size) {

  void *newMemory;

  newMemory = malloc(size);

  if (newMemory == NULL && size != 0) {

    // Report error allocating memory
    fprintf(stderr, "Error allocating %lu Megabytes: ", size >> 20);
    fprintf(stderr, "%s\n", strerror(errno));
    fprintf(stderr, "\n");
    fflush(stderr);
    exit(-1);
  }

  global_totalMalloc += size;

  return newMemory;
}

// Realloc memory, check that realloc was successful
void *global_realloc(void *ptr, size_t size) {
    ptr = realloc(ptr, size);
    if (ptr == NULL && size != 0) {
        // Report error allocating memory
        fprintf(stderr, "Error allocating %lu bytes: ", size);
        fprintf(stderr, "%s\n", strerror(errno));
        fflush(stderr);
        exit(-1);
    }

    // global_totalMalloc += size;
    __sync_fetch_and_add(&global_totalMalloc, size);

    return ptr;
}
