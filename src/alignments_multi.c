/*
* muBLASTP - A database indexed protein sequence search tool 
* Version 0.1 (beta)
*
* (c) 2010-2015 Virginia Tech
* This version of dbBLASTP is licensed for non-commercial use only,
*  as specified in LICENSE. For all other use contact vtiplicensing@vtip.org
* 
* Authors: Jing Zhang 
*
*/


#define _GNU_SOURCE
#include "blast.h"
#include<sys/time.h>

struct memSingleBlock *alignments_finalAlignments_multi[BATCH_SIZE];

void alignments_addGoodAlignment_multi(int4 highestNominalScore,
        struct alignment *alignment,
        int queryNum);

void alignments_sortFinalAlignments_multi(int queryNum);

void alignments_checkForJoin_multi(struct alignment *alignment,
        struct ungappedExtension *extension1,
        struct PSSMatrix PSSMatrix, int queryNum);

void alignments_regularGappedAlignment_multi(
        struct PSSMatrix PSSMatrix, struct ungappedExtension *ungappedExtension,
        struct alignment *alignment, int queryNum);

int alignments_expandCluster_multi(struct alignment *alignment,
        struct PSSMatrix PSSMatrix, int queryNum);

// Initialize array storing pointers to alignments
void alignments_initialize_multi2() {
    int ii;
    for (ii = 0; ii < blast_numQuery; ii++) {
        alignments_finalAlignments_multi[ii] = memSingleBlock_initialize(
                sizeof(struct finalAlignment), constants_initialAllocFinalAlignments);
    }
}

// Initialize array storing pointers to alignments
void alignments_initialize_multi() {
    int ii;
    // Initialize alignments, good alignments and final alignments blocks
    for (ii = 0; ii < blast_numQuery; ii++) {
        alignments_finalAlignments_multi[ii] = memSingleBlock_initialize(
                sizeof(struct finalAlignment), constants_initialAllocFinalAlignments);
    }
}


int4 alignments_compareFinalAlignments2(const void *alignment1,
                                       const void *alignment2) {
  const struct finalAlignment *a1, *a2;

  a1 = (struct finalAlignment *)alignment1;
  a2 = (struct finalAlignment *)alignment2;

  if (a1->alignment->gappedExtensions->eValue > a2->alignment->gappedExtensions->eValue) {
    return 1;
  } else if (a1->alignment->gappedExtensions->eValue < a2->alignment->gappedExtensions->eValue) {
    return -1;
  } else {
    // Resolve conflicts using subject length
      if (a1->highestNominalScore > a2->highestNominalScore) {
          return -1;
      } else if (a1->highestNominalScore < a2->highestNominalScore) {
          return 1;
      }
      else
          return 0;
  }
}

void alignments_sortFinalAlignments_multi2(int queryNum) {
    qsort(alignments_finalAlignments_multi[queryNum]->block,
            alignments_finalAlignments_multi[queryNum]->numEntries,
            sizeof(struct finalAlignment), alignments_compareFinalAlignments2);
}

// Add the current alignment (which contains at least one gapped extension
// scoring above cutoff) to to-be-sorted list of final alignments
struct finalAlignment * alignments_addFinalAlignment_multi(int4 highestNominalScore,
        struct alignment *alignment,
        int queryNum) {
    struct finalAlignment *finalAlignment;

    // Get a new final alignment entry
    finalAlignment =
        memSingleBlock_newEntry(alignments_finalAlignments_multi[queryNum]);

    // Insert alignment information into the new final alignment
    finalAlignment->highestNominalScore = highestNominalScore;
    finalAlignment->alignment = alignment;
    finalAlignment->description = NULL;
    return finalAlignment;
}

// Free memory used to store alignments
void alignments_free_multi2(int queryNum) {
    struct alignment *alignment;
    struct finalAlignment *finalAlignment;

    // For each final alignment, free the description
    memSingleBlock_resetCurrent(alignments_finalAlignments_multi[queryNum]);

    //while ((finalAlignment =
                //memSingleBlock_getCurrent(alignments_finalAlignments_multi[queryNum])) != NULL) {
        //free(finalAlignment->description);
    //}

    // Free final alignments
    memSingleBlock_free(alignments_finalAlignments_multi[queryNum]);
}
