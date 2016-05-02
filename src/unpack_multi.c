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
#include "blast.h"

struct memBlocks *unpack_unpackRegions_multi[BATCH_SIZE];
struct memBlocks *unpack_subjectRegions_multi[BATCH_SIZE];

// Initialize region copying/unpacking
void unpack_initialize_multi(int queryNum) {
  unpack_unpackRegions_multi[queryNum] = memBlocks_initialize(
      sizeof(struct unpackRegion), constants_initialAllocUnpackRegions);

  unpack_subjectRegions_multi[queryNum] = memBlocks_initialize(
      sizeof(struct unpackRegion), constants_initialAllocUnpackRegions);
}

// Unpack entire or sections of a subject sequence before gapped alignment
void unpack_unpackSubject_multi(struct PSSMatrix PSSMatrix,
                                struct alignment *alignment, int queryNum) {
  unsigned char *subject, *unpackedSubject, wildcard, *edits, *endEdits;
  uint4 wildcardPosition;
  struct unpackRegion *firstRegion = NULL, *lastRegion, *currentRegion,
                      *unpackRegion;
  int4 regionStart, regionEnd, numRegions;

  ASSERT(encoding_alphabetType == encoding_protein);
  // No need to unpack a protein subject, or already unpacked nucleotide subject
  if (parameters_ssearch || encoding_alphabetType == encoding_protein) {
    // Just create a single region covering the entire sequence
    firstRegion = memBlocks_newEntry(unpack_unpackRegions_multi[queryNum]);
    firstRegion->startOffset = 0;
    firstRegion->endOffset = alignment->subjectLength;
    firstRegion->subject = alignment->subject;
    firstRegion->unpackedSubject = alignment->subject;
    firstRegion->subjectLength = alignment->subjectLength;
    alignment->unpackRegions = firstRegion;
    alignment->numUnpackRegions = 1;
    return;
  }

  // Get the subject regions for this alignment
  numRegions = unpack_getRegions(PSSMatrix, alignment, 0,
                                 unpack_unpackRegions_multi[queryNum]);
  lastRegion = memBlocks_getLastEntry(unpack_unpackRegions_multi[queryNum]);
  lastRegion++;
  firstRegion = lastRegion - numRegions;

  // Sort the regions in order of start position
  qsort(firstRegion, lastRegion - firstRegion, sizeof(struct unpackRegion),
        unpack_compareUnpackRegions);

  // Unpack each region
  currentRegion = firstRegion;
  while (currentRegion < lastRegion) {
    regionEnd = currentRegion->endOffset;
    regionStart = currentRegion->startOffset;

#ifdef VERBOSE
    if (parameters_verboseDloc == alignment->descriptionLocation) {
      printf("Unpack subject region %d to %d (length=%d)\n", regionStart,
             regionEnd, alignment->subjectLength);
      fflush(stdout);
    }
#endif

    // Get the subject region to be unpacked
    if (alignment->unpackRegions == NULL) {
      subject = alignment->subject;
    } else {
      unpackRegion = unpack_selectRegion(
          alignment->unpackRegions, alignment->numUnpackRegions, regionStart);
      subject = unpackRegion->subject;
    }

    // Declare memory for the region
    unpackedSubject = (unsigned char *)global_malloc(sizeof(char) *
                                                     (regionEnd - regionStart));

    // Unpack the region of interest
    encoding_byteUnpackRegion(unpackedSubject, subject + (regionStart / 4),
                              regionEnd - regionStart);
    unpackedSubject -= regionStart;
    currentRegion->unpackedSubject = unpackedSubject;

    currentRegion->subject = subject;
    currentRegion->subjectLength = alignment->subjectLength;

    blast_totalUnpacked += (regionEnd - regionStart);

    currentRegion++;
  }

  currentRegion = firstRegion;

  // Get wildcard edits for the sequence
  edits = alignment->edits;
  endEdits = alignment->edits + alignment->encodedLength -
             ((alignment->subjectLength + 3) / 4);

  // If there are edits
  if (edits < endEdits) {
    // Read first wildcard
    wildcard = *edits;
    edits++;

    // Read its position
    vbyte_getVbyte(edits, &wildcardPosition);

    // For each region in order of position in the subject
    while (currentRegion < lastRegion) {
      // Skip past edits that are before current region
      while (edits < endEdits &&
             wildcardPosition < currentRegion->startOffset) {
        // Read wildcard
        wildcard = *edits;
        edits++;

        // Read its position
        vbyte_getVbyte(edits, &wildcardPosition);
      }

      // Process edits that are in the current region
      while (edits < endEdits && wildcardPosition < currentRegion->endOffset) {
        // Insert wildcard into sequence
        currentRegion->unpackedSubject[wildcardPosition] = wildcard;

        // Read next wildcard
        wildcard = *edits;
        edits++;

        // Read its position
        vbyte_getVbyte(edits, &wildcardPosition);
      }

      // Advance to the next region
      currentRegion++;
    }
  }

  alignment->unpackRegions = firstRegion;
  alignment->numUnpackRegions = lastRegion - firstRegion;
}

// Free memory used to store unpacked regions
void unpack_free_multi(int queryNum) {
  struct unpackRegion *region;

  // For each unpack region
  if (!parameters_ssearch && encoding_alphabetType == encoding_nucleotide) {
    memBlocks_resetCurrent(unpack_unpackRegions_multi[queryNum]);
    while ((region = memBlocks_getCurrent(
                unpack_unpackRegions_multi[queryNum])) != NULL) {
      // Free the unpacked sequence
      free(region->unpackedSubject + region->startOffset);
    }
  }

  // For each copied subject region
  memBlocks_resetCurrent(unpack_subjectRegions_multi[queryNum]);
  while ((region = memBlocks_getCurrent(
              unpack_subjectRegions_multi[queryNum])) != NULL) {
    // Free the subject
    free(region->subject + region->startOffset / 4);
  }

  memBlocks_free(unpack_unpackRegions_multi[queryNum]);
  memBlocks_free(unpack_subjectRegions_multi[queryNum]);
}

// Load a single subject into memory
int4 unpack_loadSubject_multi(struct PSSMatrix PSSMatrix,
                              struct alignment *alignment, int queryNum) {
  uint4 totalCopied = 0;
  unsigned char *subject, *edits, *endEdits;
  struct unpackRegion *firstRegion = NULL, *lastRegion, *currentRegion;
  int4 numRegions, regionStart, regionEnd;

  // If protein search
  if (encoding_alphabetType == encoding_protein) {
    // Make copy of sequence
    subject = (unsigned char *)global_malloc(sizeof(unsigned char) *
                                             alignment->encodedLength);
    subject++;
    memcpy(subject - 1, alignment->subject - 1, alignment->encodedLength);
    alignment->subject = subject;

    blast_totalCopied += alignment->encodedLength;
  }
  // If a nucleotide search
  else {
    // Get a list of regions to copy
    numRegions = unpack_getRegions(PSSMatrix, alignment, 1,
                                   unpack_subjectRegions_multi[queryNum]);
    lastRegion = memBlocks_getLastEntry(unpack_subjectRegions_multi[queryNum]);
    lastRegion++;
    firstRegion = lastRegion - numRegions;

#ifdef VERBOSE
    if (parameters_verboseDloc == alignment->descriptionLocation) {
      printf("%d regions for subject\n", lastRegion - firstRegion);
      fflush(stdout);
    }
#endif

    // Copy each region into memory
    currentRegion = firstRegion;
    while (currentRegion < lastRegion) {
#ifdef VERBOSE
      if (parameters_verboseDloc == alignment->descriptionLocation) {
        printf("Load region %d to %d into memory\n", currentRegion->startOffset,
               currentRegion->endOffset);
        fflush(stdout);
        fflush(stdout);
      }
#endif

      regionStart = currentRegion->startOffset / 4;
      regionEnd = (currentRegion->endOffset + 3) / 4;

      currentRegion->unpackedSubject = NULL;
      currentRegion->subject = (unsigned char *)global_malloc(
          sizeof(unsigned char) * (regionEnd - regionStart));

      totalCopied += regionEnd - regionStart;
      memcpy(currentRegion->subject, alignment->subject + regionStart,
             regionEnd - regionStart);
      currentRegion->subject -= regionStart;
      currentRegion->subjectLength = alignment->subjectLength;

      blast_totalCopied += (regionEnd - regionStart);

      currentRegion++;
    }

    // Store new alignment regions
    alignment->unpackRegions = firstRegion;
    alignment->numUnpackRegions = lastRegion - firstRegion;

    // If there are edits for this subject
    if (alignment->edits != NULL) {
      edits = alignment->edits;
      endEdits = alignment->subject + alignment->encodedLength;

      // Make an in-memory copy of them
      alignment->edits =
          (unsigned char *)malloc(sizeof(char) * (endEdits - edits));
      memcpy(alignment->edits, edits, endEdits - edits);
    }

    alignment->subject = NULL;
  }

  alignment->inMemorySubject = 1;

  return totalCopied;
}
