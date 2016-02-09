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

int4 *gappedExtension_matchRow_multi[BATCH_SIZE] = { NULL };
int4 *gappedExtension_insertQrow_multi[BATCH_SIZE] = { NULL };
int4 *gappedExtension_insertSrow_multi[BATCH_SIZE] = { NULL };

unsigned char **gappedExtension_traceback_multi[BATCH_SIZE] = { NULL };
unsigned char *gappedExtension_tracebackData_multi[BATCH_SIZE] = { NULL };
int4 gappedExtension_rowSizes_multi[BATCH_SIZE] = { 0 };
int4 gappedExtension_numRows_multi[BATCH_SIZE] = { 0 };
int4 gappedExtension_tracebackAlloc_multi[BATCH_SIZE] = { 0 };

// Prototypes
struct dpResults gappedExtension_dpBeforeSeed_multi(
    struct PSSMatrix PSSMatrix, int4 dropoff, struct coordinate seed,
    struct unpackRegion *unpackRegion, int queryNum);
struct dpResults gappedExtension_dpAfterSeed_multi(
    struct PSSMatrix PSSMatrix, int4 dropoff, struct unpackRegion *unpackRegion,
    int4 subjectLength, int4 seedSubjectOffset, int queryNum);
struct trace
gappedExtension_traceBeforeSeed_multi(struct dpResults beforeDpResults,
                                      struct coordinate seed, int queryNum);
struct trace
gappedExtension_traceAfterSeed_multi(struct dpResults beforeDpResults,
                                     int4 queryLength, int queryNum);

// Build a gapped extension with a trace and nominal score from the seed point
// of an ungapped
// extension using dynamic programming
struct gappedExtension *gappedExtension_build_multi(
    struct ungappedExtension *ungappedExtension, struct PSSMatrix PSSMatrix,
    int4 subjectSize, unsigned char *subject, struct unpackRegion *unpackRegion,
    int4 dropoff, int queryNum) {
  struct coordinate seed;
  unsigned char *choppedSubject;
  struct dpResults beforeDpResults, afterDpResults;
  struct trace beforeTrace, afterTrace, trace;
  struct PSSMatrix choppedPSSMatrix;
  int4 choppedSubjectSize;
  struct gappedExtension *gappedExtension;
  int4 strandOffset = 0;

  // Perform dynamic programming for points before the seed
  seed = ungappedExtension->seed;
  if (seed.queryOffset > PSSMatrix.strandLength) {
    // If query position is in the second strand, remove first strand from PSSM
    strandOffset = PSSMatrix.strandLength;
    seed.queryOffset -= PSSMatrix.strandLength;
    PSSMatrix = PSSMatrix_chop(PSSMatrix, PSSMatrix.strandLength);
  } else {
    // Otherwise remove second strand
    PSSMatrix.length = PSSMatrix.strandLength;
  }

  beforeDpResults = gappedExtension_dpBeforeSeed_multi(PSSMatrix, dropoff, seed,
                                                       unpackRegion, queryNum);

  // Trace back and create the trace which specifies the first half of the
  // alignment
  beforeTrace =
      gappedExtension_traceBeforeSeed_multi(beforeDpResults, seed, queryNum);

  // Chop the start off the query and subject so they begin at the seed
  choppedPSSMatrix = PSSMatrix_chop(PSSMatrix, seed.queryOffset);
  choppedSubject = subject + seed.subjectOffset;
  choppedSubjectSize = subjectSize - (seed.subjectOffset);

  // Perform dynamic programming for points after the seed
  afterDpResults = gappedExtension_dpAfterSeed_multi(
      choppedPSSMatrix, dropoff, unpackRegion, choppedSubjectSize,
      seed.subjectOffset, queryNum);

  // Trace back to get the trace for the seed onwards
  afterTrace = gappedExtension_traceAfterSeed_multi(
      afterDpResults, choppedPSSMatrix.length, queryNum);

  // Join afterTrace to the end of beforeTrace
  trace = gappedExtension_joinTraces(beforeTrace, afterTrace);
  free(afterTrace.traceCodes);

  // Adjust coordinates if extension was performed in the second strand
  afterDpResults.best.queryOffset += strandOffset;
  beforeDpResults.best.queryOffset += strandOffset;
  trace.queryStart += strandOffset;

  // Create gapped extension
  gappedExtension =
      (struct gappedExtension *)global_malloc(sizeof(struct gappedExtension));
  gappedExtension->trace = trace;
  gappedExtension->next = NULL;

  // Start of afterTrace is end of the gapped extension, but we need to add seed
  // position
  // to get correct offset
  gappedExtension->queryEnd =
      seed.queryOffset + afterTrace.queryStart + strandOffset;
  gappedExtension->subjectEnd = seed.subjectOffset + afterTrace.subjectStart;

  //	if (dloc == 88197331)
  //		printf("final[%d,%d,%d](%d)\n", beforeDpResults.bestScore,
  // afterDpResults.bestScore,
  //		choppedPSSMatrix.matrix[0][unpackRegion->unpackedSubject[seed.subjectOffset]],
  // seed.queryOffset);

  // Determine score by combining score from the two traces, and the match score
  // at
  // the seed position
  gappedExtension->nominalScore =
      beforeDpResults.bestScore + afterDpResults.bestScore +
      choppedPSSMatrix.matrix
          [0][unpackRegion->unpackedSubject[seed.subjectOffset]];

  // Update ungappedExtension start/end
  ungappedExtension->start.queryOffset = trace.queryStart;
  ungappedExtension->end.queryOffset = gappedExtension->queryEnd;
  ungappedExtension->start.subjectOffset = trace.subjectStart;
  ungappedExtension->end.subjectOffset = gappedExtension->subjectEnd;
  ungappedExtension->nominalScore = gappedExtension->nominalScore;

#ifdef VERBOSE
  if (parameters_verboseDloc == blast_dloc) {
    printf("Gapped Extension from %d,%d to %d,%d score %d\n", trace.queryStart,
           trace.subjectStart, gappedExtension->queryEnd,
           gappedExtension->subjectEnd, gappedExtension->nominalScore);
  }
#endif

  return gappedExtension;
}

// Given the results of dynamic programming (a matrix of trace codes and a
// highest scoring position in
// the matrix) for finding the START of the alignment, performs the simple
// operation of finding the path
// from the highest scoring point back to the seed
struct trace
gappedExtension_traceBeforeSeed_multi(struct dpResults beforeDpResults,
                                      struct coordinate seed, int queryNum) {
  int4 queryPosition, subjectPosition;
  unsigned char **traceback;
  unsigned char traceCode;
  unsigned char state = 0;
  struct trace trace;
  unsigned char *traceCodes;
  uint4 traceCount = 0;

  traceback = beforeDpResults.traceback;
  trace.queryStart = queryPosition = beforeDpResults.best.queryOffset;
  trace.subjectStart = subjectPosition = beforeDpResults.best.subjectOffset;

  // Declare memory for tracecodes; for maximum possible number of codes that
  // could
  // be generated by this trace
  traceCodes = (unsigned char *)global_malloc(
      sizeof(unsigned char) * (seed.queryOffset - queryPosition +
                               seed.subjectOffset - subjectPosition));

  while (queryPosition < seed.queryOffset &&
         subjectPosition < seed.subjectOffset) {
    // Construct the trace
    traceCodes[traceCount] = state;
    traceCount++;

    //        printf("(%p)", traceback[queryPosition]);
    //        printf("(%d,%d)", queryPosition, subjectPosition); fflush(stdout);
    traceCode = traceback[queryPosition][subjectPosition];

    // If we got to current cell through a MATCH
    if (state == 0) {
      // Move to cell before this one
      queryPosition++;
      subjectPosition++;

      // We are only interested in lowest 2 bits of tracecode
      traceCode = traceCode << 6;
      traceCode = traceCode >> 6;

      // Tracecode determines if we matched or inserted here
      state = traceCode;
    }
    // If we got to current cell through an Ix
    else if (state == 1) {
      // Move to cell before this one
      subjectPosition++;

      // We are int4erest in bits 3rd and 4th from right
      traceCode = traceCode << 4;
      traceCode = traceCode >> 6;

      // Tracecode determines if we matched or inserted here
      state = traceCode;
    }
    // If we got to current cell through an Iy
    else if (state == 2) {
      // Move to cell before this one
      queryPosition++;

      // We are int4erest in bits 5th and 6th from right
      traceCode = traceCode << 2;
      traceCode = traceCode >> 6;

      // Tracecode determines if we matched or inserted here
      state = traceCode;
    }
  }

  // End trace with insertions needed to get us back to the seed
  // (most likely none will be required)
  while (queryPosition < seed.queryOffset) {
    traceCodes[traceCount] = 2;
    traceCount++;
    queryPosition++;
  }

  while (subjectPosition < seed.subjectOffset) {
    traceCodes[traceCount] = 1;
    traceCount++;
    subjectPosition++;
  }

  trace.traceCodes = traceCodes;
  trace.length = traceCount;

  return trace;
}

// Given the results of dynamic programming (a matrix of trace codes and a
// highest scoring position in
// the matrix) for finding the END of the alignment, performs the simple
// operation of finding the path
// from the highest scoring point back to the seed
struct trace
gappedExtension_traceAfterSeed_multi(struct dpResults beforeDpResults,
                                     int4 queryLength, int queryNum) {
  int4 queryPosition, subjectPosition;
  unsigned char **traceback;
  unsigned char traceCode;
  unsigned char state = 0;
  struct trace trace;
  unsigned char *traceCodes;
  uint4 traceCount = 0;

  traceback = beforeDpResults.traceback;
  trace.queryStart = 0;
  trace.subjectStart = 0;

  // Start at the end of the alignment
  queryPosition = beforeDpResults.best.queryOffset;
  subjectPosition = beforeDpResults.best.subjectOffset;

  // Declare memory for tracecodes; for maximum possible number of codes that
  // could
  // be generated by this trace
  traceCodes = (unsigned char *)global_malloc(
      sizeof(unsigned char) * (queryPosition + subjectPosition));

  while (queryPosition > 0 && subjectPosition > 0) {
    traceCode = traceback[queryPosition][subjectPosition];
    // If we got to current cell through a MATCH
    if (state == 0) {
      // Move to cell before this one
      queryPosition--;
      subjectPosition--;

      // We are only interested in lowest 2 bits of tracecode
      traceCode = traceCode << 6;
      traceCode = traceCode >> 6;

      // Tracecode determines if we matched or inserted here
      state = traceCode;
    }
    // If we got to current cell through an Ix
    else if (state == 1) {
      // Move to cell before this one
      subjectPosition--;

      // We are int4erest in bits 3rd and 4th from right
      traceCode = traceCode << 4;
      traceCode = traceCode >> 6;

      // Tracecode determines if we matched or inserted here
      state = traceCode;
    }
    // If we got to current cell through an Iy
    else if (state == 2) {
      // Move to cell before this one
      queryPosition--;

      // We are int4erest in bits 5th and 6th from right
      traceCode = traceCode << 2;
      traceCode = traceCode >> 6;

      // Tracecode determines if we matched or inserted here
      state = traceCode;
    }

    // Construct the trace
    traceCodes[traceCount] = state;
    traceCount++;
  }

  // End trace with insertions needed to get us back to the seed
  // (most likely none will be required)
  while (queryPosition > 0) {
    traceCodes[traceCount] = 2;
    traceCount++;
    queryPosition--;
  }

  while (subjectPosition > 0) {
    traceCodes[traceCount] = 1;
    traceCount++;
    subjectPosition--;
  }

  trace.traceCodes = traceCodes;
  trace.length = traceCount;
  trace.queryStart = beforeDpResults.best.queryOffset;
  trace.subjectStart = beforeDpResults.best.subjectOffset;

  return trace;
}

// Perform dynamic programming to explore possible start points and alignments
// that end at
// the given seed
struct dpResults gappedExtension_dpBeforeSeed_multi(
    struct PSSMatrix PSSMatrix, int4 dropoff, struct coordinate seed,
    struct unpackRegion *unpackRegion, int queryNum) {
  int2 **queryPosition, **bestQueryPosition;
  int2 *matrixColumn;
  unsigned char *subject, *rowDropoff, *columnDropoff, *minRowDropoff;
  unsigned char *subjectPosition, *bestSubjectPosition;
  unsigned char **tracebackRow, *tracebackColumn;
  int4 bestScore = 0, dropoffThreshold;
  int4 *matchRow, *insertQrow, *insertSrow, rowOffset;
  int4 queryDistance, subjectDistance;
  int4 oldMatch, match, previousOldMatch, previousOldInsertS,
      previousOldInsertQ, previousInsertQ = 0;
  int4 previousMatch, previousInsertS;
  struct dpResults dpResults;
  unsigned char rightOfDropoff;
  unsigned char *tracebackDataPtr, *newTracebackData;

  subject = unpackRegion->unpackedSubject;

  if (gappedExtension_matchRow_multi[queryNum] == NULL) {
    // Declare processing rows for storing match, insert-subject and
    // insert-query values
    gappedExtension_rowSizes_multi[queryNum] =
        2 * statistics_gappedFinalNominalDropoff_multi[queryNum] /
        parameters_extendGap;

    // Malloc new rows
    gappedExtension_matchRow_multi[queryNum] = (int4 *)global_malloc(
        sizeof(int4) * gappedExtension_rowSizes_multi[queryNum]);
    gappedExtension_insertQrow_multi[queryNum] = (int4 *)global_malloc(
        sizeof(int4) * gappedExtension_rowSizes_multi[queryNum]);
    gappedExtension_insertSrow_multi[queryNum] = (int4 *)global_malloc(
        sizeof(int4) * gappedExtension_rowSizes_multi[queryNum]);
  }

  // Determine lowest score before dropoff
  dropoffThreshold = -dropoff;

  // If more rows are required
  if (seed.queryOffset > gappedExtension_numRows_multi[queryNum]) {
    // Increase number of row pointers
    gappedExtension_numRows_multi[queryNum] = seed.queryOffset;
    gappedExtension_traceback_multi[queryNum] =
        (unsigned char **)global_realloc(
            gappedExtension_traceback_multi[queryNum],
            sizeof(unsigned char *) *
                (gappedExtension_numRows_multi[queryNum]));
  }

  // If first time, initialize traceback rows
  if (gappedExtension_tracebackData_multi[queryNum] == NULL) {
    gappedExtension_tracebackAlloc_multi[queryNum] =
        constants_initialTracebackAlloc;
    gappedExtension_tracebackData_multi[queryNum] = global_malloc(
        sizeof(char) * gappedExtension_tracebackAlloc_multi[queryNum]);
  }

  tracebackDataPtr = gappedExtension_tracebackData_multi[queryNum] +
                     gappedExtension_tracebackAlloc_multi[queryNum] - 1;

  bestSubjectPosition = subjectPosition = subject + seed.subjectOffset - 1;
  bestQueryPosition = queryPosition = PSSMatrix.matrix + seed.queryOffset - 1;

  // Initialize row pointers
  rowOffset = (subjectPosition - subject);
  matchRow = gappedExtension_matchRow_multi[queryNum] +
             (rowOffset % gappedExtension_rowSizes_multi[queryNum]);
  insertQrow = gappedExtension_insertQrow_multi[queryNum] +
               (rowOffset % gappedExtension_rowSizes_multi[queryNum]);
  insertSrow = gappedExtension_insertSrow_multi[queryNum] +
               (rowOffset % gappedExtension_rowSizes_multi[queryNum]);

  // Set initial row dropoff and column dropoff
  rowDropoff = subject;
  columnDropoff = subject + seed.subjectOffset;

  // Calculate minimum row dropoff after this row is processed
  minRowDropoff = columnDropoff - ((dropoff / parameters_extendGap) + 2);
  if (minRowDropoff < subject)
    minRowDropoff = subject;

  // Unpack more of the subject if required
  unpack_extendRegionStart(minRowDropoff - subject, unpackRegion);

  // If the unpacked subject has moved in memory
  if (subject != unpackRegion->unpackedSubject) {
    // Update all pointers to point at the new subject
    subjectPosition =
        (subjectPosition - subject) + unpackRegion->unpackedSubject;
    rowDropoff = (rowDropoff - subject) + unpackRegion->unpackedSubject;
    columnDropoff = (columnDropoff - subject) + unpackRegion->unpackedSubject;
    bestSubjectPosition =
        (bestSubjectPosition - subject) + unpackRegion->unpackedSubject;
    subject = unpackRegion->unpackedSubject;
  }

  // Initialize traceback pointer for this row
  tracebackRow = gappedExtension_traceback_multi[queryNum] +
                 (queryPosition - PSSMatrix.matrix);
  *tracebackRow = tracebackDataPtr - (subjectPosition - subject);
  tracebackColumn = tracebackDataPtr;
  tracebackDataPtr -= (subjectPosition - minRowDropoff + 1);

  // -----FIRST ROW-----

  // Using first column of query matrix
  matrixColumn = *queryPosition;

  // -----FIRST CELL-----
  // Set M value for bottom-right cell
  match = matrixColumn[*subjectPosition];
  *matchRow = match;

  // Set DUMMY Ix and Iy values, which should never be used
  *insertSrow = constants_gappedExtensionDummyValue;
  *insertQrow = constants_gappedExtensionDummyValue;

  // M came from M
  *tracebackColumn = 0;

  // If this is the best-yet scoring cell
  if (match > bestScore) {
    // Update best start cell data
    bestScore = *matchRow;
    dropoffThreshold = bestScore - dropoff;
    bestQueryPosition = queryPosition;
    bestSubjectPosition = subjectPosition;
  }

  // Record match and insertS for this about-to-be-previous cell
  previousMatch = match;
  previousInsertS = *insertSrow;

  subjectDistance = 0;
  subjectPosition--;
  matchRow--;
  insertSrow--;
  insertQrow--;
  tracebackColumn--;

  // Check for scoring row wrap-around
  if (matchRow < gappedExtension_matchRow_multi[queryNum]) {
    matchRow += gappedExtension_rowSizes_multi[queryNum];
    insertSrow += gappedExtension_rowSizes_multi[queryNum];
    insertQrow += gappedExtension_rowSizes_multi[queryNum];
  }

  // ----- REMAINING CELLS -----
  // For each remaining column in the bottom row, scanning from right-to-left
  while (subjectPosition >= subject) {
    // Set value for M
    match = matrixColumn[*subjectPosition] - parameters_openGap -
            subjectDistance * parameters_extendGap;
    *matchRow = match;

    // Set value for Ix
    if (previousInsertS - parameters_extendGap >
        previousMatch - parameters_openGap) {
      *insertSrow = previousInsertS - parameters_extendGap;
      // M came from Ix and Ix came from Ix
      *tracebackColumn = 5;
    } else {
      *insertSrow = previousMatch - parameters_openGap;
      // M came from Ix and Ix came from M
      *tracebackColumn = 1;
    }

    // Set DUMMY Iy value, which should never be used
    *insertQrow = constants_gappedExtensionDummyValue;

    // If this is the best-yet scoring cell
    if (match > bestScore) {
      // Update best start cell data
      bestScore = match;
      dropoffThreshold = bestScore - dropoff;
      bestQueryPosition = queryPosition;
      bestSubjectPosition = subjectPosition;
    }

    // If score at current cell is below dropoff
    if (dropoffThreshold > match && dropoffThreshold > *insertSrow) {
      // Record dropoff position
      rowDropoff = subjectPosition;
      // And stop processing row
      break;
    }

    // Record match and insertS for this about-to-be-previous cell
    previousMatch = match;
    previousInsertS = *insertSrow;

    subjectPosition--;
    matchRow--;
    insertSrow--;
    insertQrow--;
    tracebackColumn--;

    // Check for scoring row wrap-around
    if (matchRow < gappedExtension_matchRow_multi[queryNum]) {
      matchRow += gappedExtension_rowSizes_multi[queryNum];
      insertSrow += gappedExtension_rowSizes_multi[queryNum];
      insertQrow += gappedExtension_rowSizes_multi[queryNum];
    }

    subjectDistance++;
  }

  //	print(gappedExtension_matchRow_multi[queryNum], subject, rowDropoff,
  // columnDropoff);

  queryDistance = 0;

  // -----REMAINING ROWS-----
  while (queryPosition > PSSMatrix.matrix && rowDropoff < columnDropoff) {
    queryPosition--;
    tracebackRow--;
    subjectPosition = columnDropoff - 1;

    // Calculate minimum row dropoff after this row is processed
    minRowDropoff =
        rowDropoff - ((PSSMatrix.highestValue / parameters_extendGap) + 2);
    if (minRowDropoff < subject)
      minRowDropoff = subject;

    // If not enough space in traceback data for this row, realloc
    if (subjectPosition - minRowDropoff >=
        tracebackDataPtr - gappedExtension_tracebackData_multi[queryNum]) {
      newTracebackData = (unsigned char *)global_malloc(
          sizeof(char) * gappedExtension_tracebackAlloc_multi[queryNum] * 2);

      // Move existing data to end of new block
      memcpy(newTracebackData + gappedExtension_tracebackAlloc_multi[queryNum],
             gappedExtension_tracebackData_multi[queryNum],
             gappedExtension_tracebackAlloc_multi[queryNum]);

      // Update traceback data pointer
      tracebackDataPtr =
          newTracebackData + gappedExtension_tracebackAlloc_multi[queryNum] +
          (tracebackDataPtr - gappedExtension_tracebackData_multi[queryNum]);

      // Update existing rows to point to new data block
      while (tracebackRow <
             gappedExtension_traceback_multi[queryNum] + seed.queryOffset - 1) {
        tracebackRow++;
        *tracebackRow =
            newTracebackData + gappedExtension_tracebackAlloc_multi[queryNum] +
            (*tracebackRow - gappedExtension_tracebackData_multi[queryNum]);
      }

      // Data block is now double the size
      gappedExtension_tracebackAlloc_multi[queryNum] *= 2;

      free(gappedExtension_tracebackData_multi[queryNum]);
      gappedExtension_tracebackData_multi[queryNum] = newTracebackData;
    }

    // Initialize traceback pointer for this row
    tracebackRow = gappedExtension_traceback_multi[queryNum] +
                   (queryPosition - PSSMatrix.matrix);
    *tracebackRow = tracebackDataPtr - (subjectPosition - subject);
    tracebackColumn = tracebackDataPtr;
    tracebackDataPtr -= (subjectPosition - minRowDropoff + 1);

    // Unpack more of the subject if required
    unpack_extendRegionStart(minRowDropoff - subject, unpackRegion);

    // If the unpacked subject has moved in memory
    if (subject != unpackRegion->unpackedSubject) {
      // Update all pointers to point at the new subject
      subjectPosition =
          (subjectPosition - subject) + unpackRegion->unpackedSubject;
      rowDropoff = (rowDropoff - subject) + unpackRegion->unpackedSubject;
      columnDropoff = (columnDropoff - subject) + unpackRegion->unpackedSubject;
      bestSubjectPosition =
          (bestSubjectPosition - subject) + unpackRegion->unpackedSubject;
      subject = unpackRegion->unpackedSubject;
    }

    // Reset row pointers to start of rows
    rowOffset = (subjectPosition - subject);
    matchRow = gappedExtension_matchRow_multi[queryNum] +
               (rowOffset % gappedExtension_rowSizes_multi[queryNum]);
    insertQrow = gappedExtension_insertQrow_multi[queryNum] +
                 (rowOffset % gappedExtension_rowSizes_multi[queryNum]);
    insertSrow = gappedExtension_insertSrow_multi[queryNum] +
                 (rowOffset % gappedExtension_rowSizes_multi[queryNum]);

    // Using next column of query matrix
    matrixColumn = *queryPosition;

    // -----FAR RIGHT CELL-----
    // Record some old values
    previousOldMatch = *matchRow;
    previousOldInsertQ = *insertQrow;
    previousOldInsertS = *insertSrow;

    // Set Iy value
    if (*insertQrow - parameters_extendGap > *matchRow - parameters_openGap) {
      *insertQrow = *insertQrow - parameters_extendGap;
      // Iy is derived from Iy, M is derived from Iy
      *tracebackColumn = 34;
    } else {
      *insertQrow = *matchRow - parameters_openGap;
      // Iy is derived from M, M is derived from Iy
      *tracebackColumn = 2;
    }

    // Set DUMMY values for M and Iy, which should never be used
    match = *matchRow = constants_gappedExtensionDummyValue;
    *insertSrow = constants_gappedExtensionDummyValue;

    // If score at current cell is below dropoff
    if (dropoffThreshold > *insertQrow) {
      // Record dropoff position
      columnDropoff = subjectPosition;
      rightOfDropoff = 1;
    } else {
      // We are left of the column dropoff for this row
      rightOfDropoff = 0;
    }

    // Record match and insertS for this about-to-be-previous cell
    previousMatch = match;
    previousInsertS = *insertSrow;

    subjectPosition--;
    matchRow--;
    insertSrow--;
    insertQrow--;
    tracebackColumn--;

    // Check for scoring row wrap-around
    if (matchRow < gappedExtension_matchRow_multi[queryNum]) {
      matchRow += gappedExtension_rowSizes_multi[queryNum];
      insertSrow += gappedExtension_rowSizes_multi[queryNum];
      insertQrow += gappedExtension_rowSizes_multi[queryNum];
    }

    // -----CELLS RIGHT OF ROW DROPOFF-----
    while (subjectPosition >= rowDropoff) {
      // Remember old M value (for cell below this one)
      oldMatch = *matchRow;

      // Calculate new M value
      if (previousOldMatch >= previousOldInsertQ) {
        if (previousOldMatch >= previousOldInsertS) {
          match = matrixColumn[*subjectPosition] + previousOldMatch;
          // M is derived from M
          *tracebackColumn = 0;
        } else {
          match = matrixColumn[*subjectPosition] + previousOldInsertS;
          // M is derived from Ix
          *tracebackColumn = 1;
        }
      } else {
        if (previousOldInsertQ >= previousOldInsertS) {
          match = matrixColumn[*subjectPosition] + previousOldInsertQ;
          // M is derived from Iy
          *tracebackColumn = 2;
        } else {
          match = matrixColumn[*subjectPosition] + previousOldInsertS;
          // M is derived from Ix
          *tracebackColumn = 1;
        }
      }

      *matchRow = match;

      // Record some old values
      previousOldMatch = oldMatch;
      previousOldInsertQ = *insertQrow;
      previousOldInsertS = *insertSrow;

      // Set new Iy value
      if (oldMatch - parameters_openGap >= *insertQrow - parameters_extendGap) {
        *insertQrow = oldMatch - parameters_openGap;
        // Iy is derived from M
        // No change to traceback
      } else {
        *insertQrow = *insertQrow - parameters_extendGap;
        // Iy is derived from Iy
        *tracebackColumn |= 32;
      }
      // Calculate new Ix
      if (previousMatch - parameters_openGap >=
          previousInsertS - parameters_extendGap) {
        *insertSrow = previousMatch - parameters_openGap;
        // Ix is derived from M
        // No change to traceback
      } else {
        *insertSrow = previousInsertS - parameters_extendGap;
        // Ix is derived from Ix
        *tracebackColumn |= 4;
      }

      // If this is the best-yet scoring cell
      if (match > bestScore) {
        // Update best start cell data
        bestScore = match;
        dropoffThreshold = bestScore - dropoff;
        bestQueryPosition = queryPosition;
        bestSubjectPosition = subjectPosition;
      }

      // If score at current cell (and cells to its right) are below dropoff
      if (rightOfDropoff) {
        if (dropoffThreshold > match && dropoffThreshold > *insertSrow &&
            dropoffThreshold > *insertQrow) {
          // Record dropoff position
          columnDropoff = subjectPosition;
        } else {
          // We are left of the column dropoff for this row
          rightOfDropoff = 0;
        }
      }

      // Record match and insertS for this about-to-be-previous cell
      previousMatch = match;
      previousInsertS = *insertSrow;
      previousInsertQ = *insertQrow;

      subjectPosition--;
      matchRow--;
      insertSrow--;
      insertQrow--;
      tracebackColumn--;

      // Check for scoring row wrap-around
      if (matchRow < gappedExtension_matchRow_multi[queryNum]) {
        matchRow += gappedExtension_rowSizes_multi[queryNum];
        insertSrow += gappedExtension_rowSizes_multi[queryNum];
        insertQrow += gappedExtension_rowSizes_multi[queryNum];
      }
    }

    // -----CELLS LEFT OF ROW DROPOFF -----
    if (!(dropoffThreshold > previousMatch &&
          dropoffThreshold > previousInsertS &&
          dropoffThreshold > previousInsertQ)) {
      while (subjectPosition >= subject) {
        // Set value for Ix
        *insertSrow = previousInsertS - parameters_extendGap;
        // Ix came from Ix
        *tracebackColumn = 4;

        // Set DUMMY values for M and Ix, which should never be used
        *matchRow = constants_gappedExtensionDummyValue;
        *insertQrow = constants_gappedExtensionDummyValue;

        // If score at current cell is below dropoff
        if (dropoffThreshold > *insertSrow) {
          // Stop processing row
          subjectPosition--;
          break;
        }

        // Record match and insertS for this about-to-be-previous cell
        previousInsertS = *insertSrow;

        subjectPosition--;
        matchRow--;
        insertSrow--;
        insertQrow--;
        tracebackColumn--;

        // Check for scoring row wrap-around
        if (matchRow < gappedExtension_matchRow_multi[queryNum]) {
          matchRow += gappedExtension_rowSizes_multi[queryNum];
          insertSrow += gappedExtension_rowSizes_multi[queryNum];
          insertQrow += gappedExtension_rowSizes_multi[queryNum];
        }

        subjectDistance++;
      }
    }

    // Record dropoff position
    rowDropoff = subjectPosition + 1;

    //		print(gappedExtension_matchRow_multi[queryNum], subject,
    //rowDropoff,
    // columnDropoff);

    queryDistance++;
  }

  dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
  dpResults.best.subjectOffset = bestSubjectPosition - subject;
  dpResults.bestScore = bestScore;
  dpResults.traceback = gappedExtension_traceback_multi[queryNum];

  return dpResults;
}

// Perform dynamic programming to explore possible END points and alignments
// that start at
// the given seed
struct dpResults gappedExtension_dpAfterSeed_multi(
    struct PSSMatrix PSSMatrix, int4 dropoff, struct unpackRegion *unpackRegion,
    int4 subjectLength, int4 seedSubjectOffset, int queryNum) {
  int2 **queryPosition, **bestQueryPosition, **queryEnd;
  int2 *matrixColumn;
  unsigned char *subject, *rowDropoff, *columnDropoff, *maxRowDropoff,
      *newSubject;
  unsigned char *subjectPosition, *bestSubjectPosition, *subjectEnd;
  unsigned char **tracebackRow, *tracebackColumn;
  int4 bestScore = 0, dropoffThreshold;
  int4 *matchRow, *insertQrow, *insertSrow, rowOffset, *endMatchRow;
  int4 queryDistance, subjectDistance;
  int4 oldMatch, match, previousOldMatch, previousOldInsertS,
      previousOldInsertQ;
  int4 previousMatch, previousInsertS, previousInsertQ = 0;
  struct dpResults dpResults;
  unsigned char leftOfDropoff;
  int4 queryLength;
  unsigned char *tracebackDataPtr, *newTracebackData;

  subject = unpackRegion->unpackedSubject + seedSubjectOffset;

  queryLength = PSSMatrix.length;
  subjectEnd = subject + subjectLength;
  queryEnd = PSSMatrix.matrix + queryLength;
  endMatchRow = gappedExtension_matchRow_multi[queryNum] +
                gappedExtension_rowSizes_multi[queryNum];

  // Determine lowest score before dropoff
  dropoffThreshold = -dropoff;

  // If more rows are required
  if (queryLength > gappedExtension_numRows_multi[queryNum]) {
    // Increase number of row pointers
    gappedExtension_numRows_multi[queryNum] = queryLength;
    gappedExtension_traceback_multi[queryNum] =
        (unsigned char **)global_realloc(
            gappedExtension_traceback_multi[queryNum],
            sizeof(unsigned char *) *
                (gappedExtension_numRows_multi[queryNum]));
  }

  tracebackDataPtr = gappedExtension_tracebackData_multi[queryNum];

  bestSubjectPosition = subjectPosition = subject + 1;
  bestQueryPosition = queryPosition = PSSMatrix.matrix + 1;

  // Initialize rows
  matchRow = gappedExtension_matchRow_multi[queryNum] + 1;
  insertQrow = gappedExtension_insertQrow_multi[queryNum] + 1;
  insertSrow = gappedExtension_insertSrow_multi[queryNum] + 1;

  // Set initial row dropoff and column dropoff
  rowDropoff = subject + subjectLength - 1;
  columnDropoff = subject;

  // Calculate maximum row dropoff after this row is processed
  maxRowDropoff = columnDropoff + ((dropoff / parameters_extendGap) + 2);
  if (maxRowDropoff > subjectEnd)
    maxRowDropoff = subjectEnd;

  // Initialize traceback pointer for this row
  tracebackRow = gappedExtension_traceback_multi[queryNum] +
                 (queryPosition - PSSMatrix.matrix);
  *tracebackRow = tracebackDataPtr - (subjectPosition - subject);
  tracebackColumn = tracebackDataPtr;
  tracebackDataPtr += (maxRowDropoff - subjectPosition + 1);

  // Unpack more of the subject if required
  unpack_extendRegionEnd(maxRowDropoff - subject + seedSubjectOffset,
                         unpackRegion);

  // If the unpacked subject has moved in memory
  if (subject != unpackRegion->unpackedSubject + seedSubjectOffset) {
    // Update all pointers to point at the new subject
    newSubject = unpackRegion->unpackedSubject + seedSubjectOffset;
    subjectPosition = (subjectPosition - subject) + newSubject;
    rowDropoff = (rowDropoff - subject) + newSubject;
    columnDropoff = (columnDropoff - subject) + newSubject;
    bestSubjectPosition = (bestSubjectPosition - subject) + newSubject;
    subject = newSubject;
    subjectEnd = subject + subjectLength;
  }

  // -----FIRST ROW-----

  // Using first column of the query matrix
  matrixColumn = (*queryPosition);

  // -----FIRST CELL-----
  // Set M value for top-left cell
  match = matrixColumn[*subjectPosition];
  *matchRow = match;
  // Set DUMMY Ix and Iy values, which should never be used
  *insertSrow = constants_gappedExtensionDummyValue;
  *insertQrow = constants_gappedExtensionDummyValue;
  // M came from M
  *tracebackColumn = 0;

  // If this is the best-yet scoring cell
  if (match > bestScore) {
    // Update best start cell data
    bestScore = match;
    dropoffThreshold = bestScore - dropoff;
    bestQueryPosition = queryPosition;
    bestSubjectPosition = subjectPosition;
  }

  // Record match and insertS for this about-to-be-previous cell
  previousMatch = match;
  previousInsertS = *insertSrow;

  subjectDistance = 0;
  subjectPosition++;
  matchRow++;
  insertQrow++;
  insertSrow++;
  tracebackColumn++;

  // ----- REMAINING CELLS -----
  // For each remaining columns in the top row, scanning from left-to-right
  while (subjectPosition < subjectEnd) {
    // Set value for M
    match = matrixColumn[*subjectPosition] - parameters_openGap -
            subjectDistance * parameters_extendGap;
    *matchRow = match;

    // Set value for Ix
    if (previousInsertS - parameters_extendGap >
        previousMatch - parameters_openGap) {
      *insertSrow = previousInsertS - parameters_extendGap;
      // M came from Ix and Ix came from Ix
      *tracebackColumn = 5;
    } else {
      *insertSrow = previousMatch - parameters_openGap;
      // M came from Ix and Ix came from M
      *tracebackColumn = 1;
    }

    // Set DUMMY Iy value, which should never be used
    *insertQrow = constants_gappedExtensionDummyValue;

    // If this is the best-yet scoring cell
    if (match > bestScore) {
      // Update best start cell data
      bestScore = match;
      dropoffThreshold = bestScore - dropoff;
      bestQueryPosition = queryPosition;
      bestSubjectPosition = subjectPosition;
    }

    // If score at current cell is below dropoff
    if (dropoffThreshold > match && dropoffThreshold > *insertSrow) {
      // Record dropoff position
      rowDropoff = subjectPosition;
      // And stop processing row
      break;
    }

    // Record match and insertS for this about-to-be-previous cell
    previousMatch = match;
    previousInsertS = *insertSrow;

    subjectPosition++;
    matchRow++;
    insertQrow++;
    insertSrow++;
    tracebackColumn++;
    subjectDistance++;
  }

  //    if (dloc==88197331)
  //    print2(gappedExtension_matchRow_multi[queryNum], subject, rowDropoff,
  // columnDropoff);

  queryDistance = 0;
  queryPosition++;
  tracebackRow++;

  // -----REMAINING ROWS-----
  while (queryPosition < queryEnd && rowDropoff > columnDropoff) {
    subjectPosition = columnDropoff + 1;

    // Calculate maximum row dropoff after this row is processed
    maxRowDropoff =
        rowDropoff + ((PSSMatrix.highestValue / parameters_extendGap) + 2);
    if (maxRowDropoff > subjectEnd)
      maxRowDropoff = subjectEnd;

    // If not enough space in traceback data for this row, realloc
    if (maxRowDropoff - subjectPosition >=
        gappedExtension_tracebackAlloc_multi[queryNum] -
            (tracebackDataPtr -
             gappedExtension_tracebackData_multi[queryNum])) {
      newTracebackData = (unsigned char *)global_malloc(
          sizeof(char) * gappedExtension_tracebackAlloc_multi[queryNum] * 2);

      // Move existing data to end of new block
      memcpy(newTracebackData, gappedExtension_tracebackData_multi[queryNum],
             gappedExtension_tracebackAlloc_multi[queryNum]);

      // Update traceback data pointer
      tracebackDataPtr =
          newTracebackData +
          (tracebackDataPtr - gappedExtension_tracebackData_multi[queryNum]);

      // Update existing rows to point to new data block
      while (tracebackRow > gappedExtension_traceback_multi[queryNum]) {
        tracebackRow--;
        *tracebackRow =
            newTracebackData +
            (*tracebackRow - gappedExtension_tracebackData_multi[queryNum]);
      }

      // Data block is now double the size
      gappedExtension_tracebackAlloc_multi[queryNum] *= 2;

      free(gappedExtension_tracebackData_multi[queryNum]);
      gappedExtension_tracebackData_multi[queryNum] = newTracebackData;
    }

    // Initialize traceback pointer for this row
    tracebackRow = gappedExtension_traceback_multi[queryNum] +
                   (queryPosition - PSSMatrix.matrix);
    *tracebackRow = tracebackDataPtr - (subjectPosition - subject);
    tracebackColumn = tracebackDataPtr;
    tracebackDataPtr += (maxRowDropoff - subjectPosition + 1);

    // Unpack more of the subject if required
    unpack_extendRegionEnd(maxRowDropoff - subject + seedSubjectOffset,
                           unpackRegion);

    // If the unpacked subject has moved in memory
    if (subject != unpackRegion->unpackedSubject + seedSubjectOffset) {
      // Update all pointers to point at the new subject
      newSubject = unpackRegion->unpackedSubject + seedSubjectOffset;
      subjectPosition = (subjectPosition - subject) + newSubject;
      rowDropoff = (rowDropoff - subject) + newSubject;
      columnDropoff = (columnDropoff - subject) + newSubject;
      bestSubjectPosition = (bestSubjectPosition - subject) + newSubject;
      subject = newSubject;
      subjectEnd = subject + subjectLength;
    }

    // Reset rows
    rowOffset = (subjectPosition - subject);
    matchRow = gappedExtension_matchRow_multi[queryNum] +
               (rowOffset % gappedExtension_rowSizes_multi[queryNum]);
    insertQrow = gappedExtension_insertQrow_multi[queryNum] +
                 (rowOffset % gappedExtension_rowSizes_multi[queryNum]);
    insertSrow = gappedExtension_insertSrow_multi[queryNum] +
                 (rowOffset % gappedExtension_rowSizes_multi[queryNum]);

    // Using next column of the query matrix
    matrixColumn = (*queryPosition);

    // -----FAR LEFT CELL-----
    // Record some old values
    previousOldMatch = *matchRow;
    previousOldInsertQ = *insertQrow;
    previousOldInsertS = *insertSrow;

    // Set Iy value
    if (*insertQrow - parameters_extendGap > *matchRow - parameters_openGap) {
      *insertQrow = *insertQrow - parameters_extendGap;
      // Iy is derived from Iy, M is derived from Iy
      *tracebackColumn = 34;
    } else {
      *insertQrow = *matchRow - parameters_openGap;
      // Iy is derived from M, M is derived from Iy
      *tracebackColumn = 2;
    }

    // Set DUMMY values for M and Iy, which should never be used
    match = *matchRow = constants_gappedExtensionDummyValue;
    *insertSrow = constants_gappedExtensionDummyValue;

    // If score at current cell is below dropoff
    if (dropoffThreshold > *insertQrow) {
      // Record dropoff position
      columnDropoff = subjectPosition;
      leftOfDropoff = 1;
    } else {
      // We are left of the column dropoff for this row
      leftOfDropoff = 0;
    }

    // Record match and insertS for this about-to-be-previous cell
    previousMatch = match;
    previousInsertS = *insertSrow;

    subjectPosition++;
    matchRow++;
    insertQrow++;
    insertSrow++;
    tracebackColumn++;

    // Check for scoring rows wrap-around
    if (matchRow >= endMatchRow) {
      matchRow -= gappedExtension_rowSizes_multi[queryNum];
      insertQrow -= gappedExtension_rowSizes_multi[queryNum];
      insertSrow -= gappedExtension_rowSizes_multi[queryNum];
    }

    // -----CELLS LEFT OF ROW DROPOFF-----
    while (subjectPosition <= rowDropoff) {
      // Remember old M value (for cell below this one)
      oldMatch = *matchRow;

      // Calculate new M value
      if (previousOldMatch >= previousOldInsertQ) {
        if (previousOldMatch >= previousOldInsertS) {
          match = matrixColumn[*subjectPosition] + previousOldMatch;
          // M is derived from M
          *tracebackColumn = 0;
        } else {
          match = matrixColumn[*subjectPosition] + previousOldInsertS;
          // M is derived from Ix
          *tracebackColumn = 1;
        }
      } else {
        if (previousOldInsertQ >= previousOldInsertS) {
          match = matrixColumn[*subjectPosition] + previousOldInsertQ;
          // M is derived from Iy
          *tracebackColumn = 2;
        } else {
          match = matrixColumn[*subjectPosition] + previousOldInsertS;
          // M is derived from Ix
          *tracebackColumn = 1;
        }
      }

      *matchRow = match;

      // Record some old values
      previousOldMatch = oldMatch;
      previousOldInsertQ = *insertQrow;
      previousOldInsertS = *insertSrow;

      // Set new Iy value
      if (oldMatch - parameters_openGap >= *insertQrow - parameters_extendGap) {
        *insertQrow = oldMatch - parameters_openGap;
        // Iy is derived from M
        // No change to traceback
      } else {
        *insertQrow = *insertQrow - parameters_extendGap;
        // Iy is derived from Iy
        *tracebackColumn |= 32;
      }
      // Calculate new Ix
      if (previousMatch - parameters_openGap >=
          previousInsertS - parameters_extendGap) {
        *insertSrow = previousMatch - parameters_openGap;
        // Ix is derived from M
        // No change to traceback
      } else {
        *insertSrow = previousInsertS - parameters_extendGap;
        // Ix is derived from Ix
        *tracebackColumn |= 4;
      }

      // If this is the best-yet scoring cell
      if (match > bestScore) {
        // Update best start cell data
        bestScore = match;
        dropoffThreshold = bestScore - dropoff;
        bestQueryPosition = queryPosition;
        bestSubjectPosition = subjectPosition;
      }

      // If score at current cell (and cells to its left) are below dropoff
      if (leftOfDropoff) {
        if (dropoffThreshold > match && dropoffThreshold > *insertSrow &&
            dropoffThreshold > *insertQrow) {
          // Record dropoff position
          columnDropoff = subjectPosition;
        } else {
          // We are left of the column dropoff for this row
          leftOfDropoff = 0;
        }
      }

      // Record match and insertS for this about-to-be-previous cell
      previousMatch = match;
      previousInsertS = *insertSrow;
      previousInsertQ = *insertQrow;

      subjectPosition++;
      matchRow++;
      insertQrow++;
      insertSrow++;
      tracebackColumn++;

      // Check for scoring rows wrap-around
      if (matchRow >= endMatchRow) {
        matchRow -= gappedExtension_rowSizes_multi[queryNum];
        insertQrow -= gappedExtension_rowSizes_multi[queryNum];
        insertSrow -= gappedExtension_rowSizes_multi[queryNum];
      }
    }

    // -----CELLS RIGHT OF ROW DROPOFF -----
    if (!(dropoffThreshold > previousMatch &&
          dropoffThreshold > previousInsertS &&
          dropoffThreshold > previousInsertQ)) {
      while (subjectPosition < subjectEnd) {
        *insertSrow = previousInsertS - parameters_extendGap;
        // Ix came from Ix
        *tracebackColumn = 4;

        // Set DUMMY values for M and Ix, which should never be used
        *matchRow = constants_gappedExtensionDummyValue;
        *insertQrow = constants_gappedExtensionDummyValue;

        // If score at current cell is below dropoff
        if (dropoffThreshold > *insertSrow) {
          // And stop processing row
          subjectPosition++;
          break;
        }

        // Record insertS for this about-to-be-previous cell
        previousInsertS = *insertSrow;

        subjectPosition++;
        matchRow++;
        insertQrow++;
        insertSrow++;
        tracebackColumn++;

        // Check for scoring rows wrap-around
        if (matchRow >= endMatchRow) {
          matchRow -= gappedExtension_rowSizes_multi[queryNum];
          insertQrow -= gappedExtension_rowSizes_multi[queryNum];
          insertSrow -= gappedExtension_rowSizes_multi[queryNum];
        }

        subjectDistance++;
      }
    }

    // Record dropoff position
    rowDropoff = subjectPosition - 1;

    queryDistance++;
    queryPosition++;
    tracebackRow++;

    //        if (dloc==88197331)
    //			gappedExtension_printAfterRow(gappedExtension_matchRow_multi[queryNum],
    // subject, rowDropoff,
    //                                          columnDropoff);
  }

  dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
  dpResults.best.subjectOffset = bestSubjectPosition - subject;
  dpResults.bestScore = bestScore;
  dpResults.traceback = gappedExtension_traceback_multi[queryNum];

  return dpResults;
}

// Given a gapped extension with a nominal score, calculate the normalized score
// and E-Value
void gappedExtension_score_multi(struct gappedExtension *gappedExtension,
                                 int queryNum) {
  gappedExtension->normalizedScore =
      statistics_gappedNominal2normalized(gappedExtension->nominalScore);

  gappedExtension->eValue = statistics_gappedCalculateEvalue_multi(
      gappedExtension->normalizedScore, queryNum);
}

void gappedExtension_score_ncbi_multi(struct gappedExtension *gappedExtension,
                                 int queryNum) {
  gappedExtension->normalizedScore =
      statistics_gappedNominal2normalized_ncbi(gappedExtension->nominalScore);

  gappedExtension->eValue = gappedExtension->eValue;
}

void gappedExtension_free_multi(int queryNum) {
  // Free memory used by traceback array
  free(gappedExtension_traceback_multi[queryNum]);
  free(gappedExtension_tracebackData_multi[queryNum]);

  //    printf("gappedExtension_tracebackAlloc_multi[queryNum]=%d\n",
  // gappedExtension_tracebackAlloc_multi[queryNum]);

  // Free memory used to store row scores
  free(gappedExtension_matchRow_multi[queryNum]);
  free(gappedExtension_insertQrow_multi[queryNum]);
  free(gappedExtension_insertSrow_multi[queryNum]);

  gappedExtension_rowSizes_multi[queryNum] = 0;
  gappedExtension_numRows_multi[queryNum] = 0;
  gappedExtension_tracebackAlloc_multi[queryNum] = 0;
}
