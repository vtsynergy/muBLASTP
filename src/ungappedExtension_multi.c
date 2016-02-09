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

int4 ungappedExtension_minus3reward;
int4 ungappedExtension_tableMatchesReward;
struct memBlocks *ungappedExtension_extensions_multi[BATCH_SIZE];
struct memBlocks ***ungappedExtension_extensions_multi2;

// Initialize the creation of ungapped extensions
void ungappedExtension_initialize_multi() {
  // ungappedExtension_extensions_multi = (struct memBlocks
  // **)malloc(sizeof(struct memBlocks *) * blast_numQuery);
  ungappedExtension_extensions_multi2 = (struct memBlocks ***)malloc(
      sizeof(struct memBlocks **) * blast_numQuery);

  int ii;
  for (ii = 0; ii < blast_numQuery; ii++) {
    ungappedExtension_extensions_multi2[ii] = (struct memBlocks **)malloc(
        sizeof(struct memBlocks *) * blast_numBlocks);
    int jj;
    for (jj = 0; jj < blast_numBlocks; jj++) {
      ungappedExtension_extensions_multi2[ii][jj] = memBlocks_initialize(
          sizeof(struct ungappedExtension),
          constants_initialAllocUngappedExtensions);
    }

    ungappedExtension_extensions_multi[ii] =
        memBlocks_initialize(sizeof(struct ungappedExtension),
                             constants_initialAllocUngappedExtensions);
  }
  ungappedExtension_minus3reward = parameters_matchScore * -3;
  ungappedExtension_tableMatchesReward =
      parameters_matchScore * parameters_wordTableLetters;
}

// Initialize the creation of ungapped extensions
void ungappedExtension_initialize_multi2() {
  // ungappedExtension_extensions_multi = (struct memBlocks
  // **)malloc(sizeof(struct memBlocks *) * blast_numQuery);
#if 0
  ungappedExtension_extensions_multi2 = (struct memBlocks ***)malloc(
      sizeof(struct memBlocks **) * blast_numQuery);

  int ii;
  for (ii = 0; ii < blast_numQuery; ii++) {
    ungappedExtension_extensions_multi2[ii] = (struct memBlocks **)malloc(
        sizeof(struct memBlocks *) * blast_numBlocks);
    int jj;
    for (jj = 0; jj < blast_numBlocks; jj++) {
      ungappedExtension_extensions_multi2[ii][jj] = memBlocks_initialize(
          sizeof(struct ungappedExtension),
          constants_initialAllocUngappedExtensions);
    }

    ungappedExtension_extensions_multi[ii] =
        memBlocks_initialize(sizeof(struct ungappedExtension),
                             constants_initialAllocUngappedExtensions);
  }
#endif
  ungappedExtension_minus3reward = parameters_matchScore * -3;
  ungappedExtension_tableMatchesReward =
      parameters_matchScore * parameters_wordTableLetters;
}

// Perform an ungapped extension between points queryStart,subjectStart and
// queryEnd,subjectEnd
// and extend in each direction until score drops below best score yet minus a
// dropoff parameter
struct ungappedExtension *ungappedExtension_extend_multi(
    int2 **queryHit, unsigned char *subjectHit, unsigned char *lastHit,
    struct PSSMatrix PSSMatrix, unsigned char *subject, int queryNum,
    unsigned char **ungappedExtension_subjectEndReached) {
  int2 **queryPosition;
  unsigned char *subjectPosition, *subjectStart, *subjectEnd;
  int4 changeSinceBest = 0;
  int4 dropoff, originalDropoff;

  originalDropoff = dropoff =
      -statistics_ungappedNominalDropoff_multi[queryNum];
  int ungappedExtension_bestScore = 0;

  // Start at queryEnd,subjectEnd (right/last hit position)
  queryPosition = queryHit;
  subjectPosition = subjectStart = subjectHit;

  // Extend the start of the hit backwards until dropoff
  while (changeSinceBest > dropoff) {
    changeSinceBest += (*queryPosition)[*subjectPosition];

    // If we have got a positive score
    if (changeSinceBest > 0) {
      // Keep updating best score and resetting change-since-best
      // whilst we are reading positive scores
      do {
        ungappedExtension_bestScore += changeSinceBest;
        queryPosition--;
        subjectPosition--;
        changeSinceBest = (*queryPosition)[*subjectPosition];
      } while (changeSinceBest > 0);

      subjectStart = subjectPosition;
    }
    queryPosition--;
    subjectPosition--;
  }

  // Correct for extra decrement
  subjectStart++;

  // If best start point is right of previous hit which helped trigger this
  // extension
  // then stop now
  if (subjectStart > lastHit) {
    *ungappedExtension_subjectEndReached = subjectHit;
    return NULL;
  }

  // Starting at right/last hit position again
  queryPosition = queryHit + 1;
  subjectEnd = subjectHit;
  subjectPosition = subjectHit + 1;
  changeSinceBest = 0;

  // May need to alter dropoff so we also dropoff if below zero
  if (-ungappedExtension_bestScore > originalDropoff) {
    dropoff = -ungappedExtension_bestScore;
  }

  // Extend end of alignment until dropoff
  while (changeSinceBest > dropoff) {
    changeSinceBest += (*queryPosition)[*subjectPosition];

    // If we have got a positive score
    if (changeSinceBest > 0) {
      // Keep updating best score and resetting change-since-best
      // whilst we are reading positive scores
      do {
        ungappedExtension_bestScore += changeSinceBest;
        queryPosition++;
        subjectPosition++;
        changeSinceBest = (*queryPosition)[*subjectPosition];
      } while (changeSinceBest > 0);

      subjectEnd = subjectPosition;

      // Check need for change in dropoff
      if ((dropoff = -ungappedExtension_bestScore) < originalDropoff) {
        dropoff = originalDropoff;
      }
    }
    queryPosition++;
    subjectPosition++;
  }

  // Correct for extra increment
  subjectEnd--;
  *ungappedExtension_subjectEndReached = subjectEnd;

  // If extension scored above trigger for gapping, create object and return it
  if (ungappedExtension_bestScore >=
      blast_ungappedNominalTrigger_multi[queryNum]) {
    int4 diagonal;
    struct ungappedExtension newUngappedExtension_t;
    struct ungappedExtension *newUngappedExtension = &newUngappedExtension_t;
    newUngappedExtension =
        memBlocks_newEntry(ungappedExtension_extensions_multi[queryNum]);

    // Calculate diagonal
    diagonal = (subjectHit - subject) - (queryHit - PSSMatrix.matrix);

    // Determine offsets from pointers
    newUngappedExtension->start.subjectOffset = subjectStart - subject;
    newUngappedExtension->end.subjectOffset = subjectEnd - subject;
    newUngappedExtension->start.queryOffset =
        newUngappedExtension->start.subjectOffset - diagonal;
    newUngappedExtension->end.queryOffset =
        newUngappedExtension->end.subjectOffset - diagonal;

    // Find the seed point
    newUngappedExtension->seed = ungappedExtension_findProteinSeed(
        newUngappedExtension, PSSMatrix, subject);

    // Initialize next to null
    newUngappedExtension->next = NULL;
    newUngappedExtension->nominalScore = ungappedExtension_bestScore;
    newUngappedExtension->status = ungappedExtension_UNGAPPED;

    return newUngappedExtension;
  } else {
    return NULL;
  }
}

// Perform an ungapped extension between points queryStart,subjectStart and
// queryEnd,subjectEnd
// and extend in each direction until score drops below best score yet minus a
// dropoff parameter
struct ungappedExtension *ungappedExtension_extend_multi2(
    int2 **queryHit, unsigned char *subjectHit, unsigned char *lastHit,
    struct PSSMatrix PSSMatrix, unsigned char *subject, int queryNum,
    int blockNum, unsigned char **ungappedExtension_subjectEndReached) {
  int2 **queryPosition;
  unsigned char *subjectPosition, *subjectStart, *subjectEnd;
  int4 changeSinceBest = 0;
  int4 dropoff, originalDropoff;

  originalDropoff = dropoff =
      -statistics_ungappedNominalDropoff_multi[queryNum];
  int ungappedExtension_bestScore = 0;

  // Start at queryEnd,subjectEnd (right/last hit position)
  queryPosition = queryHit;
  subjectPosition = subjectStart = subjectHit;

  // Extend the start of the hit backwards until dropoff
  while (changeSinceBest > dropoff) {
    changeSinceBest += (*queryPosition)[*subjectPosition];

    // If we have got a positive score
    if (changeSinceBest > 0) {
      // Keep updating best score and resetting change-since-best
      // whilst we are reading positive scores
      do {
        ungappedExtension_bestScore += changeSinceBest;
        queryPosition--;
        subjectPosition--;
        changeSinceBest = (*queryPosition)[*subjectPosition];
      } while (changeSinceBest > 0);

      subjectStart = subjectPosition;
    }
    queryPosition--;
    subjectPosition--;
  }

  // Correct for extra decrement
  subjectStart++;

  // If best start point is right of previous hit which helped trigger this
  // extension
  // then stop now
  if (subjectStart > lastHit) {
    *ungappedExtension_subjectEndReached = subjectHit;
    return NULL;
  }

  // Starting at right/last hit position again
  queryPosition = queryHit + 1;
  subjectEnd = subjectHit;
  subjectPosition = subjectHit + 1;
  changeSinceBest = 0;

  // May need to alter dropoff so we also dropoff if below zero
  if (-ungappedExtension_bestScore > originalDropoff) {
    dropoff = -ungappedExtension_bestScore;
  }

  // Extend end of alignment until dropoff
  while (changeSinceBest > dropoff) {
    changeSinceBest += (*queryPosition)[*subjectPosition];

    // If we have got a positive score
    if (changeSinceBest > 0) {
      // Keep updating best score and resetting change-since-best
      // whilst we are reading positive scores
      do {
        ungappedExtension_bestScore += changeSinceBest;
        queryPosition++;
        subjectPosition++;
        changeSinceBest = (*queryPosition)[*subjectPosition];
      } while (changeSinceBest > 0);

      subjectEnd = subjectPosition;

      // Check need for change in dropoff
      if ((dropoff = -ungappedExtension_bestScore) < originalDropoff) {
        dropoff = originalDropoff;
      }
    }
    queryPosition++;
    subjectPosition++;
  }

  // Correct for extra increment
  subjectEnd--;
  *ungappedExtension_subjectEndReached = subjectEnd;

  // If extension scored above trigger for gapping, create object and return it
  if (ungappedExtension_bestScore >=
      blast_ungappedNominalTrigger_multi[queryNum]) {
    int4 diagonal;
    struct ungappedExtension *newUngappedExtension;
    newUngappedExtension = memBlocks_newEntry(
        ungappedExtension_extensions_multi2[queryNum][blockNum]);

    // Calculate diagonal
    diagonal = (subjectHit - subject) - (queryHit - PSSMatrix.matrix);

    // Determine offsets from pointers
    newUngappedExtension->start.subjectOffset = subjectStart - subject;
    newUngappedExtension->end.subjectOffset = subjectEnd - subject;
    newUngappedExtension->start.queryOffset =
        newUngappedExtension->start.subjectOffset - diagonal;
    newUngappedExtension->end.queryOffset =
        newUngappedExtension->end.subjectOffset - diagonal;

    // Find the seed point
    newUngappedExtension->seed = ungappedExtension_findProteinSeed(
        newUngappedExtension, PSSMatrix, subject);

    // Initialize next to null
    newUngappedExtension->next = NULL;
    newUngappedExtension->nominalScore = ungappedExtension_bestScore;
    newUngappedExtension->status = ungappedExtension_UNGAPPED;

    return newUngappedExtension;
  } else {
    return NULL;
  }
}

// Perform an ungapped extension between points queryStart,subjectStart and
// queryEnd,subjectEnd
// and extend in each direction until score drops below best score yet minus a
// dropoff parameter
int ungappedExtension_extend_multi3(
    struct ungappedExtension *newUngappedExtension,
    int2 **queryHit, unsigned char *subjectHit, unsigned char *lastHit,
    struct PSSMatrix PSSMatrix, unsigned char *subject, int queryNum,
    int blockNum, unsigned char **ungappedExtension_subjectEndReached) {
  int2 **queryPosition;
  unsigned char *subjectPosition, *subjectStart, *subjectEnd;
  int4 changeSinceBest = 0;
  int4 dropoff, originalDropoff;

  originalDropoff = dropoff =
      -statistics_ungappedNominalDropoff_multi[queryNum];
  int ungappedExtension_bestScore = 0;

  // Start at queryEnd,subjectEnd (right/last hit position)
  queryPosition = queryHit;
  subjectPosition = subjectStart = subjectHit;

  // Extend the start of the hit backwards until dropoff
  while (changeSinceBest > dropoff) {
    changeSinceBest += (*queryPosition)[*subjectPosition];

    // If we have got a positive score
    if (changeSinceBest > 0) {
      // Keep updating best score and resetting change-since-best
      // whilst we are reading positive scores
      do {
        ungappedExtension_bestScore += changeSinceBest;
        queryPosition--;
        subjectPosition--;
        changeSinceBest = (*queryPosition)[*subjectPosition];
      } while (changeSinceBest > 0);

      subjectStart = subjectPosition;
    }
    queryPosition--;
    subjectPosition--;
  }

  // Correct for extra decrement
  subjectStart++;

  // If best start point is right of previous hit which helped trigger this
  // extension
  // then stop now
  if (subjectStart > lastHit) {
    *ungappedExtension_subjectEndReached = subjectHit;
    return 0;
  }

  // Starting at right/last hit position again
  queryPosition = queryHit + 1;
  subjectEnd = subjectHit;
  subjectPosition = subjectHit + 1;
  changeSinceBest = 0;

  // May need to alter dropoff so we also dropoff if below zero
  if (-ungappedExtension_bestScore > originalDropoff) {
    dropoff = -ungappedExtension_bestScore;
  }

  // Extend end of alignment until dropoff
  while (changeSinceBest > dropoff) {
    changeSinceBest += (*queryPosition)[*subjectPosition];

    // If we have got a positive score
    if (changeSinceBest > 0) {
      // Keep updating best score and resetting change-since-best
      // whilst we are reading positive scores
      do {
        ungappedExtension_bestScore += changeSinceBest;
        queryPosition++;
        subjectPosition++;
        changeSinceBest = (*queryPosition)[*subjectPosition];
      } while (changeSinceBest > 0);

      subjectEnd = subjectPosition;

      // Check need for change in dropoff
      if ((dropoff = -ungappedExtension_bestScore) < originalDropoff) {
        dropoff = originalDropoff;
      }
    }
    queryPosition++;
    subjectPosition++;
  }

  // Correct for extra increment
  subjectEnd--;
  *ungappedExtension_subjectEndReached = subjectEnd;

  // If extension scored above trigger for gapping, create object and return it
  if (ungappedExtension_bestScore >=
      blast_ungappedNominalTrigger_multi[queryNum]) {
    int4 diagonal;

    // Calculate diagonal
    diagonal = (subjectHit - subject) - (queryHit - PSSMatrix.matrix);

    // Determine offsets from pointers
    newUngappedExtension->start.subjectOffset = subjectStart - subject;
    newUngappedExtension->end.subjectOffset = subjectEnd - subject;
    newUngappedExtension->start.queryOffset =
        newUngappedExtension->start.subjectOffset - diagonal;
    newUngappedExtension->end.queryOffset =
        newUngappedExtension->end.subjectOffset - diagonal;

    // Find the seed point
    newUngappedExtension->seed = ungappedExtension_findProteinSeed(
        newUngappedExtension, PSSMatrix, subject);

    // Initialize next to null
    newUngappedExtension->next = NULL;
    newUngappedExtension->nominalScore = ungappedExtension_bestScore;
    newUngappedExtension->status = ungappedExtension_UNGAPPED;

    return 1;
  } else {
    return 0;
  }
}

static int4 s_BlastPSSMExtendLeft(int2 **matrix, char *subject, int4 s_off,
                                  int4 q_off, int4 dropoff, int4 *length,
                                  int4 maxscore) {
  int4 i, n, best_i;
  int4 score = maxscore;
  char *s;

  n = min(s_off, q_off);
  best_i = n + 1;
  s = subject + s_off - n;

  for (i = n; i >= 0; i--) {
    score += matrix[q_off - n + i][s[i]];

    if (score > maxscore) {
      maxscore = score;
      best_i = i;
    }
    /* The comparison below is really >= and is different than the old
       code (e.g., blast.c:BlastWordExtend_prelim). In the old code the
       loop continued as long as sum > X (X being negative).  The loop
       control here is different and we *break out* when the if statement
       below is true. */
    if ((maxscore - score) >= dropoff)
      break;
  }

  *length = n - best_i + 1;
  return maxscore;
}

static int4 s_BlastPSSMExtendRight(int2 **matrix, char *subject,
                                   int4 query_size, int4 s_off, int4 q_off,
                                   int4 dropoff, int4 *length, int4 maxscore,
                                   int4 *s_last_off, int subjectLength) {
  int4 i, n, best_i = -1;
  int4 score = maxscore;
  unsigned char *s;

  n = min(subjectLength - s_off, query_size - q_off);
  s = subject + s_off;

  for (i = 0; i < n; i++) {
    score += matrix[q_off + i][s[i]];

    if (score > maxscore) {
      maxscore = score;
      best_i = i;
    }

    /* The comparison below is really >= and is different than the old
       code (e.g., blast.c:BlastWordExtend_prelim). In the old code the
       loop continued as long as sum > X (X being negative).  The loop
       control here is different and we *break out* when the if statement
       below is true. */
    if (score <= 0 || (maxscore - score) >= dropoff)
      break;
  }

  *length = best_i + 1;
  *s_last_off = s_off + i;
  return maxscore;
}

int4 s_BlastAaExtendLeft(int2 **matrix, char *subject, char *query, int4 s_off,
                         int4 q_off, int4 dropoff, int4 *length,
                         int4 maxscore) {
  int4 i, n, best_i;
  int4 score = maxscore;

  unsigned char *s, *q;

  n = min(s_off, q_off);
  best_i = n + 1;

  s = subject + s_off - n;
  q = query + q_off - n;

  for (i = n; i >= 0; i--) {
    score += matrix[q[i]][s[i]];

    if (score > maxscore) {
      maxscore = score;
      best_i = i;
    }
    /* The comparison below is really >= and is different than the old
       code (e.g., blast.c:BlastWordExtend_prelim). In the old code the
       loop continued as long as sum > X (X being negative).  The loop
       control here is different and we *break out* when the if statement
       below is true. */
    if ((maxscore - score) >= dropoff)
      break;
  }

  *length = n - best_i + 1;
  return maxscore;
}

int4 s_BlastAaExtendRight(int2 **matrix, char *subject, char *query, int4 s_off,
                          int4 q_off, int4 dropoff, int4 *length, int4 maxscore,
                          int4 *s_last_off, int subjectLength,
                          int queryLength) {
  int4 i, n, best_i = -1;
  int4 score = maxscore;

  unsigned char *s, *q;
  n = min(subjectLength - s_off, queryLength - q_off);

  s = subject + s_off;
  q = query + q_off;

  for (i = 0; i < n; i++) {
    score += matrix[q[i]][s[i]];

    if (score > maxscore) {
      maxscore = score;
      best_i = i;
    }

    /* The comparison below is really >= and is different than the old
       code (e.g., blast.c:BlastWordExtend_prelim). In the old code the
       loop continued as long as sum > X (X being negative).  The loop
       control here is different and we *break out* when the if statement
       below is true. */
    if (score <= 0 || (maxscore - score) >= dropoff)
      break;
  }

  *length = best_i + 1;
  *s_last_off = s_off + i;
  return maxscore;
}

int4 s_BlastAaExtendTwoHit_SM(int2 **matrix, char *query, char *subject,
                              int4 s_left_off, int4 s_right_off,
                              int4 q_right_off, int4 dropoff, int4 *hsp_q,
                              int4 *hsp_s, int4 *hsp_len, int4 word_size,
                              int4 *s_last_off, int subjectLength,
                              int queryLength, int4 sequenceCount, char *rightExtend) {
  int4 left_d = 0, right_d = 0;         /* left and right displacements */
  int4 left_score = 0, right_score = 0; /* left and right scores */
  int4 i, score = 0;
  char *s = subject;
  char *q = query;

  /* find one beyond the position (up to word_size-1 letters to the right)
     that gives the largest starting score */
  /* Use "one beyond" to make the numbering consistent with how it's done
     for BlastAaExtendOneHit and the "Extend" functions called here. */
  for (i = 0; i < word_size; i++) {
    score += matrix[q[q_right_off + i]][s[s_right_off + i]];
    if (score > left_score) {
      left_score = score;
      right_d = i + 1; /* Position is one beyond the end of the
                          word. */
    }
  }

  q_right_off += right_d;
  s_right_off += right_d;

  right_d = 0;
  *s_last_off = s_right_off;

  *rightExtend = 0;

  /* first, try to extend left, from the second hit to the first hit. */
  left_score = s_BlastAaExtendLeft(matrix, subject, query, s_right_off - 1,
                                   q_right_off - 1, dropoff, &left_d, 0);

  /* Extend to the right only if left extension reached the first hit. */
  if (left_d >= (s_right_off - s_left_off)) {

      *rightExtend = 1;
    right_score = s_BlastAaExtendRight(
        matrix, subject, query, s_right_off, q_right_off, dropoff, &right_d,
        left_score, s_last_off, subjectLength, queryLength);
  }

  *hsp_q = q_right_off - left_d;
  *hsp_s = s_right_off - left_d;
  *hsp_len = left_d + right_d;

  return max(left_score, right_score);
}

int4 s_BlastAaExtendTwoHit_PSSM(int2 **matrix, char *subject, int4 s_left_off,
                                int4 s_right_off, int4 q_right_off,
                                int4 dropoff, int4 *hsp_q, int4 *hsp_s,
                                int4 *hsp_len, int4 word_size, int4 *s_last_off,
                                int subjectLength, int queryLength,
                                int4 sequenceCount, char *rightExtend) {
  int4 left_d = 0, right_d = 0;         /* left and right displacements */
  int4 left_score = 0, right_score = 0; /* left and right scores */
  int4 i, score = 0;
  char *s = subject;

  /* find one beyond the position (up to word_size-1 letters to the right)
     that gives the largest starting score */
  /* Use "one beyond" to make the numbering consistent with how it's done
     for BlastAaExtendOneHit and the "Extend" functions called here. */
  for (i = 0; i < word_size; i++) {
    score += matrix[q_right_off + i][s[s_right_off + i]];
    if (score > left_score) {
      left_score = score;
      right_d = i + 1; /* Position is one beyond the end of the
                          word. */
    }
  }

  q_right_off += right_d;
  s_right_off += right_d;

  right_d = 0;
  *s_last_off = s_right_off;

  *rightExtend = 0;

  /* first, try to extend left, from the second hit to the first hit. */
  left_score = s_BlastPSSMExtendLeft(matrix, subject, s_right_off - 1,
                                     q_right_off - 1, dropoff, &left_d, 0);

  /* Extend to the right only if left extension reached the first hit. */
  if (left_d >= (s_right_off - s_left_off)) {

      *rightExtend = 1;
    right_score = s_BlastPSSMExtendRight(
        matrix, subject, queryLength, s_right_off, q_right_off, dropoff,
        &right_d, left_score, s_last_off, subjectLength);
  }

  *hsp_q = q_right_off - left_d;
  *hsp_s = s_right_off - left_d;
  *hsp_len = left_d + right_d;

  return max(left_score, right_score);
}

struct ungappedExtension *ungappedExtension_extend_ncbi_multi2(
    struct PSSMatrix PSSMatrix, struct scoreMatrix scoreMatrix, char *subject,
    int4 lastHitOffset, int4 subjectOffset, int4 queryOffset, int subjectLength,
    int queryLength, int4 sequenceCount,
    unsigned char **ungappedExtension_subjectEndReached, int queryNum,
    struct ungappedExtension *goodExtensionBuf, int *goodExtensionCount,
    char *rightExtend) {
  int4 subjectStart, subjectEnd, queryStart, queryEnd, extendLength,
      subjectEndReached;
  int4 ungappedExtension_bestScore;

#if 0
  ungappedExtension_bestScore = s_BlastAaExtendTwoHit_PSSM(
      PSSMatrix.matrix, subject, lastHitOffset, subjectOffset, queryOffset,
      statistics_ungappedNominalDropoff_multi[queryNum], &queryStart,
      &subjectStart, &extendLength, parameters_wordSize, &subjectEndReached,
      subjectLength, queryLength, sequenceCount, rightExtend);
#else
  ungappedExtension_bestScore = s_BlastAaExtendTwoHit_SM(
      scoreMatrix.matrix, PSSMatrix.queryCodes, subject, lastHitOffset,
      subjectOffset, queryOffset,
      statistics_ungappedNominalDropoff_multi[queryNum], &queryStart,
      &subjectStart, &extendLength, parameters_wordSize, &subjectEndReached,
      subjectLength, queryLength, sequenceCount, rightExtend);
#endif



  *ungappedExtension_subjectEndReached = subject + subjectEndReached;

  //if(sequenceCount == 545626)
  //fprintf(stderr, "hit id: %d diag: %d q: %d s: %d\n", sequenceCount, (queryStart - subjectStart) & 2047, queryOffset, subjectOffset);

  if (ungappedExtension_bestScore >=
      blast_ungappedNominalTrigger_multi[queryNum]) {

      //fprintf(stderr, "ext id: %d diag: %d q: %d - %d s: %d - %d score: %d\n", sequenceCount, (queryStart - subjectStart) & 2047, queryStart, queryStart + extendLength, subjectStart, subjectStart + extendLength, ungappedExtension_bestScore);

    int4 diagonal;
    struct ungappedExtension *newUngappedExtension = goodExtensionBuf + *goodExtensionCount;
    (*goodExtensionCount)++;
    //if(*goodExtensionCount >= MAX_EXTENSIONS_PER_QUERY)
    //{
        //fprintf(stderr, "goodExtensionCount larger than MAX_EXTENSIONS_PER_QUERY\n");
        //exit(1);
    //}
    //newUngappedExtension = memBlocks_newEntry(
    //ungappedExtension_extensions_multi2[queryNum][blockNum]);

    // Determine offsets from pointers
    newUngappedExtension->start.subjectOffset = subjectStart;
    newUngappedExtension->end.subjectOffset = subjectStart + extendLength;
    newUngappedExtension->start.queryOffset = queryStart;
    newUngappedExtension->end.queryOffset = queryStart + extendLength;

    // Find the seed point
    newUngappedExtension->seed = ungappedExtension_findProteinSeed(
        newUngappedExtension, PSSMatrix, subject);
    // Initialize next to null
    newUngappedExtension->next = NULL;
    newUngappedExtension->nominalScore = ungappedExtension_bestScore;
    newUngappedExtension->status = ungappedExtension_UNGAPPED;
    newUngappedExtension->sequenceCount = sequenceCount;

    return newUngappedExtension;
  } else {
    return NULL;
  }
}

int ungappedExtension_extend_ncbi_multi3(
    struct ungappedExtension *newUngappedExtension,
    struct PSSMatrix PSSMatrix, struct scoreMatrix scoreMatrix, char *subject,
    int4 lastHitOffset, int4 subjectOffset, int4 queryOffset, int subjectLength,
    int queryLength, int4 sequenceCount,
    unsigned char **ungappedExtension_subjectEndReached, int queryNum,
    int blockNum, char *rightExtend) {
  int4 subjectStart, subjectEnd, queryStart, queryEnd, extendLength,
      subjectEndReached;
  int4 ungappedExtension_bestScore;

#ifndef SM_UNGAP_EXT
  ungappedExtension_bestScore = s_BlastAaExtendTwoHit_PSSM(
      PSSMatrix.matrix, subject, lastHitOffset, subjectOffset, queryOffset,
      statistics_ungappedNominalDropoff_multi[queryNum], &queryStart,
      &subjectStart, &extendLength, parameters_wordSize, &subjectEndReached,
      subjectLength, queryLength, sequenceCount, rightExtend);
#else
  ungappedExtension_bestScore = s_BlastAaExtendTwoHit_SM(
      scoreMatrix.matrix, PSSMatrix.queryCodes, subject, lastHitOffset,
      subjectOffset, queryOffset,
      statistics_ungappedNominalDropoff_multi[queryNum], &queryStart,
      &subjectStart, &extendLength, parameters_wordSize, &subjectEndReached,
      subjectLength, queryLength, sequenceCount, rightExtend);
#endif



  *ungappedExtension_subjectEndReached = subject + subjectEndReached;

  //if(sequenceCount == 545626)
  //fprintf(stderr, "hit id: %d diag: %d q: %d s: %d\n", sequenceCount, (queryStart - subjectStart) & 2047, queryOffset, subjectOffset);

  if (ungappedExtension_bestScore >=
      blast_ungappedNominalTrigger_multi[queryNum]) {

      //fprintf(stderr, "ext id: %d diag: %d q: %d - %d s: %d - %d score: %d\n", sequenceCount, (queryStart - subjectStart) & 2047, queryStart, queryStart + extendLength, subjectStart, subjectStart + extendLength, ungappedExtension_bestScore);

    int4 diagonal;

    // Determine offsets from pointers
    newUngappedExtension->start.subjectOffset = subjectStart;
    newUngappedExtension->end.subjectOffset = subjectStart + extendLength;
    newUngappedExtension->start.queryOffset = queryStart;
    newUngappedExtension->end.queryOffset = queryStart + extendLength;

    // Find the seed point
    newUngappedExtension->seed = ungappedExtension_findProteinSeed(
        newUngappedExtension, PSSMatrix, subject);
    // Initialize next to null
    newUngappedExtension->next = NULL;
    newUngappedExtension->nominalScore = ungappedExtension_bestScore;
    newUngappedExtension->status = ungappedExtension_UNGAPPED;
    newUngappedExtension->sequenceCount = sequenceCount;

    return 1;
  } else {
    return 0;
  }
}
