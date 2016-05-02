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

unsigned char **hitMatrix_furthest_multi[BATCH_SIZE];
int4 hitMatrix_numDiagonals;

// Initialize the hitMatrix by declaring memory for the maximum number
// of diagonals required by a subject sequence
void hitMatrix_initialize_multi(int4 queryLength, int4 maximumSubjectLength,
                                unsigned char *startAddress, int queryNum) {
  int4 numDiagonals;
  unsigned char **minOffset, **offset, **maxOffset;

  // Use more memory efficient but slower hit matrix for nucleotide
  if (encoding_alphabetType == encoding_nucleotide) {
    // Calculate number of diagonals that will be required during search
    numDiagonals = 1;
    while (numDiagonals < queryLength + parameters_wordSize) {
      numDiagonals <<= 1;
    }

    // Construct mask
    hitMatrix_diagonalMask = numDiagonals - 1;

    // Declare memory for diagonal slots
    hitMatrix_furthest_multi[queryNum] =
        (unsigned char **)global_malloc(sizeof(unsigned char *) * numDiagonals);
    minOffset = hitMatrix_furthest_multi[queryNum];
  }
  // Use less memory efficient but faster hit matrix for protein
  else {
    // Maximum number of diagonals that will be required during search
    numDiagonals = queryLength + maximumSubjectLength - parameters_wordSize + 1;
    minOffset =
        (unsigned char **)global_malloc(sizeof(unsigned char *) * numDiagonals);

    // Advance array pointer to allow offset values ranging from
    // -queryLength to subjectLength - wordSize
    hitMatrix_furthest_multi[queryNum] = minOffset + queryLength;
  }

  hitMatrix_numDiagonals = numDiagonals;

  // Record query length
  hitMatrix_queryLength = queryLength;

  // Start from smallest possible offset value and iterate through to largest
  offset = minOffset;
  maxOffset = minOffset + numDiagonals;

  // For each diagonal, reset furthest to address at start of file
  while (offset < maxOffset) {
    *offset = startAddress;
    offset++;
  }
}

// Reinitialize the hit matrix so all values point to start of new file
void hitMatrix_reinitialize_multi(int4 queryLength, int4 maximumSubjectLength,
                                  unsigned char *startAddress, int queryNum) {
  int4 numDiagonals;
  unsigned char **minOffset, **offset, **maxOffset;

  if (encoding_alphabetType == encoding_nucleotide) {
    // Calculate number of diagonals that will be required during search
    numDiagonals = 1;
    while (numDiagonals < queryLength + parameters_wordSize) {
      numDiagonals <<= 1;
    }

    minOffset = hitMatrix_furthest_multi[queryNum];
  }
  // Use less memory efficient but faster hit matrix for protein
  else {
    // Maximum number of diagonals that will be required during search
    numDiagonals = queryLength + maximumSubjectLength - parameters_wordSize + 1;

    minOffset = hitMatrix_furthest_multi[queryNum] - queryLength;
  }

  hitMatrix_numDiagonals = numDiagonals;

  // Start from smallest possible offset value and iterate through to largest
  offset = minOffset;
  maxOffset = minOffset + numDiagonals;

  // For each diagonal, reset furthest to address at start of file
  while (offset < maxOffset) {
    *offset = startAddress;
    offset++;
  }
}

// Free the matrix
void hitMatrix_free_multi(int queryNum) {
  if (encoding_alphabetType == encoding_protein)
    hitMatrix_furthest_multi[queryNum] -= hitMatrix_queryLength;
  free(hitMatrix_furthest_multi[queryNum]);
}
