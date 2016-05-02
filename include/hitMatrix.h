#ifndef _hitMatrix_
#define _hitMatrix_

unsigned char** hitMatrix_furthest;
int4 hitMatrix_queryLength;
uint4 hitMatrix_diagonalMask;

// Initialize the hitMatrix by declaring memory for the maximum number
// of diagonals required by a subject sequence
void hitMatrix_initialize(int4 queryLength, int4 maximumSubjectLength, unsigned char* startAddress);

// Reinitialize the hit matrix so all values point to start of new file
void hitMatrix_reinitialize(int4 queryLength, int4 maximumSubjectLength, unsigned char* startAddress);

// Free the matrix
void hitMatrix_free();

#endif

