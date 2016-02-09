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



#ifndef _hitMatrix_multi_
#define _hitMatrix_multi_
extern unsigned char** hitMatrix_furthest_multi[];

// Initialize the hitMatrix by declaring memory for the maximum number
// of diagonals required by a subject sequence
void hitMatrix_initialize_multi(int4 queryLength, int4 maximumSubjectLength, unsigned char* startAddress, int queryNum);

// Reinitialize the hit matrix so all values point to start of new file
void hitMatrix_reinitialize_multi(int4 queryLength, int4 maximumSubjectLength, unsigned char* startAddress, int queryNum);

// Free the matrix
void hitMatrix_free_multi(int queryNum);


#endif
