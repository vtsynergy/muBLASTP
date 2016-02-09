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


#ifndef _wordLookupDFA_multi__
#define _wordLookupDFA_multi__

//Multi
extern uint2* wordLookupDFA_additionalQueryPositions_multi[];
extern int4 wordLookupDFA_numAdditionalQueryPositions_multi[];
extern struct group *wordLookupDFA_groups_multi[];
extern int4 wordLookupDFA_numBlocks_multi[];


// Build the word-lookup structure
void wordLookupDFA_assign(int queryNum);


// Free memory used by the word lookup table
void wordLookupDFA_free_multi(int queryNum);


void wordLookupDFA_print_multi(int queryNum);

#endif
