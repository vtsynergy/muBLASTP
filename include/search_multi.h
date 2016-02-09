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



#ifndef _search_multi_
#define _search_multi_


// Search a protein database using 2-hit extension mode
void search_protein2hit_multithread(struct PSSMatrix *PSSMatrix_arr, struct sequenceData* sequenceData, uint4 numSequences, int numQuery);
void search_protein2hit_multi(struct PSSMatrix *PSSMatrix_arr, struct sequenceData* sequenceData, uint4 numSequences, int queryNum);

void search_protein2hit_lookup_multi(struct PSSMatrix *PSSMatrix_arr,
                              struct sequenceData *sequenceData,
                              uint4 numSequences, int numQuery, 
                              hit_t *firstBin, hit_t *secondBin, 
                              uint4 *binOffset,
                              uint4 *numExtHit,
                              uint4 *numHit,
                              hit_t *lastHits,
                              struct scoreMatrix scoreMatrix,
                              struct ungappedExtension *goodExtensionBuf,
                              int *goodExtensionCount,
                              struct alignment *goodAlignBuf, 
                              int *goodAlignCount,
                              struct ungappedExtension **ungappedExtension_new, 
                              BlastGapDP *dp_mem,
                              BlastIntervalTree *tree,
                              BlastIntervalTree *private_tree,
                              BlastHSP *BlastHSP
                              );

#endif
