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



#ifndef _search_queryIdx_
#define _search_queryIdx_

void search_protein2hit_queryIdx(
        int thread_id,
        int subSequenceId,
        int block_id,
        struct PSSMatrix *PSSMatrix_arr,
        struct sequenceData *sequenceData,
        uint4 numSequences, int numQuery, 
        hit_t *secondBin, 
        uint4 *numExtHit,
        uint2 *lastHits,
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
