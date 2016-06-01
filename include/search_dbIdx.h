#ifndef _search_db_
#define _search_db_

typedef struct hit_pair_struct
{
    uint4 hitIndex;
    uint2 queryOffset:10, distance:6;
    /*uint2 queryOffset, distance;*/
}HitPair;


void search_protein2hit_dbIdx_lasthit_radix(
        int tid, struct PSSMatrix *PSSMatrix_arr, 
        struct scoreMatrix scoreMatrix,
        int queryNum, 
        struct sequenceData *sequenceData, 
        int bid,
        uint2 *lastHits, 
        HitPair *selectHits1, HitPair *selectHits2,
        struct alignment **goodAlignBuf, 
        int4 *goodAlignCount, 
        size_t *goodAlignBufSize,
        struct ungappedExtension **goodExtensionBuf, 
        int4 *goodExtensionCount, 
        size_t *goodExtensionBufSize,
        struct ungappedExtension **ungappedExtension_new, BlastGapDP *dp_mem,
        BlastIntervalTree *tree, BlastIntervalTree *private_tree,
        uint4 *blast_numHits_t, uint4 *blast_numUngappedExtensions_t,
        uint4 *blast_numTriggerExtensions_t,
        uint4 *blast_numTriggerSequences_t, BlastHSP *BlastHSP); 


void search_protein2hit_dbIdx_hitDetect(
        int tid, struct PSSMatrix *PSSMatrix_arr, 
        int queryNum, 
        int bid,
        uint2 *lastHits, 
        HitPair **extHitBuf,
        int *numExtHit,
        size_t *extHitsBufSize); 


void search_protein2hit_dbIdx_ungapExt(
        int tid, struct PSSMatrix *PSSMatrix_arr, 
        struct scoreMatrix scoreMatrix,
        int queryNum, 
        struct sequenceData *sequenceData, 
        int bid,
        HitPair *extHitBuf,
        int numExtHit,
        HitPair *selectHits2,
        struct alignment **goodAlignBuf, 
        int4 *goodAlignCount, 
        size_t *goodAlignBufSize,
        struct ungappedExtension **goodExtensionBuf, 
        int4 *goodExtensionCount, 
        size_t *goodExtensionBufSize,
        struct ungappedExtension **ungappedExtension_new, BlastGapDP *dp_mem,
        BlastIntervalTree *tree, BlastIntervalTree *private_tree,
        uint4 *blast_numHits_t, uint4 *blast_numUngappedExtensions_t,
        uint4 *blast_numTriggerExtensions_t,
        uint4 *blast_numTriggerSequences_t, BlastHSP *BlastHSP); 

#endif

