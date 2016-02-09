#ifndef _search_db_
#define _search_db_

typedef struct hit_pair_struct
{
    uint4 hitIndex;
    uint2 queryOffset;
    uint2 distance;
}HitPair;

void search_protein2hit_lookup_db(
        int tid,
        struct PSSMatrix *PSSMatrix_arr,
        struct scoreMatrix scoreMatrix,
        int queryNum,
        struct sequenceData *sequenceData,
        int dbIdxBlockNum,
        uint4 *binOffset,
        hit_t *hits,
        uint4 *numHitsPerQPos,
        hit_t *binnedHits,
        hit_t *lastHits,
        hit_t *seqHits,
        uint4 *numExtHit,
        /*uint4 *allocExtHit,*/
        struct alignment *goodAlignBuf, int *goodAlignCount,
        struct ungappedExtension *goodExtensionBuf, int *goodExtensionCount,
        struct ungappedExtension **ungappedExtension_new, BlastGapDP *dp_mem,
        BlastIntervalTree *tree, BlastIntervalTree *private_tree,
        BlastHSP *BlastHSP,
        uint4 *blast_numHits_multi_t,
        uint4 *blast_numUngappedExtensions_multi_t,
        uint4 *blast_numTriggerExtensions_multi_t,
        uint4 *blast_numTriggerSequences_multi_t,
        Profile *profile
        ); 

#endif

