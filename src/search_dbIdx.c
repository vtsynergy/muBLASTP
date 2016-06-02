#define _GNU_SOURCE
#include "blast.h"
#include <pthread.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include <limits.h>
#include <immintrin.h>

#define LOG_NUM_BINS 8
#define NUM_BINS (1 << LOG_NUM_BINS)
#define NUM_BINS_1 (NUM_BINS - 1)

inline int4 getCodeword_protein_query(unsigned char *codes,
        int4 wordLength) {
    int4 codeword;
    uint4 codeCount;

    codeword = 0;
    codeCount = 0;
    while (codeCount < wordLength) {
        codeword *= wordLookupDFA_numCodes;
        int code = codes[codeCount];
        if (code < wordLookupDFA_numCodes)
            codeword += code;
        codeCount++;
    }
    return codeword;
}

inline HitPair *hit_sort_radix(HitPair *selectHits1, HitPair *selectHits2, int numSecondBins, int numExtHit)
{

    int8 bitCount = 0, bitValue = 1;

    while(bitValue < numSecondBins)
    {
        bitValue *= 2;
        bitCount++;
    }

#if defined(__ICC) || defined(__INTEL_COMPILER)
    __declspec(align(64)) unsigned int BO1[NUM_BINS + 1];
    __declspec(align(64)) unsigned int BO2[NUM_BINS + 1];
#else
    unsigned int BO1[NUM_BINS + 1];
    unsigned int BO2[NUM_BINS + 1];
#endif

    int binset = 0;
    HitPair *from_bin, *to_bin;
    unsigned int *to_bin_num;

    int bit, jj;
    for(bit = 0; bit < ceil((float)bitCount/LOG_NUM_BINS); bit++)
    {
        if(binset == 0)
        {
            to_bin_num = BO2;
            from_bin = selectHits1;
            to_bin = selectHits2;
            binset = 1;
        }
        else
        {
            to_bin_num = BO1;
            from_bin = selectHits2;
            to_bin = selectHits1;
            binset = 0;
        }


        for(jj = 0; jj < NUM_BINS + 1; jj++)
        {
            to_bin_num[jj] = 0;
        }

        for(jj = 0; jj < numExtHit; jj++)
        {
            HitPair hp = from_bin[jj];
            uint4 hitIndex = hp.hitIndex;
            uint4 binId = (hitIndex >> (bit * LOG_NUM_BINS)) & NUM_BINS_1;
            to_bin_num[binId + 1]++;
        }

        for(jj = 1; jj < NUM_BINS + 1; jj++)
        {
            to_bin_num[jj] = to_bin_num[jj - 1] + to_bin_num[jj];
        }

        for(jj = 0; jj < numExtHit; jj++)
        {
            HitPair hp = from_bin[jj];
            uint4 hitIndex = hp.hitIndex;
            uint4 binId = (hitIndex >> (bit * LOG_NUM_BINS)) & NUM_BINS_1;
            to_bin[to_bin_num[binId]++] = hp;
        }
    }

    return to_bin;
}

void search_protein2hit_dbIdx_hitDetect(
        int tid, struct PSSMatrix *PSSMatrix_arr, 
        int queryNum, 
        int bid,
        uint2 *lastHits, 
        HitPair **extHitBuf,
        int *numExtHit,
        size_t *extHitsBufSize) 
{
    uint2 subjectOffset, wordLengthMinusOne;
    uint2 queryOffset;
    struct ungappedExtension *ungappedExtension;
    int4 diagonal;

    uint4 codeword;
    uint8 queryPosition;
    wordLengthMinusOne = parameters_wordSize - 1;
    unsigned char *ungappedExtension_subjectEndReached;


    int4 length = PSSMatrix_arr[queryNum].length;
    int4 longestSeqLen = proteinLookup_db_b[bid].dbIdxblk_longestSeq;
    int numHits = 0, numUngappedExt = 0, numTriggerExt = 0, numTriggerSeq = 0;
    int numSeqBlk = proteinLookup_db_b[bid].numSeqBlk;

    int maxDiag = longestSeqLen + length;

    subPos_t *subSequencePositions =
        proteinLookup_db_b[bid].subSequencePositions;
    queryPosition = 0;

    int ii, kk, jj;

    int4 numSecondBins = (((numSeqBlk + 1) * maxDiag)) + 1;

    for(ii = 0; ii < numSecondBins; ii++)
    {
        lastHits[ii] = 0xFFFF;
    }

    int numExtHit_t = 0;
    int4 numNeighbours;

    while (queryPosition < length - wordLengthMinusOne) {

        codeword = getCodeword_protein_query(
                PSSMatrix_arr[queryNum].queryCodes + queryPosition,
                wordLookupDFA_wordLength);
        numNeighbours = neighborLookup[codeword].numNeighbours;

        while (numNeighbours > 0) {
            numNeighbours--;

            int4 codeword_neighbours = neighborLookup[codeword].neighbours[numNeighbours];

            uint4 *subPositionOffset = proteinLookup_db_b[bid].subPositionOffset;

            for (jj = subPositionOffset[codeword_neighbours]; 
                    jj < subPositionOffset[codeword_neighbours + 1]; 
                    jj++) {

                subPos_t hit = subSequencePositions[jj];
                subjectOffset = hit.subOff;
                diagonal =
                    subjectOffset - queryPosition + (length - wordLengthMinusOne);
                int seqId = hit.seqId;

                uint4 secondBinId = seqId * maxDiag + diagonal;

                uint4 currHit = queryPosition;
                uint4 lastHit = lastHits[secondBinId];

                int distance = currHit - lastHit;
                if(distance >= parameters_A || lastHit == 0xFFFF)
                {
                    lastHits[secondBinId] = currHit;
                }
                else if(distance >= parameters_overlap)
                {
                    lastHits[secondBinId] = currHit;
                    HitPair hp;
                    hp.queryOffset = queryPosition;
                    hp.hitIndex = seqId * maxDiag + diagonal;
                    hp.distance = distance;

                    if(numExtHit_t + 1 >= *extHitsBufSize)
                    {
                        *extHitsBufSize *= 2;
                        *extHitBuf = (HitPair *)realloc(*extHitBuf, 
                                sizeof(HitPair) * *extHitsBufSize);
                        fprintf(stderr, "extHitBuf resized to %d\n", *extHitsBufSize);
                    }

                    (*extHitBuf)[numExtHit_t] = hp;
                    numExtHit_t++;
                }

                numHits++;
            }
        }
        queryPosition++;
    }

    *numExtHit = numExtHit_t;
}

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
        struct ungappedExtension **ungappedExtension_new, 
        BlastGapDP *dp_mem,
        BlastIntervalTree *tree, 
        BlastIntervalTree *private_tree,
        uint4 *blast_numHits_t, 
        uint4 *blast_numUngappedExtensions_t,
        uint4 *blast_numTriggerExtensions_t,
        uint4 *blast_numTriggerSequences_t, 
        BlastHSP *BlastHSP) 
{
    uint2 subjectOffset, wordLengthMinusOne;
    uint2 queryOffset;
    struct ungappedExtension *ungappedExtension;
    int4 diagonal;

    uint4 codeword;
    uint8 queryPosition;
    wordLengthMinusOne = parameters_wordSize - 1;
    unsigned char *ungappedExtension_subjectEndReached;


    int4 length = PSSMatrix_arr[queryNum].length;
    int4 longestSeqLen = proteinLookup_db_b[bid].dbIdxblk_longestSeq;
    int numHits = 0, numUngappedExt = 0, numTriggerExt = 0, numTriggerSeq = 0;
    int numSeqBlk = proteinLookup_db_b[bid].numSeqBlk;

    int maxDiag = longestSeqLen + length;

    subPos_t *subSequencePositions =
        proteinLookup_db_b[bid].subSequencePositions;
    queryPosition = 0;
    int ii, kk, jj;

    int4 numSecondBins = (((numSeqBlk + 1) * maxDiag)) + 1;

    int4 numGoodAlign = *goodAlignCount;
    int4 numGoodExtensions = *goodExtensionCount;


    HitPair *to_bin = hit_sort_radix(extHitBuf, selectHits2, numSecondBins, numExtHit);

#if 1
    int prev_hitIndex = -1;
    int prev_seqId = -1;
    int numUngappedExtSeq = 0;
    struct alignment *alignment = NULL;
    int4 lastExt = -1;
    for(jj = 0; jj < numExtHit; jj++)
    {
        HitPair hp = to_bin[jj];
        int hitIndex = hp.hitIndex;
        if(hitIndex < prev_hitIndex)
        {
            fprintf(stderr, "ERROR! sort failed prev_hitIdex = %d, hitIndex = %d\n", prev_hitIndex, hitIndex);
            exit(0);
        }


        int seqId = hp.hitIndex/maxDiag;
        int diagonal = hp.hitIndex%maxDiag;
        int queryOffset = hp.queryOffset + wordLengthMinusOne;
        int subjectOffset = queryOffset + diagonal - (length - wordLengthMinusOne);

        if(prev_seqId != seqId)
        {
#if 1
            if(numUngappedExtSeq > 0)
            {
                if(alignments_findGoodAlignments_ncbi(
                            alignment,
                            *goodExtensionBuf,
                            PSSMatrix_arr[queryNum], scoreMatrix,
                            queryNum, ungappedExtension_new, dp_mem,
                            tree, private_tree, BlastHSP))
                {
                    numGoodAlign++;
                    numGoodExtensions += numUngappedExtSeq;
                    numTriggerExt += numUngappedExtSeq;
                    numTriggerSeq++;
                }
            }
#endif
            numUngappedExtSeq = 0;
            alignment = NULL;
            lastExt = -1;
        }


        prev_seqId = seqId;

        int4 currExt = (diagonal << 16) + subjectOffset; 

        if(lastExt > currExt)
        {
            continue;
        }

        int distance = hp.distance;;

        int4 sequenceCount =
            seqId + proteinLookup_db_b[bid].seqOffset;
        unsigned char *subject = sequenceData[sequenceCount].sequence;
        unsigned char *address = subject + subjectOffset;
        unsigned char *lastHit_addr = address - distance;
        int4 subjectLength = sequenceData[sequenceCount].sequenceLength;

        char rightExtend = 0;

        int4 lastHitOffset = subjectOffset - distance;
        numUngappedExt++;

#if 1
        if(numGoodExtensions + numUngappedExtSeq >= *goodExtensionBufSize)
        {
            *goodExtensionBufSize *= 2;
            //fprintf(stderr, "goodExtensionBuf resize to %d\n", *goodAlignBufSize);
            *goodExtensionBuf = (struct ungappedExtension *)global_realloc(
                    *goodExtensionBuf, 
                    sizeof(struct ungappedExtension) * (*goodExtensionBufSize));
        }
#endif

#if 1
        struct ungappedExtension *ungappedExtension =
            ungappedExtension_extend_ncbi_multi2(
                    PSSMatrix_arr[queryNum], scoreMatrix, subject,
                    lastHitOffset + 1, subjectOffset - wordLengthMinusOne,
                    queryOffset - wordLengthMinusOne, subjectLength,
                    length, sequenceCount,
                    &ungappedExtension_subjectEndReached, queryNum,
                    *goodExtensionBuf + numGoodExtensions, &numUngappedExtSeq,
                    &rightExtend);
#endif

        // If extension scores high enough to trigger gapping
        if (ungappedExtension) {

            if(rightExtend)
            {
                subjectOffset = (ungappedExtension_subjectEndReached -
                        subject);

                lastExt = (diagonal << 16) + subjectOffset;
            }

            if (alignment == NULL) {

                uint4 descriptionStart =
                    sequenceData[sequenceCount].descriptionStart;
                uint4 encodedLength = 
                    sequenceData[sequenceCount].encodedLength;
                uint4 descriptionLength =
                    sequenceData[sequenceCount].descriptionLength;
                int4 subjectLength = 
                    sequenceData[sequenceCount].sequenceLength;

#if 1
                if(numGoodAlign >= *goodAlignBufSize)
                {
                    *goodAlignBufSize *= 2;
                    //fprintf(stderr, "goodAlignBuf resize to %d\n", *goodAlignBufSize);
                    *goodAlignBuf = (struct alignment *)global_realloc(
                            *goodAlignBuf, 
                            sizeof(struct alignment) * (*goodAlignBufSize));
                }
#endif

                alignment = *goodAlignBuf + numGoodAlign;
                alignment->descriptionLocation = descriptionStart;
                alignment->descriptionLength = descriptionLength;
                alignment->subject = subject;
                alignment->subjectLength = subjectLength;
                alignment->encodedLength = encodedLength;
                alignment->inMemorySubject = 0;
                alignment->sequenceCount = sequenceCount;
                alignment->ungappedExtensionOffset = numGoodExtensions;
                alignment->ungappedExtensions = NULL;
                alignment->gappedExtensions = NULL;
                alignment->gappedExtensionOffset = -1;
                alignment->volumnNumber = readdb_volume;
                alignment->numExtensions = 0;
                alignment->queryCount = queryNum;
            }
            alignment->numExtensions++;
        }
        prev_hitIndex = hitIndex;
    }

    if(numUngappedExtSeq > 0)
    {
#if 1
        if(alignments_findGoodAlignments_ncbi(
                    alignment,
                    *goodExtensionBuf,
                    PSSMatrix_arr[queryNum], scoreMatrix,
                    queryNum, ungappedExtension_new, dp_mem,
                    tree, private_tree, BlastHSP))
        {
            numGoodAlign++;
            numGoodExtensions += numUngappedExtSeq;

            numTriggerExt += numUngappedExtSeq;
            numTriggerSeq++;
        }
#endif
    }
#endif

    (*goodAlignCount) = numGoodAlign;
    (*goodExtensionCount) = numGoodExtensions;

    (*blast_numHits_t) += numHits;
    (*blast_numUngappedExtensions_t) += numUngappedExt;
    (*blast_numTriggerExtensions_t) += numTriggerExt;
    (*blast_numTriggerSequences_t) += numTriggerSeq;

}

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
        uint4 *blast_numTriggerSequences_t, BlastHSP *BlastHSP) 
{
    uint2 subjectOffset, wordLengthMinusOne;
    uint2 queryOffset;
    struct ungappedExtension *ungappedExtension;
    int4 diagonal;

    uint4 codeword;
    uint8 queryPosition;
    wordLengthMinusOne = parameters_wordSize - 1;
    unsigned char *ungappedExtension_subjectEndReached;


    int4 length = PSSMatrix_arr[queryNum].length;
    int4 longestSeqLen = proteinLookup_db_b[bid].dbIdxblk_longestSeq;
    int numHits = 0, numUngappedExt = 0, numTriggerExt = 0, numTriggerSeq = 0;
    int numSeqBlk = proteinLookup_db_b[bid].numSeqBlk;

    int maxDiag = longestSeqLen + length;

    subPos_t *subSequencePositions =
        proteinLookup_db_b[bid].subSequencePositions;
    queryPosition = 0;

    int ii, kk, jj;

    int4 numSecondBins = (((numSeqBlk + 1) * maxDiag)) + 1;

    int4 numGoodAlign = *goodAlignCount;
    int4 numGoodExtensions = *goodExtensionCount;

    for(ii = 0; ii < numSecondBins; ii++)
    {
        lastHits[ii] = 0xFFFF;
    }

    int numExtHit = 0;
    int4 numNeighbours;

    while (queryPosition < length - wordLengthMinusOne) {

        codeword = getCodeword_protein_query(
                PSSMatrix_arr[queryNum].queryCodes + queryPosition,
                wordLookupDFA_wordLength);
        numNeighbours = neighborLookup[codeword].numNeighbours;

        while (numNeighbours > 0) {
            numNeighbours--;

            int4 codeword_neighbours = neighborLookup[codeword].neighbours[numNeighbours];

            uint4 *subPositionOffset = proteinLookup_db_b[bid].subPositionOffset;

            for (jj = subPositionOffset[codeword_neighbours]; 
                    jj < subPositionOffset[codeword_neighbours + 1]; 
                    jj++) {
                subPos_t hit = subSequencePositions[jj];
                subjectOffset = hit.subOff;
                diagonal =
                    subjectOffset - queryPosition + (length - wordLengthMinusOne);
                int seqId = hit.seqId;

                uint4 secondBinId = seqId * maxDiag + diagonal;

                uint4 currHit = queryPosition;
                uint4 lastHit = lastHits[secondBinId];

                int distance = currHit - lastHit;
                if(distance >= parameters_A || lastHit == 0xFFFF)
                {
                    lastHits[secondBinId] = currHit;
                }
                else if(distance >= parameters_overlap)
                {
                    lastHits[secondBinId] = currHit;
                    HitPair hp;
                    hp.queryOffset = queryPosition;
                    hp.hitIndex = seqId * maxDiag + diagonal;
                    hp.distance = distance;
                    selectHits1[numExtHit] = hp;
                    numExtHit++;
                }

                numHits++;
            }
        }
        queryPosition++;
    }

    //if(numExtHit > MAX_NUM_UNGAPPED_EXT)
    //{
        //fprintf(stderr, "%d, %d] ERROR! numExtHit = %d\n", tid, queryNum, numExtHit);
        //exit(0);
    //}

    HitPair *to_bin = hit_sort_radix(selectHits1, selectHits2, numSecondBins, numExtHit);

#if 1
    int prev_hitIndex = -1;
    int prev_seqId = -1;
    int numUngappedExtSeq = 0;
    struct alignment *alignment = NULL;
    int4 lastExt = -1;
    for(jj = 0; jj < numExtHit; jj++)
    {
        HitPair hp = to_bin[jj];
        int hitIndex = hp.hitIndex;
        if(hitIndex < prev_hitIndex)
        {
            fprintf(stderr, "ERROR! sort failed prev_hitIdex = %d, hitIndex = %d\n", prev_hitIndex, hitIndex);
            exit(0);
        }


        int seqId = hp.hitIndex/maxDiag;
        int diagonal = hp.hitIndex%maxDiag;
        int queryOffset = hp.queryOffset + wordLengthMinusOne;
        int subjectOffset = queryOffset + diagonal - (length - wordLengthMinusOne);

        if(prev_seqId != seqId)
        {
#if 1
            if(numUngappedExtSeq > 0)
            {
                if(alignments_findGoodAlignments_ncbi(
                            alignment,
                            *goodExtensionBuf,
                            PSSMatrix_arr[queryNum], scoreMatrix,
                            queryNum, ungappedExtension_new, dp_mem,
                            tree, private_tree, BlastHSP))
                {
                    numGoodAlign++;
                    numGoodExtensions += numUngappedExtSeq;
                    numTriggerExt += numUngappedExtSeq;
                    numTriggerSeq++;
                }
            }
#endif
            numUngappedExtSeq = 0;
            alignment = NULL;
            lastExt = -1;
        }


        prev_seqId = seqId;

        int4 currExt = (diagonal << 16) + subjectOffset; 

        if(lastExt > currExt)
        {
            continue;
        }

        int distance = hp.distance;;

        int4 sequenceCount =
            seqId + proteinLookup_db_b[bid].seqOffset;
        unsigned char *subject = sequenceData[sequenceCount].sequence;
        unsigned char *address = subject + subjectOffset;
        unsigned char *lastHit_addr = address - distance;
        int4 subjectLength = sequenceData[sequenceCount].sequenceLength;

        char rightExtend = 0;

        int4 lastHitOffset = subjectOffset - distance;
        numUngappedExt++;

#if 1
        if(numGoodExtensions + numUngappedExtSeq >= *goodExtensionBufSize)
        {
            *goodExtensionBufSize *= 2;
            //fprintf(stderr, "goodExtensionBuf resize to %d\n", *goodAlignBufSize);
            *goodExtensionBuf = (struct ungappedExtension *)global_realloc(
                    *goodExtensionBuf, 
                    sizeof(struct ungappedExtension) * (*goodExtensionBufSize));
        }
#endif

#if 1
        struct ungappedExtension *ungappedExtension =
            ungappedExtension_extend_ncbi_multi2(
                    PSSMatrix_arr[queryNum], scoreMatrix, subject,
                    lastHitOffset + 1, subjectOffset - wordLengthMinusOne,
                    queryOffset - wordLengthMinusOne, subjectLength,
                    length, sequenceCount,
                    &ungappedExtension_subjectEndReached, queryNum,
                    *goodExtensionBuf + numGoodExtensions, &numUngappedExtSeq,
                    &rightExtend);
#endif

        // If extension scores high enough to trigger gapping
        if (ungappedExtension) {

            if(rightExtend)
            {
                subjectOffset = (ungappedExtension_subjectEndReached -
                        subject);

                lastExt = (diagonal << 16) + subjectOffset;
            }

            if (alignment == NULL) {

                uint4 descriptionStart =
                    sequenceData[sequenceCount].descriptionStart;
                uint4 encodedLength = 
                    sequenceData[sequenceCount].encodedLength;
                uint4 descriptionLength =
                    sequenceData[sequenceCount].descriptionLength;
                int4 subjectLength = 
                    sequenceData[sequenceCount].sequenceLength;

#if 1
                if(numGoodAlign >= *goodAlignBufSize)
                {
                    *goodAlignBufSize *= 2;
                    //fprintf(stderr, "goodAlignBuf resize to %d\n", *goodAlignBufSize);
                    *goodAlignBuf = (struct alignment *)global_realloc(
                            *goodAlignBuf, 
                            sizeof(struct alignment) * (*goodAlignBufSize));
                }
#endif

                alignment = *goodAlignBuf + numGoodAlign;
                alignment->descriptionLocation = descriptionStart;
                alignment->descriptionLength = descriptionLength;
                alignment->subject = subject;
                alignment->subjectLength = subjectLength;
                alignment->encodedLength = encodedLength;
                alignment->inMemorySubject = 0;
                alignment->sequenceCount = sequenceCount;
                alignment->ungappedExtensionOffset = numGoodExtensions;
                alignment->ungappedExtensions = NULL;
                alignment->gappedExtensions = NULL;
                alignment->gappedExtensionOffset = -1;
                alignment->volumnNumber = readdb_volume;
                alignment->numExtensions = 0;
                alignment->queryCount = queryNum;
            }
            alignment->numExtensions++;
        }
        prev_hitIndex = hitIndex;
    }

    if(numUngappedExtSeq > 0)
    {
#if 1
        if(alignments_findGoodAlignments_ncbi(
                    alignment,
                    *goodExtensionBuf,
                    PSSMatrix_arr[queryNum], scoreMatrix,
                    queryNum, ungappedExtension_new, dp_mem,
                    tree, private_tree, BlastHSP))
        {
            numGoodAlign++;
            numGoodExtensions += numUngappedExtSeq;

            numTriggerExt += numUngappedExtSeq;
            numTriggerSeq++;
        }
#endif
    }
#endif

    (*goodAlignCount) = numGoodAlign;
    (*goodExtensionCount) = numGoodExtensions;

    (*blast_numHits_t) += numHits;
    (*blast_numUngappedExtensions_t) += numUngappedExt;
    (*blast_numTriggerExtensions_t) += numTriggerExt;
    (*blast_numTriggerSequences_t) += numTriggerSeq;

}

void search_protein2hit_dbIdx_lasthit_radix_cp(
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
        uint4 *blast_numTriggerSequences_t, BlastHSP *BlastHSP) 
{
    uint2 subjectOffset, wordLengthMinusOne;
    uint2 queryOffset;
    struct ungappedExtension *ungappedExtension;
    int4 diagonal;

    uint4 codeword;
    uint8 queryPosition;
    wordLengthMinusOne = parameters_wordSize - 1;
    unsigned char *ungappedExtension_subjectEndReached;


    int4 length = PSSMatrix_arr[queryNum].length;
    int4 longestSeqLen = proteinLookup_db_blk_cp[bid].dbIdxblk_longestSeq;
    int numHits = 0, numUngappedExt = 0, numTriggerExt = 0, numTriggerSeq = 0;
    int numSeqBlk = proteinLookup_db_blk_cp[bid].numSeqBlk;

    int maxDiag = longestSeqLen + length;

    //subPos_t *subSequencePositions =
    //proteinLookup_db_blk_cp[bid].subSequencePositions;
    queryPosition = 0;
    int ii, kk, jj;

    int4 numSecondBins = (((numSeqBlk + 1) * maxDiag)) + 1;

    int4 numGoodAlign = *goodAlignCount;
    int4 numGoodExtensions = *goodExtensionCount;

    for(ii = 0; ii < numSecondBins; ii++)
    {
        lastHits[ii] = 0xFFFF;
    }

    int numExtHit = 0;
    int4 numNeighbours;

    while (queryPosition < length - wordLengthMinusOne) {

        codeword = getCodeword_protein_query(
                PSSMatrix_arr[queryNum].queryCodes + queryPosition,
                wordLookupDFA_wordLength);
        numNeighbours = neighborLookup[codeword].numNeighbours;

        while (numNeighbours > 0) {
            numNeighbours--;

            int4 codeword_neighbours = neighborLookup[codeword].neighbours[numNeighbours];

            //uint4 *subPositionOffset = proteinLookup_db_b[bid].subPositionOffset;

            struct initialWord_protein_db_cp *initialWord;
            initialWord = proteinLookup_db_blk_cp[bid].proteinLookup_db_cp + codeword_neighbours;
            int subOff = 0;
            int numPosWord = 0;
            uint2 *seqIds = initialWord->seqId;
            for (jj = 0; jj < initialWord->numSubOff; jj++) {
                subOff += initialWord->subOff[jj].subOff;
                subjectOffset = subOff;
                diagonal =
                    subjectOffset - queryPosition + (length - wordLengthMinusOne);
                int kk;
                int numPos = initialWord->subOff[jj].numPos;

                for(kk = 0; kk < numPos; kk++)
                {
                    //uint4 hit = (subjectOffset << H_BITS) + seqIds[kk]; 
                    diagonal =
                        subjectOffset - queryPosition + (length - wordLengthMinusOne);

                    int seqId = seqIds[kk];

                    uint4 secondBinId = seqId * maxDiag + diagonal;

                    uint4 currHit = queryPosition;
                    uint4 lastHit = lastHits[secondBinId];

                    int distance = currHit - lastHit;
                    if(distance >= parameters_A || lastHit == 0xFFFF)
                    {
                        lastHits[secondBinId] = currHit;
                    }
                    else if(distance >= parameters_overlap)
                    {
                        lastHits[secondBinId] = currHit;
                        HitPair hp;
                        hp.queryOffset = queryPosition;
                        hp.hitIndex = seqId * maxDiag + diagonal;
                        hp.distance = distance;
                        selectHits1[numExtHit] = hp;
                        numExtHit++;
                    }
                }
                seqIds += numPos;
                numPosWord += numPos;
                //binOffset[diagonal + 1]+=numPos;
            }
        }
        queryPosition++;
    }

    //if(numExtHit > MAX_NUM_UNGAPPED_EXT)
    //{
        //fprintf(stderr, "%d, %d] ERROR! numExtHit = %d\n", tid, queryNum, numExtHit);
        //exit(0);
    //}

    HitPair *to_bin = hit_sort_radix(selectHits1, selectHits2, numSecondBins, numExtHit);

#if 1
    int prev_hitIndex = -1;
    int prev_seqId = -1;
    int numUngappedExtSeq = 0;
    struct alignment *alignment = NULL;
    int4 lastExt = -1;
    for(jj = 0; jj < numExtHit; jj++)
    {
        HitPair hp = to_bin[jj];
        int hitIndex = hp.hitIndex;
        if(hitIndex < prev_hitIndex)
        {
            fprintf(stderr, "ERROR! sort failed prev_hitIdex = %d, hitIndex = %d\n", prev_hitIndex, hitIndex);
            exit(0);
        }


        int seqId = hp.hitIndex/maxDiag;
        int diagonal = hp.hitIndex%maxDiag;
        int queryOffset = hp.queryOffset + wordLengthMinusOne;
        int subjectOffset = queryOffset + diagonal - (length - wordLengthMinusOne);

        if(prev_seqId != seqId)
        {
#if 1
            if(numUngappedExtSeq > 0)
            {
                if(alignments_findGoodAlignments_ncbi(
                            alignment,
                            *goodExtensionBuf,
                            PSSMatrix_arr[queryNum], scoreMatrix,
                            queryNum, ungappedExtension_new, dp_mem,
                            tree, private_tree, BlastHSP))
                {
                    numGoodAlign++;
                    numGoodExtensions += numUngappedExtSeq;
                    numTriggerExt += numUngappedExtSeq;
                    numTriggerSeq++;
                }
            }
#endif
            numUngappedExtSeq = 0;
            alignment = NULL;
            lastExt = -1;
        }


        prev_seqId = seqId;

        int4 currExt = (diagonal << 16) + subjectOffset; 

        if(lastExt > currExt)
        {
            continue;
        }

        int distance = hp.distance;;

        int4 sequenceCount =
            seqId + proteinLookup_db_blk_cp[bid].seqOffset;
        unsigned char *subject = sequenceData[sequenceCount].sequence;
        unsigned char *address = subject + subjectOffset;
        unsigned char *lastHit_addr = address - distance;
        int4 subjectLength = sequenceData[sequenceCount].sequenceLength;

        char rightExtend = 0;

        int4 lastHitOffset = subjectOffset - distance;
        numUngappedExt++;

#if 1
        if(numGoodExtensions + numUngappedExtSeq >= *goodExtensionBufSize)
        {
            *goodExtensionBufSize *= 2;
            //fprintf(stderr, "goodExtensionBuf resize to %d\n", *goodAlignBufSize);
            *goodExtensionBuf = (struct ungappedExtension *)global_realloc(
                    *goodExtensionBuf, 
                    sizeof(struct ungappedExtension) * (*goodExtensionBufSize));
        }
#endif

#if 1
        struct ungappedExtension *ungappedExtension =
            ungappedExtension_extend_ncbi_multi2(
                    PSSMatrix_arr[queryNum], scoreMatrix, subject,
                    lastHitOffset + 1, subjectOffset - wordLengthMinusOne,
                    queryOffset - wordLengthMinusOne, subjectLength,
                    length, sequenceCount,
                    &ungappedExtension_subjectEndReached, queryNum,
                    *goodExtensionBuf + numGoodExtensions, &numUngappedExtSeq,
                    &rightExtend);
#endif

        // If extension scores high enough to trigger gapping
        if (ungappedExtension) {

            if(rightExtend)
            {
                subjectOffset = (ungappedExtension_subjectEndReached -
                        subject);

                lastExt = (diagonal << 16) + subjectOffset;
            }

            if (alignment == NULL) {

                uint4 descriptionStart =
                    sequenceData[sequenceCount].descriptionStart;
                uint4 encodedLength = 
                    sequenceData[sequenceCount].encodedLength;
                uint4 descriptionLength =
                    sequenceData[sequenceCount].descriptionLength;
                int4 subjectLength = 
                    sequenceData[sequenceCount].sequenceLength;

#if 1
                if(numGoodAlign >= *goodAlignBufSize)
                {
                    *goodAlignBufSize *= 2;
                    //fprintf(stderr, "goodAlignBuf resize to %d\n", *goodAlignBufSize);
                    *goodAlignBuf = (struct alignment *)global_realloc(
                            *goodAlignBuf, 
                            sizeof(struct alignment) * (*goodAlignBufSize));
                }
#endif

                alignment = *goodAlignBuf + numGoodAlign;
                alignment->descriptionLocation = descriptionStart;
                alignment->descriptionLength = descriptionLength;
                alignment->subject = subject;
                alignment->subjectLength = subjectLength;
                alignment->encodedLength = encodedLength;
                //alignment->joinChecked = 0;
                alignment->inMemorySubject = 0;
                //alignment->numUnpackRegions = 0;
                //alignment->cluster = 0;
                alignment->sequenceCount = sequenceCount;
                alignment->ungappedExtensionOffset = numGoodExtensions;
                alignment->ungappedExtensions = NULL;
                alignment->gappedExtensions = NULL;
                alignment->gappedExtensionOffset = -1;
                alignment->volumnNumber = readdb_volume;
                //alignment->unpackRegions = NULL;
                //alignment->edits = NULL;
                alignment->numExtensions = 0;
                alignment->queryCount = queryNum;
            }
            alignment->numExtensions++;
        }
        prev_hitIndex = hitIndex;
    }

    if(numUngappedExtSeq > 0)
    {
#if 1
        if(alignments_findGoodAlignments_ncbi(
                    alignment,
                    *goodExtensionBuf,
                    PSSMatrix_arr[queryNum], scoreMatrix,
                    queryNum, ungappedExtension_new, dp_mem,
                    tree, private_tree, BlastHSP))
        {
            numGoodAlign++;
            numGoodExtensions += numUngappedExtSeq;

            numTriggerExt += numUngappedExtSeq;
            numTriggerSeq++;
        }
#endif
    }
#endif

    (*goodAlignCount) = numGoodAlign;
    (*goodExtensionCount) = numGoodExtensions;

    (*blast_numHits_t) += numHits;
    (*blast_numUngappedExtensions_t) += numUngappedExt;
    (*blast_numTriggerExtensions_t) += numTriggerExt;
    (*blast_numTriggerSequences_t) += numTriggerSeq;

}
