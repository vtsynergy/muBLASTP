#define _GNU_SOURCE
#include "blast.h"
#include <pthread.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include <limits.h>
#include <immintrin.h>
//#include <zmmintrin.h>

//#define PREFETCH
#define H_BITS 15
#define L_MASK 0x7fff

pthread_barrier_t barrier;

inline int4 wordLookupDFA_getCodeword_protein_query(unsigned char *codes,
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


void search_protein2hit_lookup_db_cp(
        int tid, struct PSSMatrix *PSSMatrix_arr, struct scoreMatrix scoreMatrix,
        int queryNum, struct sequenceData *sequenceData, int dbIdxBlockNum,
        uint4 *binOffset, hit_t *hits, uint4 *numHitsPerQPos, hit_t *binnedHits,
        hit_t *lastHits, hit_t *seqHits, uint4 *numExtHit, 
        struct alignment *goodAlignBuf, int *goodAlignCount,
        struct ungappedExtension *goodExtensionBuf, int *goodExtensionCount,
        struct ungappedExtension **ungappedExtension_new, BlastGapDP *dp_mem,
        BlastIntervalTree *tree, BlastIntervalTree *private_tree,
        BlastHSP *BlastHSP,
        uint4 *blast_numHits_multi_t, uint4 *blast_numUngappedExtensions_multi_t,
        uint4 *blast_numTriggerExtensions_multi_t,
        uint4 *blast_numTriggerSequences_multi_t, Profile *profile) {
    int4 subjectOffset, wordLengthMinusOne, count = 0;
    uint2 queryOffset;
    struct ungappedExtension *ungappedExtension;
    int4 diagonal;

    uint4 codeword, queryPosition;
    int4 numNeighbours;
    struct initialWord_protein_db_cp *initialWord;
    wordLengthMinusOne = parameters_wordSize - 1;
    unsigned char *ungappedExtension_subjectEndReached_t;

#ifdef PROFILE
    RDTSC_INIT;
#endif

#ifndef NEIGHBOR_INDEX
    struct neighbour *neighbours = NULL;
    neighbours = (struct neighbour *)global_malloc(sizeof(struct neighbour) *
            proteinLookup_numWords);
#endif

    int4 numGoodAlign = *goodAlignCount;
    int4 numGoodExtensions = *goodExtensionCount;

    int4 length = PSSMatrix_arr[queryNum].length;
    int numHits = 0, numUngappedExt = 0, numTriggerExt = 0, numTriggerSeq = 0;
    //uint32_t *subSequencePositions =
    //proteinLookup_db_blk_cp[dbIdxBlockNum].subSequencePositions;
    int numSeqBlk = proteinLookup_db_blk_cp[dbIdxBlockNum].numSeqBlk;
    int4 maxNumHitPerSeq = MAX_NUM_HIT_PER_BLK / numSeqBlk; 
    queryPosition = 0;
    // printf("longestSeqBlk: %d\n", dbIdxblk_longestSeq[dbIdxBlockNum]);
    int numBinBlk =
        proteinLookup_db_blk_cp[dbIdxBlockNum].dbIdxblk_longestSeq + length + 1;

#ifdef PROFILE
    // long long start = __rdtsc();
    RDTSC_START;
#endif
    memset(binOffset, 0, sizeof(uint4) * numBinBlk);
    memset(numHitsPerQPos, 0, sizeof(uint4) * length);
    while (queryPosition < length - wordLengthMinusOne) {

#ifdef NEIGHBOR_INDEX
        codeword = wordLookupDFA_getCodeword_protein_query(
                PSSMatrix_arr[queryNum].queryCodes + queryPosition,
                wordLookupDFA_wordLength);
        numNeighbours = neighborLookup[codeword].numNeighbours;
#else
        numNeighbours = 0;
        wordLookupDFA_getNeighbours(PSSMatrix_arr[queryNum], queryPosition,
                &numNeighbours, neighbours);
#endif

        while (numNeighbours > 0) {
            numNeighbours--;

#ifdef NEIGHBOR_INDEX
            int4 codeword_t = neighborLookup[codeword].neighbours[numNeighbours];
#else
            int4 codeword_t = neighbours[numNeighbours].codeword;
#endif

            initialWord = proteinLookup_db_blk_cp[dbIdxBlockNum].proteinLookup_db_cp + codeword_t;
            uint4 jj;
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
                    hits[numHits++] = (subjectOffset << H_BITS) + seqIds[kk]; 
                }
                seqIds += numPos;
                numPosWord += numPos;
                binOffset[diagonal + 1]+=numPos;
            }
            numHitsPerQPos[queryPosition] += initialWord->numPos;
        }
        queryPosition++;
    }

    int jj;
    for (jj = 1; jj < numBinBlk; jj++) {
        binOffset[jj] = binOffset[jj - 1] + binOffset[jj];
    }

#ifdef PROFILE
    RDTSC_END;
    GET_RDTSC_ATOMIC(&profile->blast_pseudoCycle);
    //maxBinSize = max(maxBinSize, binOffset[numBinBlk - 1]);
    if (binOffset[numBinBlk - 1] > MAX_HITS_PER_SEQ) {
        printf("Error: bin overflow: %d\n", binOffset[numBinBlk - 1]);
        exit(1);
    }
    RDTSC_START;
#endif

#define PFD 48
    numHits = 0;
    for (queryPosition = 0; queryPosition < (length - wordLengthMinusOne);
            queryPosition++) {
        uint4 hitsPerPos = numHitsPerQPos[queryPosition];
        for (jj = 0; jj < hitsPerPos; jj++) {
            subjectOffset = hits[numHits] >> H_BITS;
            diagonal =
                subjectOffset - queryPosition + (length - wordLengthMinusOne);
            //int seqId = hits[numHits] & 0xffff;
            //fprintf(stderr, "queryOff: %d seqId: %d subOffset: %d diag: %d binOffset: %d\n", queryPosition, seqId, subjectOffset, diagonal, binOffset[diagonal]);
            binnedHits[binOffset[diagonal]] = hits[numHits];
            if ((binOffset[diagonal] % 16) == 0) {
                _mm_prefetch((char *)(&binnedHits[binOffset[diagonal] + PFD]),
                        _MM_HINT_NTA);
            }
            binOffset[diagonal]++;
            numHits++;
        }
    }

#ifdef PROFILE
    RDTSC_END;
    GET_RDTSC_ATOMIC(&profile->blast_hitDetectCycle);

    RDTSC_START;
#endif

    memset(numExtHit, 0, sizeof(uint4) * numSeqBlk);
    memset(lastHits, 255, sizeof(hit_t) * numSeqBlk);

    int kk = 0;
    for (jj = 0; jj < numBinBlk - 1; jj++) {
        //fprintf(stderr, "binOffset: %d\n", binOffset[jj]);
        for (; kk < binOffset[jj]; kk++) {
            hit_t hit = binnedHits[kk];
            //fprintf(stderr, "diag: %d seqId: %d subOffset: %d\n", jj, hit[1], hit[0]);
            //uint2 *hit = (uint2 *)(binnedHits + kk);
            // printf("%d %d %d\n", jj, hit[0], hit[1]);
            uint2 seqId = hit & L_MASK;
#ifdef PREFETCH
            uint32_t *h0 = lastHits + seqId;
            _mm_prefetch((char *)h0, _MM_HINT_NTA);
#endif

            hit_t lastHit = lastHits[seqId];
            subjectOffset = hit >> H_BITS;
            queryOffset = subjectOffset - jj + length; 
            hit_t currHit = queryOffset + (jj << H_BITS);
            //hit_t currHit = (hit >> 16) + (jj << 16);
            int distance = currHit - lastHit;
            if (distance >= parameters_A || lastHit == 0xffffffff) {
                lastHits[seqId] = currHit;
            } else if (distance >= parameters_overlap) {
                numExtHit[seqId] += 2;

                if (numExtHit[seqId] > maxNumHitPerSeq) {
                    fprintf(stderr, "numExtHit: %d > maxNumHitPerSeq: %d\n", numExtHit[seqId], maxNumHitPerSeq );
                    exit(1);
                }

                seqHits[seqId * maxNumHitPerSeq + numExtHit[seqId] - 2] = lastHit;
                seqHits[seqId * maxNumHitPerSeq + numExtHit[seqId] - 1] = currHit;
                lastHits[seqId] = currHit;
            }
        }
    }

#ifdef PROFILE
    RDTSC_END;
    GET_RDTSC_ATOMIC(&profile->blast_sortCycle);
#endif

#ifndef NO_STAGE2

#ifdef PROFILE
    RDTSC_START;
#endif
    for (jj = 0; jj < numSeqBlk; jj++) {
        int kk;

        if (numExtHit[jj] < 2)
            continue;

        struct alignment *alignment = NULL;

        hit_t lastExt = 0;
        //char extFlag = 0;
        hit_t prevHit = 0;
        int numUngappedExtSeq = 0;
        for (kk = 0; kk < numExtHit[jj]; kk += 2) {
            hit_t lastHit = seqHits[jj * maxNumHitPerSeq + kk];
            hit_t currHit = seqHits[jj * maxNumHitPerSeq + kk + 1];
            int distance = currHit - lastHit;

            if(currHit > lastExt)
            {

                // Not overlaping - extension triggered
                // Increment tally number of extensions
                // blast_numUngappedExtensions_multi[queryNum]++;
                numUngappedExt++;

                int seqId = jj;
                //uint2 *hit = (uint2 *)(&currHit);
                int diagonal = currHit >> H_BITS;
                queryOffset = currHit & L_MASK;
                subjectOffset = queryOffset + diagonal - (length - wordLengthMinusOne);
                //subjectOffset = currHit & 0xffff;
                int4 sequenceCount =
                    seqId + proteinLookup_db_blk_cp[dbIdxBlockNum].seqOffset;
                //subjectOffset = subjectOffset + wordLengthMinusOne;
                unsigned char *subject = sequenceData[sequenceCount].sequence;

                //queryOffset =
                //subjectOffset - diagonal + (length - wordLengthMinusOne);

                //fprintf(stderr, "%d %d %d\n", sequenceCount, subjectOffset, queryOffset);

                unsigned char *address = subject + subjectOffset;
                unsigned char *lastHit_addr = address - distance;

                int4 subjectLength = sequenceData[sequenceCount].sequenceLength;

                // Perform ungapped extension start between query/subject start/end
                // and extending outwards in each direction
                char rightExtend = 0;
#ifndef NCBI_BLAST
                struct ungappedExtension *ungappedExtension =
                    ungappedExtension_extend_multi2(
                            PSSMatrix_arr[queryNum].matrix + queryOffset, address,
                            lastHit_addr, PSSMatrix_arr[queryNum], subject, queryNum,
                            dbIdxBlockNum, &ungappedExtension_subjectEndReached_t);
#else

                int4 lastHitOffset = subjectOffset - distance;
                struct ungappedExtension *ungappedExtension =
                    ungappedExtension_extend_ncbi_multi2(
                            PSSMatrix_arr[queryNum], scoreMatrix, subject,
                            lastHitOffset + 1, subjectOffset - wordLengthMinusOne,
                            queryOffset - wordLengthMinusOne, subjectLength,
                            PSSMatrix_arr[queryNum].length, sequenceCount,
                            &ungappedExtension_subjectEndReached_t, queryNum,
                            goodExtensionBuf + numGoodExtensions, &numUngappedExtSeq,
                            &rightExtend);
#endif
                // If extension scores high enough to trigger gapping
                if (ungappedExtension) {
                    if(rightExtend)
                    {
                        subjectOffset = (ungappedExtension_subjectEndReached_t -
                                subject - wordLengthMinusOne);
                        queryOffset = subjectOffset - diagonal + (PSSMatrix_arr[queryNum].length - wordLengthMinusOne);

                        lastExt = (diagonal << H_BITS) + queryOffset;

                        //fprintf(stderr, "%d %d %d\n", sequenceCount, subjectOffset, queryOffset);
                    }
                    else
                    {
                        lastExt = currHit;
                    }

                    //numUngappedExtSeq++;
                    //extFlag = 1;
                    numTriggerExt++;
                    // blast_numTriggerExtensions_multi[queryNum]++;
                    if (alignment == NULL) {
                        uint4 descriptionStart =
                            sequenceData[sequenceCount].descriptionStart;
                        uint4 encodedLength = sequenceData[sequenceCount].encodedLength;
                        uint4 descriptionLength =
                            sequenceData[sequenceCount].descriptionLength;
                        int4 subjectLength = sequenceData[sequenceCount].sequenceLength;

                        alignment = goodAlignBuf + numGoodAlign;
                        alignment->descriptionLocation = descriptionStart;
                        alignment->descriptionLength = descriptionLength;
                        alignment->subject = subject;
                        alignment->subjectLength = subjectLength;
                        alignment->encodedLength = encodedLength;
                        alignment->joinChecked = 0;
                        alignment->inMemorySubject = 0;
                        alignment->numUnpackRegions = 0;
                        alignment->cluster = 0;
                        alignment->sequenceCount = sequenceCount;
                        alignment->ungappedExtensions = goodExtensionBuf + numGoodExtensions;
                        alignment->gappedExtensions = NULL;
                        alignment->unpackRegions = NULL;
                        alignment->edits = NULL;
                        alignment->numExtensions = 0;
                        alignment->queryCount = queryNum;
                        numTriggerSeq++;
                        //alignments_currentAlignment_multi_t =
                        //alignments_createNew_multi_db2(
                        //goodAlignBuf, &numGoodAlign,
                        //descriptionStart, descriptionLength, subject,
                        //subjectLength, encodedLength,
                        //sequenceCount);
                        //numTriggerSeq++;
                        //alignments_currentAlignment_multi_t->ungappedExtensions = goodExtensionBuf + numGoodExtensions;
                    }

                    alignment->numExtensions++;
                    //alignments_addUngappedExtension_multi_db(
                    //ungappedExtension, queryNum,
                    //alignments_currentAlignment_multi_t);
                }
                //else{extFlag = 0;}
            }
        }

#ifdef NO_STAGE3
        if(numUngappedExtSeq > 0)
        {
            numGoodAlign++;
            numGoodExtensions += numUngappedExtSeq;
        }
#else
        if(numUngappedExtSeq > 0)
        {
            if(alignments_findGoodAlignments_ncbi_multi3(
                        alignment,
                        PSSMatrix_arr[queryNum], scoreMatrix,
                        queryNum, ungappedExtension_new, dp_mem,
                        tree, private_tree, BlastHSP))
            {
                //fprintf(stderr, "numGoodAlign: %d numUngappedExtSeq: %d\n", numGoodAlign, numUngappedExtSeq);
                numGoodAlign++;
                numGoodExtensions += numUngappedExtSeq;
            }
        }
#endif
    }

#ifdef PROFILE
    RDTSC_END;
    GET_RDTSC_ATOMIC(&profile->blast_ungappedExtCycle);
#endif
#endif


    (*goodAlignCount) = numGoodAlign;
    (*goodExtensionCount) = numGoodExtensions;

    (*blast_numHits_multi_t) += numHits;
    (*blast_numUngappedExtensions_multi_t) += numUngappedExt;
    (*blast_numTriggerExtensions_multi_t) += numTriggerExt;
    (*blast_numTriggerSequences_multi_t) += numTriggerSeq;

#ifndef NEIGHBOR_INDEX
    free(neighbours);
#endif

}

void search_protein2hit_lookup_db(
        int tid, struct PSSMatrix *PSSMatrix_arr, struct scoreMatrix scoreMatrix,
        int queryNum, struct sequenceData *sequenceData, int dbIdxBlockNum,
        uint4 *binOffset, hit_t *hits, uint4 *numHitsPerQPos, hit_t *binnedHits,
        hit_t *lastHits, hit_t *seqHits, uint4 *numExtHit, 
        struct alignment *goodAlignBuf, int *goodAlignCount,
        struct ungappedExtension *goodExtensionBuf, int *goodExtensionCount,
        struct ungappedExtension **ungappedExtension_new, BlastGapDP *dp_mem,
        BlastIntervalTree *tree, BlastIntervalTree *private_tree,
        BlastHSP *BlastHSP,
        uint4 *blast_numHits_multi_t, uint4 *blast_numUngappedExtensions_multi_t,
        uint4 *blast_numTriggerExtensions_multi_t,
        uint4 *blast_numTriggerSequences_multi_t, Profile *profile) {
    uint2 subjectOffset, wordLengthMinusOne, count = 0;
    uint2 queryOffset;
    struct ungappedExtension *ungappedExtension;
    int4 diagonal;

    uint4 codeword, queryPosition;
    int4 numNeighbours;
    struct initialWord_protein_db *initialWord;
    wordLengthMinusOne = parameters_wordSize - 1;
    unsigned char *ungappedExtension_subjectEndReached_t;

#ifdef PROFILE
    RDTSC_INIT;
    int maxBinSize = 0;
#endif

#ifndef NEIGHBOR_INDEX
    struct neighbour *neighbours = NULL;
    neighbours = (struct neighbour *)global_malloc(sizeof(struct neighbour) *
            proteinLookup_numWords);
#endif

    int4 numGoodAlign = *goodAlignCount;
    int4 numGoodExtensions = *goodExtensionCount;

    int4 length = PSSMatrix_arr[queryNum].length;
    int numHits = 0, numUngappedExt = 0, numTriggerExt = 0, numTriggerSeq = 0;
    subPos_t *subSequencePositions =
        proteinLookup_db_b[dbIdxBlockNum].subSequencePositions;
    int numSeqBlk = proteinLookup_db_b[dbIdxBlockNum].numSeqBlk;
    int4 maxNumHitPerSeq = MAX_NUM_HIT_PER_BLK / numSeqBlk; 
    queryPosition = 0;
    // printf("longestSeqBlk: %d\n", dbIdxblk_longestSeq[dbIdxBlockNum]);
    int numBinBlk =
        proteinLookup_db_b[dbIdxBlockNum].dbIdxblk_longestSeq + length + 1;

#ifdef PROFILE
    // long long start = __rdtsc();
    RDTSC_START;
#endif
    memset(binOffset, 0, sizeof(uint4) * numBinBlk);
    memset(numHitsPerQPos, 0, sizeof(uint4) * length);
    while (queryPosition < length - wordLengthMinusOne) {

#ifdef NEIGHBOR_INDEX
        codeword = wordLookupDFA_getCodeword_protein_query(
                PSSMatrix_arr[queryNum].queryCodes + queryPosition,
                wordLookupDFA_wordLength);
        numNeighbours = neighborLookup[codeword].numNeighbours;
#else
        numNeighbours = 0;
        wordLookupDFA_getNeighbours(PSSMatrix_arr[queryNum], queryPosition,
                &numNeighbours, neighbours);
#endif

        while (numNeighbours > 0) {
            numNeighbours--;

#ifdef NEIGHBOR_INDEX
            int4 codeword_t = neighborLookup[codeword].neighbours[numNeighbours];
#else
            int4 codeword_t = neighbours[numNeighbours].codeword;
#endif
            uint4 jj;
            for (jj = proteinLookup_db_b[dbIdxBlockNum].subPositionOffset[codeword_t]; 
                    jj < proteinLookup_db_b[dbIdxBlockNum].subPositionOffset[codeword_t + 1]; jj++) {
                subjectOffset =
                    subSequencePositions[jj].subOff;
                diagonal =
                    subjectOffset - queryPosition + (length - wordLengthMinusOne);

                hits[numHits++] =
                    (subSequencePositions[jj].subOff << H_BITS) + subSequencePositions[jj].seqId;
                binOffset[diagonal + 1]++;
            }
            numHitsPerQPos[queryPosition] += 
                proteinLookup_db_b[dbIdxBlockNum].subPositionOffset[codeword_t + 1]
                - proteinLookup_db_b[dbIdxBlockNum].subPositionOffset[codeword_t];
        }
        queryPosition++;
    }

    int jj;
    for (jj = 1; jj < numBinBlk; jj++) {
        binOffset[jj] = binOffset[jj - 1] + binOffset[jj];
    }

#ifdef PROFILE
    RDTSC_END;
    GET_RDTSC_ATOMIC(&profile->blast_pseudoCycle);
    maxBinSize = MAX(maxBinSize, binOffset[numBinBlk - 1]);
    if (maxBinSize > MAX_HITS_PER_SEQ) {
        printf("Error: bin overflow: %d\n", maxBinSize);
        exit(1);
    }
    RDTSC_START;
#endif

#define PFD 48
    numHits = 0;
    for (queryPosition = 0; queryPosition < (length - wordLengthMinusOne);
            queryPosition++) {
        uint4 hitsPerPos = numHitsPerQPos[queryPosition];
        for (jj = 0; jj < hitsPerPos; jj++) {
            subjectOffset = hits[numHits] >> H_BITS;
            diagonal =
                subjectOffset - queryPosition + (length - wordLengthMinusOne);
            //int seqId = hits[numHits] & L_MASK;
            //fprintf(stderr, "seqId %d diag: %d subOffset: %d\n", seqId, diagonal, subjectOffset);
            binnedHits[binOffset[diagonal]] = hits[numHits];
            if ((binOffset[diagonal] % 16) == 0) {
                _mm_prefetch((char *)(&binnedHits[binOffset[diagonal] + PFD]),
                        _MM_HINT_NTA);
            }
            binOffset[diagonal]++;
            numHits++;
        }
    }

#ifdef PROFILE
    RDTSC_END;
    GET_RDTSC_ATOMIC(&profile->blast_hitDetectCycle);

    RDTSC_START;
#endif

    memset(numExtHit, 0, sizeof(uint4) * numSeqBlk);
    memset(lastHits, 255, sizeof(hit_t) * numSeqBlk);

    int kk = 0;
    for (jj = 0; jj < numBinBlk - 1; jj++) {
        for (; kk < binOffset[jj]; kk++) {
            hit_t hit = binnedHits[kk];
            uint2 seqId = hit & L_MASK;

#ifdef PREFETCH
            uint32_t *h0 = lastHits + seqId;
            _mm_prefetch((char *)h0, _MM_HINT_NTA);
#endif

            hit_t lastHit = lastHits[seqId];
            subjectOffset = hit >> H_BITS;
            queryOffset = subjectOffset - jj + length; 
            hit_t currHit = queryOffset + (jj << H_BITS);
            int distance = currHit - lastHit;
            if (distance >= parameters_A || lastHit == 0xffffffff) {
                lastHits[seqId] = currHit;
            } else if (distance >= parameters_overlap) {
                numExtHit[seqId] += 2;
                if (numExtHit[seqId] > maxNumHitPerSeq) {
                    fprintf(stderr, "numExtHit: %d > maxNumHitPerSeq: %d\n", numExtHit[seqId], maxNumHitPerSeq );
                    exit(1);
                }

                seqHits[seqId * maxNumHitPerSeq + numExtHit[seqId] - 2] = lastHit;
                seqHits[seqId * maxNumHitPerSeq + numExtHit[seqId] - 1] = currHit;
                lastHits[seqId] = currHit;
            }
        }
    }

#ifdef PROFILE
    RDTSC_END;
    GET_RDTSC_ATOMIC(&profile->blast_sortCycle);
#endif

#ifndef NO_STAGE2

#ifdef PROFILE
    RDTSC_START;
#endif
    for (jj = 0; jj < numSeqBlk; jj++) {
        int kk;

        if (numExtHit[jj] < 2)
            continue;

        struct alignment *alignment = NULL;

        hit_t lastExt = 0;
        //char extFlag = 0;
        hit_t prevHit = 0;
        int numUngappedExtSeq = 0;
        for (kk = 0; kk < numExtHit[jj]; kk += 2) {
            hit_t lastHit = seqHits[jj * maxNumHitPerSeq + kk];
            hit_t currHit = seqHits[jj * maxNumHitPerSeq + kk + 1];
            int distance = currHit - lastHit;

            if(currHit > lastExt)
            {

                // Not overlaping - extension triggered
                // Increment tally number of extensions
                // blast_numUngappedExtensions_multi[queryNum]++;
                numUngappedExt++;

                int seqId = jj;
                //uint2 *hit = (uint2 *)(&currHit);
                int diagonal = currHit >> H_BITS;
                queryOffset = currHit & L_MASK;
                subjectOffset = queryOffset + diagonal - (length - wordLengthMinusOne);
                int4 sequenceCount =
                    seqId + proteinLookup_db_b[dbIdxBlockNum].seqOffset;
                //subjectOffset = subjectOffset + wordLengthMinusOne;
                unsigned char *subject = sequenceData[sequenceCount].sequence;

                //queryOffset =
                //subjectOffset - diagonal + (length - wordLengthMinusOne);

                //fprintf(stderr, "%d %d %d\n", sequenceCount, subjectOffset, queryOffset);

                unsigned char *address = subject + subjectOffset;
                unsigned char *lastHit_addr = address - distance;

                int4 subjectLength = sequenceData[sequenceCount].sequenceLength;

                // Perform ungapped extension start between query/subject start/end
                // and extending outwards in each direction
                char rightExtend = 0;
#ifndef NCBI_BLAST
                struct ungappedExtension *ungappedExtension =
                    ungappedExtension_extend_multi2(
                            PSSMatrix_arr[queryNum].matrix + queryOffset, address,
                            lastHit_addr, PSSMatrix_arr[queryNum], subject, queryNum,
                            dbIdxBlockNum, &ungappedExtension_subjectEndReached_t);
#else

                int4 lastHitOffset = subjectOffset - distance;
                struct ungappedExtension *ungappedExtension =
                    ungappedExtension_extend_ncbi_multi2(
                            PSSMatrix_arr[queryNum], scoreMatrix, subject,
                            lastHitOffset + 1, subjectOffset - wordLengthMinusOne,
                            queryOffset - wordLengthMinusOne, subjectLength,
                            PSSMatrix_arr[queryNum].length, sequenceCount,
                            &ungappedExtension_subjectEndReached_t, queryNum,
                            goodExtensionBuf + numGoodExtensions, &numUngappedExtSeq,
                            &rightExtend);
#endif
                // If extension scores high enough to trigger gapping
                if (ungappedExtension) {
                    if(rightExtend)
                    {
                        subjectOffset = (ungappedExtension_subjectEndReached_t -
                                subject - wordLengthMinusOne);
                        queryOffset = subjectOffset - diagonal 
                            + (PSSMatrix_arr[queryNum].length - wordLengthMinusOne);

                        lastExt = (diagonal << H_BITS) + queryOffset;

                        //fprintf(stderr, "%d %d %d\n", sequenceCount, subjectOffset, queryOffset);
                    }
                    else
                    {
                        lastExt = currHit;
                    }

                    //numUngappedExtSeq++;
                    //extFlag = 1;
                    numTriggerExt++;
                    // blast_numTriggerExtensions_multi[queryNum]++;
                    if (alignment == NULL) {
                        uint4 descriptionStart =
                            sequenceData[sequenceCount].descriptionStart;
                        uint4 encodedLength = sequenceData[sequenceCount].encodedLength;
                        uint4 descriptionLength =
                            sequenceData[sequenceCount].descriptionLength;
                        int4 subjectLength = sequenceData[sequenceCount].sequenceLength;

                        alignment = goodAlignBuf + numGoodAlign;
                        alignment->descriptionLocation = descriptionStart;
                        alignment->descriptionLength = descriptionLength;
                        alignment->subject = subject;
                        alignment->subjectLength = subjectLength;
                        alignment->encodedLength = encodedLength;
                        alignment->joinChecked = 0;
                        alignment->inMemorySubject = 0;
                        alignment->numUnpackRegions = 0;
                        alignment->cluster = 0;
                        alignment->sequenceCount = sequenceCount;
                        alignment->ungappedExtensions = goodExtensionBuf + numGoodExtensions;
                        alignment->gappedExtensions = NULL;
                        alignment->unpackRegions = NULL;
                        alignment->edits = NULL;
                        alignment->numExtensions = 0;
                        alignment->queryCount = queryNum;
                        numTriggerSeq++;
                    }

                    alignment->numExtensions++;
                }
            }
        }

#if 0
        if(numUngappedExtSeq > 0)
        {
            numGoodAlign++;
            numGoodExtensions += numUngappedExtSeq;
        }
#else
        if(numUngappedExtSeq > 0)
        {
            if(alignments_findGoodAlignments_ncbi_multi3(
                        alignment,
                        PSSMatrix_arr[queryNum], scoreMatrix,
                        queryNum, ungappedExtension_new, dp_mem,
                        tree, private_tree, BlastHSP))
            {
                numGoodAlign++;
                numGoodExtensions += numUngappedExtSeq;
            }
        }
#endif
    }

#ifdef PROFILE
    RDTSC_END;
    GET_RDTSC_ATOMIC(&profile->blast_ungappedExtCycle);
#endif
#endif


    (*goodAlignCount) = numGoodAlign;
    (*goodExtensionCount) = numGoodExtensions;

    (*blast_numHits_multi_t) += numHits;
    (*blast_numUngappedExtensions_multi_t) += numUngappedExt;
    (*blast_numTriggerExtensions_multi_t) += numTriggerExt;
    (*blast_numTriggerSequences_multi_t) += numTriggerSeq;

#ifndef NEIGHBOR_INDEX
    free(neighbours);
#endif

}
