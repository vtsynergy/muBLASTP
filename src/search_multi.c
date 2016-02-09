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

#define _GNU_SOURCE
#include "blast.h"
#include <pthread.h>

int numQueryProcessedDetection;

typedef struct {
    int tid;
    struct PSSMatrix *PSSMatrix_arr;
    int numThreads;
    uint4 numSequences;
    struct sequenceData *sequenceData;
    struct scoreMatrix *scoreMatrix;
    int numQuery;
} thread_aux_search_protein2hit;


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
        ) {
    uint4 sequenceCount = 0;
    uint4 descriptionStart = 0, descriptionLength = 0, encodedLength;
    unsigned char *subject, *sequenceEnd;
    int4 subjectLength, subjectOffset, wordLengthMinusOne, count = 0;
    unsigned char currentWord, *currentBlock, *address;
    struct group *currentGroup;
    uint2 *queryOffsets, queryOffset;
    struct ungappedExtension *ungappedExtension;
    int4 diagonal, distance;
    unsigned char **lastHit;

    wordLengthMinusOne = parameters_wordSize - 1;

    //hit_t *firstBin = (hit_t *)malloc(sizeof(hit_t) * MAX_NUM_HIT_PER_QUERY * numQuery);
    //hit_t *secondBin = (hit_t *)malloc(sizeof(hit_t) * MAX_NUM_HIT_PER_QUERY * numQuery);

    //fprintf(stderr, "bin size: %d MB\n", (sizeof(hit_t) * MAX_NUM_HIT_PER_QUERY * numQuery)>>20);
    int4 numGoodExtensions = *goodExtensionCount;
    int4 numGoodAlign = *goodAlignCount;

    int numSecondBin = numQuery;
    //int4 *binOffset = (int4 *)malloc(sizeof(int4) * maxNumFirstBin);
    //uint4 *numExtHit = (uint4 *)malloc(sizeof(uint4) * numSecondBin);
    //uint4 *numHit = (uint4 *)malloc(sizeof(uint4) * numSecondBin);
    //hit_t *lastHits = (hit_t *)malloc(sizeof(hit_t) * numSecondBin);

    sequenceCount = 0;
    //int seqPrecent = numSequences/100;
    while (sequenceCount < numSequences) {

        //if(sequenceCount%seqPrecent == 0)
        //{
            //fprintf(stderr, "finished %d\%\n", sequenceCount/seqPrecent);
        //}
        descriptionLength = sequenceData[sequenceCount].descriptionLength;
        descriptionStart = sequenceData[sequenceCount].descriptionStart;
        subjectLength = sequenceData[sequenceCount].sequenceLength;
        //encodedLength = sequenceData[sequenceCount].encodedLength;
        address = subject = sequenceData[sequenceCount].sequence;

        int numFirstBin = subjectLength + longestQueryLength + 1;
        // New sequence, new possible alignment
        //alignments_currentAlignment_multi[queryNum] = NULL;

        // Only process sequence if at least as long as the word length
        if (subjectLength >= parameters_wordSize) {

            memset(binOffset, 0, sizeof(int4) * numFirstBin);
            memset(numHit, 0, sizeof(uint4) * numQuery);

            sequenceEnd = subject + subjectLength - wordLengthMinusOne;
            int4 hitCount = 0;
            while (address < sequenceEnd) {
                int codeWord = wordLookupDFA_getCodeword(address, parameters_wordSize);
                if (PV_TEST(pv, codeWord, PV_ARRAY_BTS)) 
                {
                    subjectOffset = address - subject;

                    int ii;
                    for(ii = 0; ii < proteinLookup_query[codeWord].numQueryPositions; ii++) {
                        queryOffset = proteinLookup_query[codeWord].querySequencePositions[ii] >> 16;
                        int queryNum = proteinLookup_query[codeWord].querySequencePositions[ii] & 0xffff;
                        diagonal =
                            subjectOffset - queryOffset + (PSSMatrix_arr[queryNum].length - wordLengthMinusOne);
                        secondBin[hitCount] = proteinLookup_query[codeWord].querySequencePositions[ii]; 
                        binOffset[diagonal + 1]++;
                        numHit[queryNum]++;
                        hitCount++;
                        //blast_numHits_multi[queryNum]++;
                    }
                }
                address++;
            }

            if(hitCount >= MAX_NUM_HIT_PER_QUERY * numQuery)
            {
                fprintf(stderr, "hitCount > MAX_NUM_HIT_PER_QUERY * numQuery\n");
                exit(1);
            }

            int jj;
            for (jj = 1; jj < numFirstBin; jj++) {
                binOffset[jj] = binOffset[jj - 1] + binOffset[jj];
            }

            int ii;
            address = sequenceData[sequenceCount].sequence;
            hitCount = 0;
            while (address < sequenceEnd) {
                subjectOffset = address - subject;
                int codeWord = wordLookupDFA_getCodeword(address, parameters_wordSize);
                if(PV_TEST(pv, codeWord, PV_ARRAY_BTS)) 
                {
                    for(ii = 0; ii < proteinLookup_query[codeWord].numQueryPositions; ii++) {
                        queryOffset = secondBin[hitCount] >> 16; 
                        int queryNum = secondBin[hitCount] & 0xffff;
                        diagonal =
                            subjectOffset - queryOffset + (PSSMatrix_arr[queryNum].length - wordLengthMinusOne);
                        //secondBin[hitCount] = proteinLookup_query[codeWord].querySequencePositions[ii]; 
                        firstBin[binOffset[diagonal]] = secondBin[hitCount];
                        binOffset[diagonal]++;
                        hitCount++;
                        //blast_numHits_multi[queryNum]++;
                    }}
                address++;
            }

            memset(numExtHit, 0, sizeof(uint4) * numSecondBin);
            memset(lastHits, 255, sizeof(hit_t) * numSecondBin);

            int kk = 0;
            for (jj = 0; jj < numFirstBin - 1; jj++)
            {
                for(; kk < binOffset[jj]; kk++)
                {
                    uint2 diagNum = jj;
                    hit_t hit = firstBin[kk];
                    uint2 queryNum = hit & 0xffff;
                    hit_t lastHit = lastHits[queryNum];
                    hit_t currHit = (hit >> 16) + (diagNum << 16);
                    int distance = currHit - lastHit;
                    if (distance >= parameters_A || lastHit == 0xffffffff) {
                        lastHits[queryNum] = currHit;
                    } else if (distance >= parameters_overlap) {
                        numExtHit[queryNum] += 2;
                        if (numExtHit[queryNum] > MAX_NUM_HIT_PER_QUERY) {
                            fprintf(stderr, "numExtHit > MAX_NUM_HIT_PER_QUERY\n");
                        }

                        secondBin[queryNum * MAX_NUM_HIT_PER_QUERY + numExtHit[queryNum] - 2] = lastHit;
                        secondBin[queryNum * MAX_NUM_HIT_PER_QUERY + numExtHit[queryNum] - 1] = currHit;
                        lastHits[queryNum] = currHit;
                    }
                }
            }

            for(jj = 0; jj < numQuery; jj++)
            {
                blast_numHits_multi[jj] += numHit[jj];
            }

            for(jj = 0; jj < numSecondBin; jj++)
            {
                if(numExtHit[jj] < 2)
                    continue;

                struct alignment *alignment = NULL;

                hit_t lastExt = 0;

                unsigned char *ungappedExtension_subjectEndReached_t;

                int numUngappedExtSeq = 0;
                int queryNum = jj;
                for (kk = 0; kk < numExtHit[jj]; kk += 2){
                    hit_t lastHit = secondBin[jj * MAX_NUM_HIT_PER_QUERY + kk];
                    hit_t currHit = secondBin[jj * MAX_NUM_HIT_PER_QUERY + kk + 1];
                    int distance = currHit - lastHit;

                    if(currHit > lastExt)
                    {
                        blast_numUngappedExtensions_multi[queryNum]++;
                        hit_t hit = currHit;
                        int diagonal = hit >> 16;
                        queryOffset = hit & 0xffff;

                        subjectOffset =
                            queryOffset  + (diagonal - (PSSMatrix_arr[queryNum].length - wordLengthMinusOne));

                        subjectOffset += 2;
                        queryOffset += 2;

                        //fprintf(stderr, "%d %d %d\n", sequenceCount, subjectOffset, queryOffset);

                        unsigned char *address = subject + subjectOffset;
                        unsigned char *lastHit_addr = address - distance;

                        ASSERT(lastHit_addr < address && lastHit_addr >= subject);

                        char rightExtend = 0;
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

                        if (ungappedExtension) {

                            if(rightExtend)
                            {
                                subjectOffset = (ungappedExtension_subjectEndReached_t -
                                        subject - wordLengthMinusOne);
                                queryOffset = subjectOffset - diagonal + (PSSMatrix_arr[queryNum].length - wordLengthMinusOne);
                                lastExt =
                                    (currHit & 0xffff0000) + queryOffset;

                                //fprintf(stderr, "%d %d %d\n", sequenceCount, subjectOffset, queryOffset);

                            }
                            else
                            {
                                lastExt = currHit;
                            }

                            blast_numTriggerExtensions_multi[queryNum]++;
                            if (alignment == NULL) {
                                //uint4 descriptionStart =
                                //sequenceData[sequenceCount].descriptionStart;
                                //uint4 encodedLength = sequenceData[sequenceCount].encodedLength;
                                //uint4 descriptionLength =
                                //sequenceData[sequenceCount].descriptionLength;
                                //int4 subjectLength = sequenceData[sequenceCount].sequenceLength;

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
                                blast_numTriggerSequences_multi[queryNum]++;
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
                    }
                }

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
            }

        }
        sequenceCount++;
    }

    (*goodAlignCount) = numGoodAlign;
    (*goodExtensionCount) = numGoodExtensions;
    //free(numExtHit);
    //free(numHit);
    //free(lastHits);
    //free(binOffset);
    //free(firstBin);
    //free(secondBin);
}

// Search a protein database using 2-hit extension mode
void search_protein2hit_multi(struct PSSMatrix *PSSMatrix_arr,
        struct sequenceData *sequenceData,
        uint4 numSequences, int queryNum) {
    uint4 sequenceCount = 0;
    uint4 descriptionStart = 0, descriptionLength = 0, encodedLength;
    unsigned char *subject, *sequenceEnd;
    int4 subjectLength, subjectOffset, wordLengthMinusOne, count = 0;
    unsigned char currentWord, *currentBlock, *address;
    struct group *currentGroup;
    uint2 *queryOffsets, queryOffset;
    struct ungappedExtension *ungappedExtension;
    int4 diagonal, distance;
    unsigned char **lastHit;

    wordLengthMinusOne = parameters_wordSize - 1;

    sequenceCount = 0;
    while (sequenceCount < numSequences) {
        descriptionLength = sequenceData[sequenceCount].descriptionLength;
        descriptionStart = sequenceData[sequenceCount].descriptionStart;
        subjectLength = sequenceData[sequenceCount].sequenceLength;
        encodedLength = sequenceData[sequenceCount].encodedLength;
        address = subject = sequenceData[sequenceCount].sequence;

        // New sequence, new possible alignment
        alignments_currentAlignment_multi[queryNum] = NULL;

        // Only process sequence if at least as long as the word length
        if (subjectLength >= parameters_wordSize) {
            // Start at 000 state in Finite State Automata
            currentGroup = wordLookupDFA_groups_multi[queryNum];

            // Read first wordLength - 1 chars and advance
            count = wordLengthMinusOne;
            while (count > 0) {
                if (*address < wordLookupDFA_numCodes)
                    currentGroup = currentGroup->nextGroups + *address;
                else
                    currentGroup = currentGroup->nextGroups;

                address++;
                count--;
            }

            // Read the rest of the codes, using the Finite State Automata defined
            // by wordLookupDFA to get the query positions of int4erest
            sequenceEnd = subject + subjectLength;
            while (address < sequenceEnd) {
                currentBlock = currentGroup->nextWords;

                // If current code is a regular letter
                if (*address < wordLookupDFA_numCodes) {
                    // Use it
                    currentWord = currentBlock[*address];
                    currentGroup = currentGroup->nextGroups + *address;
                } else {
                    // Else check if we've reached end of the file
                    if (address >= sequenceEnd)
                        break;

                    // If not, we've read a wild code. Use first code instead
                    currentWord = currentBlock[0];
                    currentGroup = currentGroup->nextGroups;
                }

                if (currentWord) {
                    // Calculate subject offset
                    subjectOffset = address - subject;

                    // If at least one query position, stored at an extenal address
                    queryOffsets = ((uint2 *)currentBlock) - currentWord;

                    // If the zero flag is stored at the first query position
                    if (!*queryOffsets) {
                        // Go to an outside address for additional positions
                        queryOffsets =
                            wordLookupDFA_additionalQueryPositions_multi[queryNum] +
                            (*(queryOffsets + 1) * constants_max_int2) +
                            *(queryOffsets + 2);
                    }

                    do {
                        queryOffset = *queryOffsets;

#ifndef NO_STAGE2
                        // Calculate the diagonal this hit is on
                        diagonal = subjectOffset - queryOffset;

                        // Calculate distance since last hit
                        lastHit = hitMatrix_furthest_multi[queryNum] + diagonal;
                        distance = address - *lastHit;

                        if (distance >= parameters_A) {
                            // Too far apart, update furthest
                            *lastHit = address;
                        } else if (distance >= parameters_overlap) {

                            // Not overlaping - extension triggered
                            // Increment tally number of extensions
                            blast_numUngappedExtensions_multi[queryNum]++;

                            // Perform ungapped extension start between query/subject
                            // start/end
                            // and extending outwards in each direction
#if 1
                            ungappedExtension = ungappedExtension_extend_multi(
                                    PSSMatrix_arr[queryNum].matrix + queryOffset, address,
                                    *lastHit, PSSMatrix_arr[queryNum], subject, queryNum,
                                    &ungappedExtension_subjectEndReached);

                            // Update furthest reached value for the diagonal
                            *lastHit = ungappedExtension_subjectEndReached;

                            //// If extension scores high enough to trigger gapping
                            if (ungappedExtension) {
                                // Increment count of number of trigger extensions
                                blast_numTriggerExtensions_multi[queryNum]++;

                                // Create new alignment if needed
                                if (alignments_currentAlignment_multi[queryNum] == NULL) {
                                    // Create new alignment object using subject with wilds
                                    alignments_createNew_multi(
                                            descriptionStart, descriptionLength, subject,
                                            subjectLength, encodedLength, queryNum);
                                }

                                // Add this extension to the alignment
                                alignments_addUngappedExtension_multi(ungappedExtension,
                                        queryNum);
                            }
#endif
                        }
#endif

                        queryOffsets++;
                        blast_numHits_multi[queryNum]++;
                    } while (*queryOffsets);
                }

                address++;
            }
        }
        sequenceCount++;
    }
    // printf("%d %d %d\n", queryNum, blast_numHits_multi[queryNum],
    // blast_numUngappedExtensions_multi[queryNum]);
}

// Search a protein database using 2-hit extension mode
void search_protein2hit_core(int tid, struct PSSMatrix *PSSMatrix_arr,
        struct sequenceData *sequenceData,
        uint4 numSequences, int numQuery) {
    uint4 sequenceCount = 0;
    uint4 descriptionStart = 0, descriptionLength = 0, encodedLength;
    unsigned char *subject, *sequenceEnd;
    int4 subjectLength, subjectOffset, wordLengthMinusOne, count = 0;
    unsigned char currentWord, *currentBlock, *address;
    struct group *currentGroup;
    uint2 *queryOffsets, queryOffset;
    struct ungappedExtension *ungappedExtension;
    int4 diagonal, distance;
    unsigned char **lastHit;

    wordLengthMinusOne = parameters_wordSize - 1;

    int queryNum;
    while ((queryNum = __sync_fetch_and_add(&numQueryProcessedDetection, 1)) <
            numQuery) {
        sequenceCount = 0;
        while (sequenceCount < numSequences) {
            descriptionLength = sequenceData[sequenceCount].descriptionLength;
            descriptionStart = sequenceData[sequenceCount].descriptionStart;
            subjectLength = sequenceData[sequenceCount].sequenceLength;
            encodedLength = sequenceData[sequenceCount].encodedLength;
            address = subject = sequenceData[sequenceCount].sequence;

            // New sequence, new possible alignment
            alignments_currentAlignment_multi[queryNum] = NULL;

            // Only process sequence if at least as long as the word length
            if (subjectLength >= parameters_wordSize) {
                // Start at 000 state in Finite State Automata
                currentGroup = wordLookupDFA_groups_multi[queryNum];

                // Read first wordLength - 1 chars and advance
                count = wordLengthMinusOne;
                while (count > 0) {
                    if (*address < wordLookupDFA_numCodes)
                        currentGroup = currentGroup->nextGroups + *address;
                    else
                        currentGroup = currentGroup->nextGroups;

                    address++;
                    count--;
                }

                // Read the rest of the codes, using the Finite State Automata defined
                // by wordLookupDFA to get the query positions of int4erest
                sequenceEnd = subject + subjectLength;
                while (address < sequenceEnd) {
                    currentBlock = currentGroup->nextWords;

                    // If current code is a regular letter
                    if (*address < wordLookupDFA_numCodes) {
                        // Use it
                        currentWord = currentBlock[*address];
                        currentGroup = currentGroup->nextGroups + *address;
                    } else {
                        // Else check if we've reached end of the file
                        if (address >= sequenceEnd)
                            break;

                        // If not, we've read a wild code. Use first code instead
                        currentWord = currentBlock[0];
                        currentGroup = currentGroup->nextGroups;
                    }

                    if (currentWord) {
                        // Calculate subject offset
                        subjectOffset = address - subject;

                        // If at least one query position, stored at an extenal address
                        queryOffsets = ((uint2 *)currentBlock) - currentWord;

                        // If the zero flag is stored at the first query position
                        if (!*queryOffsets) {
                            // Go to an outside address for additional positions
                            queryOffsets =
                                wordLookupDFA_additionalQueryPositions_multi[queryNum] +
                                (*(queryOffsets + 1) * constants_max_int2) +
                                *(queryOffsets + 2);
                        }

                        do {
                            queryOffset = *queryOffsets;

#ifndef NO_STAGE2
                            // Calculate the diagonal this hit is on
                            diagonal = subjectOffset - queryOffset;

                            // Calculate distance since last hit
                            lastHit = hitMatrix_furthest_multi[queryNum] + diagonal;
                            distance = address - *lastHit;

                            if (distance >= parameters_A) {
                                // Too far apart, update furthest
                                *lastHit = address;
                            } else if (distance >= parameters_overlap) {

                                // Not overlaping - extension triggered
                                // Increment tally number of extensions
                                blast_numUngappedExtensions_multi[queryNum]++;

                                // Perform ungapped extension start between query/subject
                                // start/end
                                // and extending outwards in each direction
                                ungappedExtension = ungappedExtension_extend_multi(
                                        PSSMatrix_arr[queryNum].matrix + queryOffset, address,
                                        *lastHit, PSSMatrix_arr[queryNum], subject, queryNum,
                                        &ungappedExtension_subjectEndReached);

                                // Update furthest reached value for the diagonal
                                *lastHit = ungappedExtension_subjectEndReached;

                                //// If extension scores high enough to trigger gapping
                                if (ungappedExtension) {
                                    // Increment count of number of trigger extensions
                                    blast_numTriggerExtensions_multi[queryNum]++;

                                    // Create new alignment if needed
                                    if (alignments_currentAlignment_multi[queryNum] == NULL) {
                                        // Create new alignment object using subject with wilds
                                        alignments_createNew_multi(
                                                descriptionStart, descriptionLength, subject,
                                                subjectLength, encodedLength, queryNum);
                                    }

                                    // Add this extension to the alignment
                                    alignments_addUngappedExtension_multi(ungappedExtension,
                                            queryNum);
                                }
                            }
#endif

                            queryOffsets++;
                            blast_numHits_multi[queryNum]++;
                        } while (*queryOffsets);
                    }

                    address++;
                }
            }
            sequenceCount++;
        }
        // printf("%d %d %d\n", queryNum, blast_numHits_multi[queryNum],
        // blast_numUngappedExtensions_multi[queryNum]);
    }
}

static void *worker_dfa_query(void *data) {
    thread_aux_search_protein2hit *d = (thread_aux_search_protein2hit *)data;

    search_protein2hit_core(d->tid, d->PSSMatrix_arr, d->sequenceData,
            d->numSequences, d->numQuery);
    return 0;
}

// Search a protein database using 2-hit extension mode
void search_protein2hit_multithread(struct PSSMatrix *PSSMatrix_arr,
        struct sequenceData *sequenceData,
        uint4 numSequences, int numQuery) {
    uint4 queryCount = 0;
    int ii = 0, jj = 0;
    int n_threads = parameters_num_threads;
    thread_aux_search_protein2hit *data = (thread_aux_search_protein2hit *)malloc(
            sizeof(thread_aux_search_protein2hit) * n_threads);
    pthread_t *tid = (pthread_t *)malloc(sizeof(pthread_t) * n_threads);
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    int cpu_map[MAX_NUM_THREADS];
    get_affinity(cpu_map);

    numQueryProcessedDetection = 0;

    for (ii = 0; ii < n_threads; ii++) {
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(cpu_map[ii], &cpuset);
        data[ii].tid = ii;
        data[ii].PSSMatrix_arr = PSSMatrix_arr;
        data[ii].numSequences = numSequences;
        data[ii].sequenceData = sequenceData;
        data[ii].numQuery = numQuery;
        int err;
        pthread_create(&tid[ii], &attr, worker_dfa_query, data + ii);
        err = pthread_setaffinity_np(tid[ii], sizeof(cpu_set_t), &cpuset);
        if (err != 0)
            fprintf(stderr, "pthread_setaffinity_np error: %d\n", err);
        // printf("thread %d mapped to core %d\n", ii, cpu_map[ii]);
    }

    for (ii = 0; ii < n_threads; ++ii) {
        pthread_join(tid[ii], 0);
    }

    // printf("total blast_numHits: %d\n", blast_numHits);
    free(data);
}
