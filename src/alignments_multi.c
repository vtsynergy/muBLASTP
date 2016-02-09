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
#include <omp.h>
#ifndef INT2_MIN
/** smallest (most negative) number represented by signed (two byte) short */
#define INT2_MIN    (-32768)
#endif

#include<sys/time.h>

#ifndef INT2_MAX
/** smallest (most negative) number represented by signed (two byte) short */
#define INT2_MAX    (32767)
#endif



struct memBlocks *alignments_alignments_multi[BATCH_SIZE];
struct memBlocks ***alignments_alignments_multi2;
struct memSingleBlock *alignments_goodAlignments_multi[BATCH_SIZE];
struct memSingleBlock *alignments_finalAlignments_multi[BATCH_SIZE];
struct alignment *alignments_currentAlignment_multi[BATCH_SIZE];

struct memBlocks *alignments_volumeAlignments_multi[BATCH_SIZE]
[constants_maxNumVolumes];
struct memBlocks ***alignments_volumeAlignments_multi2;

uint4 alignments_numVolumes_multi[BATCH_SIZE];
int alignments_finalAlignmentsSorted_multi[BATCH_SIZE];

void alignments_addGoodAlignment_multi(int4 highestNominalScore,
        struct alignment *alignment,
        int queryNum);
void alignments_sortFinalAlignments_multi(int queryNum);
void alignments_checkForJoin_multi(struct alignment *alignment,
        struct ungappedExtension *extension1,
        struct PSSMatrix PSSMatrix, int queryNum);

void alignments_regularGappedAlignment_multi(
        struct PSSMatrix PSSMatrix, struct ungappedExtension *ungappedExtension,
        struct alignment *alignment, int queryNum);

int alignments_expandCluster_multi(struct alignment *alignment,
        struct PSSMatrix PSSMatrix, int queryNum);

// Initialize array storing pointers to alignments
void alignments_initialize_multi2() {
    int ii;
    // Initialize alignments, good alignments and final alignments blocks
    //alignments_alignments_multi2 = (struct memBlocks ***)malloc(
    //sizeof(struct memBlocks **) * blast_numQuery);
    //alignments_volumeAlignments_multi2 = (struct memBlocks ***)malloc(
    //sizeof(struct memBlocks **) * blast_numQuery);
    for (ii = 0; ii < blast_numQuery; ii++) {
        //alignments_alignments_multi2[ii] = (struct memBlocks **)malloc(
        //sizeof(struct memBlocks *) * blast_numBlocks);
        //alignments_volumeAlignments_multi2[ii] = (struct memBlocks **)malloc(
        //sizeof(struct memBlocks *) * blast_numBlocks * constants_maxNumVolumes);
        //int jj;
        //for (jj = 0; jj < blast_numBlocks; jj++) {
        //alignments_alignments_multi2[ii][jj] = memBlocks_initialize(
        //sizeof(struct alignment),
        //constants_initialAllocAlignments);
        //}
        //alignments_alignments_multi[ii] = memBlocks_initialize(
        //sizeof(struct alignment),
        //constants_initialAllocAlignments);
        alignments_goodAlignments_multi[ii] = memSingleBlock_initialize(
                sizeof(struct finalAlignment), constants_initialAllocGoodAlignments);
        alignments_finalAlignments_multi[ii] = memSingleBlock_initialize(
                sizeof(struct finalAlignment), constants_initialAllocFinalAlignments);
        alignments_currentAlignment_multi[ii] = NULL;
        alignments_numVolumes_multi[ii] = 0;
        alignments_finalAlignmentsSorted_multi[ii] = 0;

        // Initialize code for processing regions of a subject
        unpack_initialize_multi(ii);
    }
}

// Initialize array storing pointers to alignments
void alignments_initialize_multi() {
    int ii;
    // Initialize alignments, good alignments and final alignments blocks
    for (ii = 0; ii < blast_numQuery; ii++) {
        alignments_alignments_multi[ii] = memBlocks_initialize(
                sizeof(struct alignment),
                constants_initialAllocAlignments);
        alignments_goodAlignments_multi[ii] = memSingleBlock_initialize(
                sizeof(struct finalAlignment), constants_initialAllocGoodAlignments);
        alignments_finalAlignments_multi[ii] = memSingleBlock_initialize(
                sizeof(struct finalAlignment), constants_initialAllocFinalAlignments);
        alignments_currentAlignment_multi[ii] = NULL;
        alignments_numVolumes_multi[ii] = 0;
        alignments_finalAlignmentsSorted_multi[ii] = 0;

        // Initialize code for processing regions of a subject
        unpack_initialize_multi(ii);
    }
}

// Use the next available slot in the alignment object array
void alignments_createNew_multi(uint4 descriptionLocation,
        uint4 descriptionLength, unsigned char *subject,
        int4 subjectLength, int4 encodedLength,
        int queryNum) {
    struct alignment *alignment;

    // Get slot for new alignment
    alignment = (struct alignment *)memBlocks_newEntry(
            alignments_alignments_multi[queryNum]);

    // Create/initialize contents
    alignment->ungappedExtensions = NULL;
    alignment->gappedExtensions = NULL;
    alignment->descriptionLocation = descriptionLocation;
    alignment->descriptionLength = descriptionLength;
    alignment->subject = subject;
    alignment->subjectLength = subjectLength;
    alignment->joinChecked = 0;
    alignment->encodedLength = encodedLength;
    alignment->inMemorySubject = 0;
    alignment->unpackRegions = NULL;
    alignment->numUnpackRegions = 0;
    alignment->edits = NULL;
    alignment->cluster = 0;

    // Record pointer to wilcard edits if there are any for this subject
    if (encoding_alphabetType == encoding_nucleotide) {
        alignment->edits = subject + ((alignment->subjectLength + 3) / 4);
        if (alignment->edits == subject + encodedLength) {
            alignment->edits = NULL;
        }
    }

    alignments_currentAlignment_multi[queryNum] = alignment;
    blast_numTriggerSequences_multi[queryNum]++;
}

// Use the next available slot in the alignment object array
struct alignment *alignments_createNew_multi_db2(
        struct alignment *goodAlignBuf, int *goodAlignCount,
        uint4 descriptionLocation, uint4 descriptionLength, unsigned char *subject,
        int4 subjectLength, int4 encodedLength,
        uint4 sequenceCount) {
    struct alignment *alignment = goodAlignBuf + *goodAlignCount;
    //(*goodAlignCount)++;

    //if(*goodAlignCount >= MAX_ALIGNMENTS_PER_QUERY)
    //{
    //fprintf(stderr, "goodAlignCount larger than MAX_ALIGNMENTS_PER_QUERY\n");
    //exit(1);
    //}

    // Get slot for new alignment
    //alignment = (struct alignment *)memBlocks_newEntry(
    //alignments_alignments_multi2[queryNum][blockNum]);

    // Create/initialize contents
    alignment->ungappedExtensions = NULL;
    alignment->gappedExtensions = NULL;
    alignment->descriptionLocation = descriptionLocation;
    alignment->descriptionLength = descriptionLength;
    alignment->subject = subject;
    alignment->subjectLength = subjectLength;
    alignment->joinChecked = 0;
    alignment->encodedLength = encodedLength;
    alignment->inMemorySubject = 0;
    alignment->unpackRegions = NULL;
    alignment->numUnpackRegions = 0;
    alignment->edits = NULL;
    alignment->cluster = 0;
    alignment->sequenceCount = sequenceCount;

    // Record pointer to wilcard edits if there are any for this subject
    if (encoding_alphabetType == encoding_nucleotide) {
        alignment->edits = subject + ((alignment->subjectLength + 3) / 4);
        if (alignment->edits == subject + encodedLength) {
            alignment->edits = NULL;
        }
    }

    return alignment;
}

// Use the next available slot in the alignment object array
struct alignment *alignments_createNew_multi_db(
        uint4 descriptionLocation, uint4 descriptionLength, unsigned char *subject,
        int4 subjectLength, int4 encodedLength, int queryNum, uint4 sequenceCount) {
    struct alignment *alignment;

    // Get slot for new alignment
    alignment = (struct alignment *)memBlocks_newEntry(
            alignments_alignments_multi[queryNum]);

    // Create/initialize contents
    alignment->ungappedExtensions = NULL;
    alignment->gappedExtensions = NULL;
    alignment->descriptionLocation = descriptionLocation;
    alignment->descriptionLength = descriptionLength;
    alignment->subject = subject;
    alignment->subjectLength = subjectLength;
    alignment->joinChecked = 0;
    alignment->encodedLength = encodedLength;
    alignment->inMemorySubject = 0;
    alignment->unpackRegions = NULL;
    alignment->numUnpackRegions = 0;
    alignment->edits = NULL;
    alignment->cluster = 0;
    alignment->sequenceCount = sequenceCount;

    // Record pointer to wilcard edits if there are any for this subject
    if (encoding_alphabetType == encoding_nucleotide) {
        alignment->edits = subject + ((alignment->subjectLength + 3) / 4);
        if (alignment->edits == subject + encodedLength) {
            alignment->edits = NULL;
        }
    }

    return alignment;
}

// Add an ungapped extension to an alignment's list of ungapped extensions,
// which
// are ordered highest nominal score first, lowest score last
void
alignments_addUngappedExtension_multi(struct ungappedExtension *newExtension,
        int queryNum) {
    struct ungappedExtension *currentExtension, *previousExtension;

    // Start at beginning of ungapped extensions (one with highest score)
    currentExtension =
        alignments_currentAlignment_multi[queryNum]->ungappedExtensions;

    // If there are none, add first/sole extension
    if (currentExtension == NULL) {
        alignments_currentAlignment_multi[queryNum]->ungappedExtensions =
            newExtension;
    } else {
        previousExtension = NULL;

        // Else move through list of existing extensions until we either reach the
        // end or reach one with a score less than newExtension
        while (currentExtension != NULL &&
                (currentExtension->nominalScore > newExtension->nominalScore)) {
            previousExtension = currentExtension;
            currentExtension = currentExtension->next;
        }

        if (previousExtension == NULL) {
            // This is the highest scoring alignment, insert at front of
            // the queue
            alignments_currentAlignment_multi[queryNum]->ungappedExtensions =
                newExtension;
            newExtension->next = currentExtension;
        } else {
            // Insert between higher and lower scoring extensions
            previousExtension->next = newExtension;
            newExtension->next = currentExtension;
        }
    }
}

// Add an ungapped extension to an alignment's list of ungapped extensions,
// which
// are ordered highest nominal score first, lowest score last
void alignments_addUngappedExtension_multi_db(
        struct ungappedExtension *newExtension, int queryNum,
        struct alignment *alignments_currentAlignment_multi_t) {
    struct ungappedExtension *currentExtension, *previousExtension;

    // Start at beginning of ungapped extensions (one with highest score)
    currentExtension = alignments_currentAlignment_multi_t->ungappedExtensions;

    // If there are none, add first/sole extension
    if (currentExtension == NULL) {
        alignments_currentAlignment_multi_t->ungappedExtensions = newExtension;
    } else {
        previousExtension = NULL;

        // Else move through list of existing extensions until we either reach the
        // end or reach one with a score less than newExtension
        while (currentExtension != NULL &&
                (currentExtension->nominalScore > newExtension->nominalScore)) {
            previousExtension = currentExtension;
            currentExtension = currentExtension->next;
        }

        if (previousExtension == NULL) {
            // This is the highest scoring alignment, insert at front of
            // the queue
            alignments_currentAlignment_multi_t->ungappedExtensions = newExtension;
            newExtension->next = currentExtension;
        } else {
            // Insert between higher and lower scoring extensions
            previousExtension->next = newExtension;
            newExtension->next = currentExtension;
        }
    }
}


// Perform initial scoring of all ungapped extensions to find "good" alignments
// that may
// score above the cutoff
int alignments_findGoodAlignments_multi3(struct alignment *alignment,
        struct PSSMatrix PSSMatrix,
        int queryNum) {
    struct ungappedExtension *ungappedExtension;
    int4 bestScore, nExt, hasChildren;
    uint4 numExtensions = alignment->numExtensions;

    bestScore = 0;
    blast_dloc_multi[queryNum] = alignment->descriptionLocation;
    // Record if subject has children
    if (encoding_alphabetType == encoding_protein &&
            alignment->encodedLength > alignment->subjectLength + 2)
    {
        hasChildren = 1;
        fprintf(stderr, "hasChildren: %d\n", hasChildren);
    }
    else
        hasChildren = 0;

    // Sort the ungappedExtensions
    int ii, jj;
    for(ii = 0; ii < (numExtensions - 1); ii++)
    {
        int iMax = ii;

        for(jj = ii + 1; jj < numExtensions; jj++)
        {
            struct ungappedExtension *curUG = alignment->ungappedExtensions + jj;
            struct ungappedExtension *maxUG = alignment->ungappedExtensions + iMax;

            if(curUG->nominalScore >= maxUG->nominalScore)
                iMax = jj;
        }

        if(iMax != ii)
        {
            struct ungappedExtension temp = alignment->ungappedExtensions[iMax];
            for(jj = iMax; jj > ii; jj--)
            {
                alignment->ungappedExtensions[jj] = alignment->ungappedExtensions[jj - 1];
            }
            alignment->ungappedExtensions[ii] = temp;
        }
    }


    for(ii = 0; ii < numExtensions; ii++)
    {
        ungappedExtension = alignment->ungappedExtensions + ii;
        ungappedExtension->next = alignment->ungappedExtensions + ii + 1;
    }
    ungappedExtension->next = NULL;
    nExt = 0;
    // For each ungapped extension (in descending order of score)
    for(ii = 0; ii < numExtensions; ii++) {
        ungappedExtension = alignment->ungappedExtensions + ii;

        if (ungappedExtension->status != ungappedExtension_DELETED) {
            // Find the seed
            ungappedExtension_findSeed(ungappedExtension, PSSMatrix,
                    alignment->subject);

            // Semi-gapped scoring
            if (parameters_semiGappedScoring) {
                ungappedExtension->nominalScore = semiGappedScoring_score_multi(
                        ungappedExtension, PSSMatrix, alignment->subjectLength,
                        alignment->subject,
                        statistics_gappedNominalDropoff_multi[queryNum], queryNum);
            }
            // Regular gapped scoring
            else {
                ungappedExtension->nominalScore = gappedScoring_score_multi(
                        ungappedExtension, PSSMatrix, alignment->subjectLength,
                        alignment->subject,
                        statistics_gappedNominalDropoff_multi[queryNum], queryNum);

                // Mark as gapped
                ungappedExtension->status = ungappedExtension_GAPPED;
            }

            // If alignment scores above R1 cutoff
            if (ungappedExtension->nominalScore >=
                    blast_nominalR1cutoff_multi[queryNum]) {
                if (hasChildren) {
                    // Subject has children so perform stage1 and 2 on children
                    alignments_expandCluster_multi(alignment, PSSMatrix, queryNum);
                    bestScore = 0;
                    break;
                } else if (ungappedExtension->nominalScore > bestScore) {
                    // Update best score for the alignment
                    bestScore = ungappedExtension->nominalScore;
                }
            } else {
                // Else mark it as deleted
                ungappedExtension->status = ungappedExtension_DELETED;
            }


            // Remove any ungapped extensions in this alignment that are in the
            // area covered
            // by the gapped scoring just performed
            alignments_pruneRegion_multi(alignment, ungappedExtension, queryNum);
            nExt++;
        }
    }
    // If this alignment contains gapped extensions that could score above
    // cutoff
    if (bestScore >= blast_nominalR1cutoff_multi[queryNum]) {
        // If a single sequence add to list of "good" alignments
        alignments_addGoodAlignment_multi(bestScore, alignment, queryNum);

        blast_numGoodExtensions_multi[queryNum] += nExt;
        blast_numGoodAlignments_multi[queryNum]++;
        return 1;
    }
    else {
        return 0;
    }
}


// Perform initial scoring of all ungapped extensions to find "good" alignments
// that may
// score above the cutoff
void alignments_findGoodAlignments_multi2(struct PSSMatrix PSSMatrix, struct scoreMatrix scoreMatrix,
        int queryNum) {
    struct alignment *alignment;
    struct ungappedExtension *ungappedExtension;
    int4 bestScore, numExtensions, hasChildren;
    //printf("alignments_findGoodAlignments_multi2\n");
    // For each alignment
    int blockNum;
    for (blockNum = 0; blockNum < blast_numBlocks; blockNum++) {
        memBlocks_resetCurrent(alignments_alignments_multi2[queryNum][blockNum]);
        while ((alignment = (struct alignment *)memBlocks_getCurrent(
                        alignments_alignments_multi2[queryNum][blockNum])) != NULL) {
            bestScore = 0;
            blast_dloc_multi[queryNum] = alignment->descriptionLocation;

            // Record if subject has children
            if (encoding_alphabetType == encoding_protein &&
                    alignment->encodedLength > alignment->subjectLength + 2)
            {
                hasChildren = 1;
                fprintf(stderr, "hasChildren: %d\n", hasChildren);
            }
            else
                hasChildren = 0;


            // For each ungapped extension (in descending order of score)
            numExtensions = 0;

            ungappedExtension = alignment->ungappedExtensions;

            while (ungappedExtension != NULL) {
                if (ungappedExtension->status != ungappedExtension_DELETED) {
                    // Find the seed
                    ungappedExtension_findSeed(ungappedExtension, PSSMatrix,
                            alignment->subject);

                    // Byte-packed scoring
                    if (parameters_bytepackedScoring) {
                        ungappedExtension->nominalScore = bytepackGappedScoring_score(
                                ungappedExtension, PSSMatrix, alignment->subjectLength,
                                alignment->subject,
                                statistics_gappedNominalDropoff_multi[queryNum]);

                        // Mark as semigapped
                        ungappedExtension->status = ungappedExtension_SEMIGAPPED;
                    }
                    // Table driven scoring
                    else if (parameters_tableScoring) {
                        ungappedExtension->nominalScore = tableGappedScoring_score(
                                ungappedExtension, PSSMatrix, alignment->subjectLength,
                                alignment->subject,
                                statistics_gappedNominalDropoff_multi[queryNum]);

                        // Mark as semigapped
                        ungappedExtension->status = ungappedExtension_SEMIGAPPED;
                    }
                    // Semi-gapped scoring
                    else if (parameters_semiGappedScoring) {

                        ungappedExtension->nominalScore = semiGappedScoring_score_multi(
                                ungappedExtension, PSSMatrix, alignment->subjectLength,
                                alignment->subject,
                                statistics_gappedNominalDropoff_multi[queryNum], queryNum);
                    }
                    // Regular gapped scoring
                    else {
                        ungappedExtension->nominalScore = gappedScoring_score_multi(
                                ungappedExtension, PSSMatrix, alignment->subjectLength,
                                alignment->subject,
                                statistics_gappedNominalDropoff_multi[queryNum], queryNum);

                        // Mark as gapped
                        ungappedExtension->status = ungappedExtension_GAPPED;
                    }

                    // If alignment scores above R1 cutoff
                    if (ungappedExtension->nominalScore >=
                            blast_nominalR1cutoff_multi[queryNum]) {
                        if (hasChildren) {
                            // Subject has children so perform stage1 and 2 on children
                            alignments_expandCluster_multi(alignment, PSSMatrix, queryNum);
                            bestScore = 0;
                            break;
                        } else if (ungappedExtension->nominalScore > bestScore) {
                            // Update best score for the alignment
                            bestScore = ungappedExtension->nominalScore;
                        }
                    } else {
                        // Else mark it as deleted
                        ungappedExtension->status = ungappedExtension_DELETED;
                    }


                    // Remove any ungapped extensions in this alignment that are in the
                    // area covered
                    // by the gapped scoring just performed
                    alignments_pruneRegion_multi(alignment, ungappedExtension, queryNum);

                    numExtensions++;
                }

                ungappedExtension = ungappedExtension->next;
            }

            // If this alignment contains gapped extensions that could score above
            // cutoff
            if (bestScore >= blast_nominalR1cutoff_multi[queryNum]) {
                // If a single sequence add to list of "good" alignments
                alignments_addGoodAlignment_multi(bestScore, alignment, queryNum);

                blast_numGoodExtensions_multi[queryNum] += numExtensions;
                blast_numGoodAlignments_multi[queryNum]++;
            }
        }

        // Record point to list of alignments for this volume
        alignments_volumeAlignments_multi2
            [queryNum][alignments_numVolumes_multi[queryNum]] =
            alignments_alignments_multi2[queryNum][blockNum];
        alignments_numVolumes_multi[queryNum]++;

        // Construct new list for next volume (if there is one)
        //alignments_alignments_multi2[queryNum][blockNum] = memBlocks_initialize(
        //    sizeof(struct alignment),
        //    constants_initialAllocAlignments);
    }
}

// Perform initial scoring of all ungapped extensions to find "good" alignments
// that may
// score above the cutoff
void alignments_findGoodAlignments_multi(struct PSSMatrix PSSMatrix,
        int queryNum) {
    struct alignment *alignment;
    struct ungappedExtension *ungappedExtension;
    int4 bestScore, numExtensions, hasChildren;
    //printf("alignments_findGoodAlignments_multi\n");
    // For each alignment
    memBlocks_resetCurrent(alignments_alignments_multi[queryNum]);
    while ((alignment = (struct alignment *)memBlocks_getCurrent(
                    alignments_alignments_multi[queryNum])) != NULL) {
        bestScore = 0;
        blast_dloc_multi[queryNum] = alignment->descriptionLocation;

        // Record if subject has children
        if (encoding_alphabetType == encoding_protein &&
                alignment->encodedLength > alignment->subjectLength + 2)
            hasChildren = 1;
        else
            hasChildren = 0;

        // For each ungapped extension (in descending order of score)
        numExtensions = 0;
        ungappedExtension = alignment->ungappedExtensions;
        while (ungappedExtension != NULL) {
            if (ungappedExtension->status != ungappedExtension_DELETED) {
                // Find the seed
                ungappedExtension_findSeed(ungappedExtension, PSSMatrix,
                        alignment->subject);

                // Byte-packed scoring
                if (parameters_bytepackedScoring) {
                    ungappedExtension->nominalScore = bytepackGappedScoring_score(
                            ungappedExtension, PSSMatrix, alignment->subjectLength,
                            alignment->subject,
                            statistics_gappedNominalDropoff_multi[queryNum]);

                    // Mark as semigapped
                    ungappedExtension->status = ungappedExtension_SEMIGAPPED;
                }
                // Table driven scoring
                else if (parameters_tableScoring) {
                    ungappedExtension->nominalScore = tableGappedScoring_score(
                            ungappedExtension, PSSMatrix, alignment->subjectLength,
                            alignment->subject,
                            statistics_gappedNominalDropoff_multi[queryNum]);

                    // Mark as semigapped
                    ungappedExtension->status = ungappedExtension_SEMIGAPPED;
                }
                // Semi-gapped scoring
                else if (parameters_semiGappedScoring) {
                    ungappedExtension->nominalScore = semiGappedScoring_score_multi(
                            ungappedExtension, PSSMatrix, alignment->subjectLength,
                            alignment->subject,
                            statistics_gappedNominalDropoff_multi[queryNum], queryNum);

                    // Mark as semigapped
                    ungappedExtension->status = ungappedExtension_SEMIGAPPED;
                }
                // Regular gapped scoring
                else {
                    ungappedExtension->nominalScore = gappedScoring_score_multi(
                            ungappedExtension, PSSMatrix, alignment->subjectLength,
                            alignment->subject,
                            statistics_gappedNominalDropoff_multi[queryNum], queryNum);

                    // Mark as gapped
                    ungappedExtension->status = ungappedExtension_GAPPED;
                }

                // If alignment scores above R1 cutoff
                if (ungappedExtension->nominalScore >=
                        blast_nominalR1cutoff_multi[queryNum]) {
                    if (hasChildren) {
                        // Subject has children so perform stage1 and 2 on children
                        alignments_expandCluster_multi(alignment, PSSMatrix, queryNum);
                        bestScore = 0;
                        break;
                    } else if (ungappedExtension->nominalScore > bestScore) {
                        // Update best score for the alignment
                        bestScore = ungappedExtension->nominalScore;
                    }
                } else {
                    // Else mark it as deleted
                    ungappedExtension->status = ungappedExtension_DELETED;
                }

                // Remove any ungapped extensions in this alignment that are in the area
                // covered
                // by the gapped scoring just performed
                alignments_pruneRegion_multi(alignment, ungappedExtension, queryNum);

                numExtensions++;
            }

            ungappedExtension = ungappedExtension->next;
        }

        // If this alignment contains gapped extensions that could score above
        // cutoff
        if (bestScore >= blast_nominalR1cutoff_multi[queryNum]) {
            // If a single sequence add to list of "good" alignments
            alignments_addGoodAlignment_multi(bestScore, alignment, queryNum);

            blast_numGoodExtensions_multi[queryNum] += numExtensions;
            blast_numGoodAlignments_multi[queryNum]++;
        }
    }

    // Record point to list of alignments for this volume
    alignments_volumeAlignments_multi[queryNum]
        [alignments_numVolumes_multi[queryNum]] =
        alignments_alignments_multi[queryNum];
    alignments_numVolumes_multi[queryNum]++;

    // Construct new list for next volume (if there is one)
    alignments_alignments_multi[queryNum] = memBlocks_initialize(
            sizeof(struct alignment), constants_initialAllocAlignments);
}

// Add the current alignment (which contains at least one gapped extension
// with semi-gapped score above semi-gapped cutoff) to to-be-sorted list of good
// alignments
void alignments_addGoodAlignment_multi(int4 highestNominalScore,
        struct alignment *alignment,
        int queryNum) {
    struct finalAlignment *goodAlignment;

    // Get a new good alignment entry
    goodAlignment =
        memSingleBlock_newEntry(alignments_goodAlignments_multi[queryNum]);

    // Insert alignment information into the new good alignment
    goodAlignment->highestNominalScore = highestNominalScore;
    goodAlignment->alignment = alignment;
}

// Sort the array of good alignments in order of score
void alignments_sortGoodAlignments_multi(int queryNum) {
    qsort(alignments_goodAlignments_multi[queryNum]->block,
            alignments_goodAlignments_multi[queryNum]->numEntries,
            sizeof(struct finalAlignment), alignments_compareFinalAlignments);
}

// Remove any ungapped extensions in this alignment that are in the area covered
// by the given gapped scored extension
void alignments_pruneRegion_multi(struct alignment *alignment,
        struct ungappedExtension *ungappedExtension,
        int queryNum) {
    struct ungappedExtension *currentExtension;

#ifdef VERBOSE
    if (alignment->descriptionLocation == parameters_verboseDloc)
        printf("pruneRegion ungappedExtension query=%d to %d subject=%d to %d "
                "score=%d dloc=%d status=%d\n",
                ungappedExtension->start.queryOffset,
                ungappedExtension->end.queryOffset,
                ungappedExtension->start.subjectOffset,
                ungappedExtension->end.subjectOffset,
                ungappedExtension->nominalScore, alignment->descriptionLocation,
                ungappedExtension->status);
#endif

    currentExtension = alignment->ungappedExtensions;
    while (currentExtension != NULL) {
        // If currentExtension start or end is contained within area bound by
        // ungappedExtension
        // (after gappedScoring ungappedExtension), then remove it
        if (currentExtension->status != ungappedExtension_DELETED &&
                alignments_contains(ungappedExtension, currentExtension)) {
#ifdef VERBOSE
            if (alignment->descriptionLocation == parameters_verboseDloc)
                printf("REMOVING %d query %d to %d subject %d to %d\n",
                        currentExtension->nominalScore,
                        currentExtension->start.queryOffset,
                        currentExtension->end.queryOffset,
                        currentExtension->start.subjectOffset,
                        currentExtension->end.subjectOffset);
#endif

            blast_numExtensionsPruned_multi[queryNum]++;

            currentExtension->status = ungappedExtension_DELETED;
        }

        // Advance to next ungapped extension
        currentExtension = currentExtension->next;
    }
}

// Find the top N final alignments
void alignments_findTopFinalAlignments_multi(struct PSSMatrix PSSMatrix,
        int queryNum) {
    struct finalAlignment *goodAlignment, *finalAlignment;
    struct ungappedExtension *ungappedExtension;
    struct alignment *alignment;
    int4 position, bestScore;
    int4 maximumDynamicCutoff;

    //    alignments_printGoodAlignments();

    // Get the Nth alignment
    goodAlignment =
        memSingleBlock_getEntry(alignments_goodAlignments_multi[queryNum],
                parameters_numDisplayAlignments);

    // Calculate the maximum possibly dynamic cutoff
    maximumDynamicCutoff =
        goodAlignment->highestNominalScore / parameters_semiGappedR1;

    // Alignments above that maximum cutoff are definately final alignments
    position = 0;
    while (position < parameters_numDisplayAlignments) {
        // Get each alignment in descending order of score
        goodAlignment = memSingleBlock_getEntry(
                alignments_goodAlignments_multi[queryNum], position);
        alignment = goodAlignment->alignment;

        // Stop when score before maximum cutoff
        if (goodAlignment->highestNominalScore < maximumDynamicCutoff)
            break;

        // Add to final alignments
        alignments_addFinalAlignment_multi(goodAlignment->highestNominalScore,
                alignment, queryNum);

        position++;
    }

    //    alignments_printFinalAlignments();

    // For the rest of the top N alignments
    while (position < parameters_numDisplayAlignments) {
        // Continue getting each alignment in descending order of score
        goodAlignment = memSingleBlock_getEntry(
                alignments_goodAlignments_multi[queryNum], position);
        alignment = goodAlignment->alignment;

        //        printf("Alignment score=%d\n",
        // goodAlignment->highestNominalScore);

        // Perform regular gapped scoring
        bestScore = 0;
        ungappedExtension = alignment->ungappedExtensions;
        while (ungappedExtension != NULL) {
            // Find the highest scoring semigapped alignment
            if (ungappedExtension->status != ungappedExtension_DELETED &&
                    ungappedExtension->nominalScore ==
                    goodAlignment->highestNominalScore) {
                // Perform gapped alignment
                alignments_regularGappedAlignment_multi(PSSMatrix, ungappedExtension,
                        alignment, queryNum);

                // Check for join to another extension
                alignments_checkForJoin_multi(alignment, ungappedExtension, PSSMatrix,
                        queryNum);

                // Update the alignment's best score
                bestScore = ungappedExtension->nominalScore;
            }

            ungappedExtension = ungappedExtension->next;
        }

        // Update the alignment's score
        goodAlignment->highestNominalScore = bestScore;

        // Add to final alignments
        alignments_addFinalAlignment_multi(goodAlignment->highestNominalScore,
                alignment, queryNum);

        position++;
    }

    //    alignments_printFinalAlignments();

    // Initialize dynamic cutoffs
    blast_dynamicGappedNominalCutoff_multi[queryNum] = 0;
    blast_dynamicNominalR1cutoff_multi[queryNum] = 0;

    // For the remaining good alignments
    while (position < alignments_goodAlignments_multi[queryNum]->numEntries) {
        // Continue getting each alignment in descending order of score
        goodAlignment = memSingleBlock_getEntry(
                alignments_goodAlignments_multi[queryNum], position);
        alignment = goodAlignment->alignment;

        // Stop when below dynamic cutoff
        if (goodAlignment->highestNominalScore <
                blast_dynamicNominalR1cutoff_multi[queryNum])
            break;

        // Perform regular gapped scoring
        bestScore = 0;
        ungappedExtension = alignment->ungappedExtensions;
        while (ungappedExtension != NULL) {
            // Find the highest scoring semigapped alignment
            if (ungappedExtension->status != ungappedExtension_DELETED &&
                    ungappedExtension->nominalScore ==
                    goodAlignment->highestNominalScore) {
                // Perform gapped alignment
                alignments_regularGappedAlignment_multi(PSSMatrix, ungappedExtension,
                        alignment, queryNum);

                // Check for join to another extension
                alignments_checkForJoin_multi(alignment, ungappedExtension, PSSMatrix,
                        queryNum);

                // Update the alignment's best score
                bestScore = ungappedExtension->nominalScore;
            }

            ungappedExtension = ungappedExtension->next;
        }

        //        printf("Position=%d OldScore=%d NewScore=%d DynamicCutoff=%d
        // FastCutoff=%d\n", position,
        //               goodAlignment->highestNominalScore, bestScore,
        // blast_dynamicGappedNominalCutoff_multi[queryNum],
        //               blast_dynamicNominalR1cutoff_multi[queryNum]);

        // Update the alignment's score
        goodAlignment->highestNominalScore = bestScore;

        // If it scores high enough to join the top N final alignments
        if (bestScore > blast_dynamicGappedNominalCutoff_multi[queryNum])

            // If it score above cutoff
            //		if (bestScore >
            // blast_gappedNominalCutoff_multi[queryNum])
        {
            // Add to final alignments
            alignments_addFinalAlignment_multi(goodAlignment->highestNominalScore,
                    alignment, queryNum);

            // Re-sort final alignments and final minimum score of top N
            alignments_sortFinalAlignments_multi(queryNum);
            finalAlignment =
                memSingleBlock_getEntry(alignments_finalAlignments_multi[queryNum],
                        parameters_numDisplayAlignments - 1);

            // Calculate dynamic cutoffs
            blast_dynamicGappedNominalCutoff_multi[queryNum] =
                finalAlignment->highestNominalScore;
            blast_dynamicNominalR1cutoff_multi[queryNum] =
                blast_dynamicGappedNominalCutoff_multi[queryNum] *
                parameters_semiGappedR1;
        }

        position++;
    }
}

// Given a collection of good alignments, find the final top N alignments above
// cutoff
void alignments_findFinalAlignments_multi(struct PSSMatrix PSSMatrix,
        int queryNum) {
    struct finalAlignment *goodAlignment;
    struct alignment *alignment;
    struct ungappedExtension *ungappedExtension;
    int4 bestScore;

    // Sort good alignments by score
    // TODO: How long does this take?
    alignments_sortGoodAlignments_multi(queryNum);

    // If we have more good alignments than we can display
    if (parameters_numDisplayAlignments > 0 &&
            alignments_goodAlignments_multi[queryNum]->numEntries >
            parameters_numDisplayAlignments) {
        goodAlignment =
            memSingleBlock_getEntry(alignments_goodAlignments_multi[queryNum],
                    parameters_numDisplayAlignments - 1);

        // If it is clear we will have more final alignments than we can display
        if (goodAlignment->highestNominalScore >
                blast_nominalR2cutoff_multi[queryNum]) {
            // Find the top N alignments
            alignments_findTopFinalAlignments_multi(PSSMatrix, queryNum);
            return;
        }
    }

    // Further score all good alignments to see if they score above cutoff
    memSingleBlock_resetCurrent(alignments_goodAlignments_multi[queryNum]);
    while ((goodAlignment = memSingleBlock_getCurrent(
                    alignments_goodAlignments_multi[queryNum])) != NULL) {
        alignment = goodAlignment->alignment;

        bestScore = 0;
        // For each ungapped extension that hasn't been deleted
        ungappedExtension = alignment->ungappedExtensions;
        while (ungappedExtension != NULL) {
#ifdef VERBOSE
            if (parameters_verboseDloc == alignment->descriptionLocation)
                ungappedExtension_print(ungappedExtension);
#endif

            if (ungappedExtension->status != ungappedExtension_DELETED) {
                // If it isn't clear if the extension is above cutoff
                if (ungappedExtension->nominalScore >=
                        blast_nominalR1cutoff_multi[queryNum] &&
                        ungappedExtension->nominalScore <=
                        blast_nominalR2cutoff_multi[queryNum]) {
                    // Perform regular gapped scoring
                    alignments_regularGappedAlignment_multi(PSSMatrix, ungappedExtension,
                            alignment, queryNum);
                }

                // Update the alignment's best score
                if (ungappedExtension->nominalScore > bestScore)
                    bestScore = ungappedExtension->nominalScore;
            }
            ungappedExtension = ungappedExtension->next;
        }

        goodAlignment->highestNominalScore = bestScore;

        // If the alignment scores above cutoff, include in list of final alignments
        if (goodAlignment->highestNominalScore >=
                blast_gappedNominalCutoff_multi[queryNum]) {
            alignments_addFinalAlignment_multi(goodAlignment->highestNominalScore,
                    alignment, queryNum);
        }
    }
}

// Get the tracebacks for all of the final alignments
void alignments_getTracebacks_multi(struct PSSMatrix PSSMatrix, int queryNum) {
    struct finalAlignment *finalAlignment;
    struct gappedExtension *gappedExtension;
    struct ungappedExtension *ungappedExtension, *highestScoringExtension,
                             oldUngappedExtension;
    struct alignment *alignment;
    struct unpackRegion *unpackRegion;
    int4 numProcessed = 0, numAboveCutoff = 0, repeatComputeTracebacks = 1;
    int4 alignmentCount;
    unsigned char *subject;

    // Sort final alignments by score
    alignments_sortFinalAlignments_multi(queryNum);

    //    alignments_printFinalAlignments();

    // Only keep alignments above cutoff
    while (alignments_finalAlignments_multi[queryNum]->numEntries > 0) {
        // Get last alignment in list
        finalAlignment =
            memSingleBlock_getLastEntry(alignments_finalAlignments_multi[queryNum]);

        // Stop if above cutoff
        if (finalAlignment->highestNominalScore != 0)
            break;

        // Otherwise remove it from the list
        alignments_finalAlignments_multi[queryNum]->numEntries--;
    }

    // For each alignment that is in the top numDisplayAlignments but not the top
    // numDisplayTracebacks
    alignmentCount = parameters_numDisplayTracebacks;
    while (alignmentCount < parameters_numDisplayAlignments &&
            alignmentCount <
            alignments_finalAlignments_multi[queryNum]->numEntries) {
        finalAlignment = memSingleBlock_getEntry(
                alignments_finalAlignments_multi[queryNum], alignmentCount);
        alignment = finalAlignment->alignment;
        blast_dloc_multi[queryNum] = alignment->descriptionLocation;

        // Get the highest scoring ungapped extension
        ungappedExtension = alignment->ungappedExtensions;
        highestScoringExtension = ungappedExtension;
        while (ungappedExtension != NULL) {
            if (ungappedExtension->nominalScore >
                    highestScoringExtension->nominalScore &&
                    ungappedExtension->status != ungappedExtension_DELETED) {
                highestScoringExtension = ungappedExtension;
            }
            ungappedExtension = ungappedExtension->next;
        }

        if (highestScoringExtension != NULL) {
            subject = alignments_selectRegion(alignment, highestScoringExtension);

            // Perform gapped scoring with higher dropoff
            highestScoringExtension->nominalScore = gappedScoring_score_multi(
                    highestScoringExtension, PSSMatrix, alignment->subjectLength, subject,
                    statistics_gappedFinalNominalDropoff_multi[queryNum], queryNum);

            finalAlignment->highestNominalScore =
                highestScoringExtension->nominalScore;
        }

        //        printf("Rescore with larger dropoff num %d: Score=%d\n",
        // alignmentCount, highestScoringExtension->nominalScore);
        alignmentCount++;
    }

    while (repeatComputeTracebacks) {
        numAboveCutoff = 0;
        numProcessed = 0;
        memSingleBlock_resetCurrent(alignments_finalAlignments_multi[queryNum]);
        while (((finalAlignment = memSingleBlock_getCurrent(
                            alignments_finalAlignments_multi[queryNum])) != NULL) &&
                (numAboveCutoff < parameters_numDisplayTracebacks ||
                 parameters_numDisplayTracebacks == 0)) {
            alignment = finalAlignment->alignment;
            blast_dloc_multi[queryNum] = alignment->descriptionLocation;

            // If traceback haven't been computed for this alignment
            if (alignment->gappedExtensions == NULL) {
                // Unpack part or all of the subject (for nucleotide)
                // blast_gappedExtendTime += clock();
                // blast_unpackTime -= clock();
                unpack_unpackSubject_multi(PSSMatrix, alignment, queryNum);
                // blast_unpackTime += clock();
                // blast_gappedExtendTime -= clock();

                // For each ungapped extension that hasn't been deleted
                ungappedExtension = alignment->ungappedExtensions;
                while (ungappedExtension != NULL) {
                    // If extension scores above cutoff
                    if (ungappedExtension->status != ungappedExtension_DELETED) {
                        // Make copy of ungapped extension
                        oldUngappedExtension = *ungappedExtension;

                        // If subject and query are short enough and sequence does not have
                        // multiple
                        // unpack regions, use faster but less memory efficient gapped
                        // alignment with traceback
                        // if (((uint8)PSSMatrix.length * (uint8)alignment->subjectLength <
                        //(uint8)constants_maximumTracebackSize) &&
                        // unpack_entireSubjectUnpacked(alignment))
                        if (1) {
                            gappedExtension = fasterGappedExtension_build_multi(
                                    ungappedExtension, PSSMatrix, alignment->subjectLength,
                                    alignment->unpackRegions[0].unpackedSubject,
                                    statistics_gappedFinalNominalDropoff_multi[queryNum],
                                    queryNum);
                        }
                        // Otherwise use slower but more memory-efficient gapped alignment
                        else {
                            unpackRegion = unpack_selectRegion(
                                    alignment->unpackRegions, alignment->numUnpackRegions,
                                    ungappedExtension->seed.subjectOffset);

                            gappedExtension = gappedExtension_build_multi(
                                    ungappedExtension, PSSMatrix, alignment->subjectLength,
                                    alignment->subject, unpackRegion,
                                    statistics_gappedFinalNominalDropoff_multi[queryNum],
                                    queryNum);
                        }

                        // Calculate normalized score and e-value
                        gappedExtension_score_multi(gappedExtension, queryNum);

                        // Add it to the current alignment
                        if (gappedExtension->nominalScore >=
                                blast_gappedNominalCutoff_multi[queryNum])
                            alignments_addGappedExtension(alignment, gappedExtension);
                        else {
                            free(gappedExtension->trace.traceCodes);
                            free(gappedExtension);
                        }

                        // Check for ungapped extensions that were mistakenly pruned
                        if (parameters_tableScoring)
                            alignments_unpruneRegion_multi(alignment, &oldUngappedExtension,
                                    ungappedExtension, PSSMatrix,
                                    queryNum);

                        // Remove any ungapped extensions in this alignment that are in the
                        // area covered
                        // by the gapped scoring just performed
                        alignments_pruneRegion_multi(alignment, ungappedExtension,
                                queryNum);
                    }
                    ungappedExtension = ungappedExtension->next;
                }

                // Finally prune extensions that share a start or end point
                alignments_pruneOverlappingExtensions(alignment);

                //                printf("Was %d Now %d\n",
                // finalAlignment->highestNominalScore,
                // alignment->gappedExtensions->nominalScore);

                // Update final alignment's high-score to that of the first (and
                // highest-scoring)
                // gapped extension in the list
                if (alignment->gappedExtensions != NULL)
                    finalAlignment->highestNominalScore =
                        alignment->gappedExtensions->nominalScore;
                else
                    finalAlignment->highestNominalScore = 0;

                numProcessed++;
                //	            printf("Computed alignment score=%d\n",
                // finalAlignment->highestNominalScore);
            }

            // Tally number of final alignments above cutoff
            if (finalAlignment->highestNominalScore > 0)
                numAboveCutoff++;
        }

        //    	printf("Traceback alignments performed=%d\n", numProcessed);

        // Sort final alignments by score
        alignments_sortFinalAlignments_multi(queryNum);

        //        printf("repeatComputeTracebacks:");

        // If the first numDisplayTracebacks alignments have traceback computed,
        // stop
        numProcessed = 0;
        repeatComputeTracebacks = 0;
        memSingleBlock_resetCurrent(alignments_finalAlignments_multi[queryNum]);
        while ((finalAlignment = memSingleBlock_getCurrent(
                        alignments_finalAlignments_multi[queryNum])) != NULL &&
                numProcessed < parameters_numDisplayTracebacks) {
            alignment = finalAlignment->alignment;
            if (alignment->gappedExtensions == NULL &&
                    finalAlignment->highestNominalScore != 0) {
                //            	printf("1");
                repeatComputeTracebacks = 1;
                break;
            }
            //            printf("0");

            numProcessed++;
        }
    }

    // Only keep top N alignments
    if (parameters_numDisplayAlignments != 0 &&
            alignments_finalAlignments_multi[queryNum]->numEntries >
            parameters_numDisplayAlignments) {
        alignments_finalAlignments_multi[queryNum]->numEntries =
            parameters_numDisplayAlignments;
    }

    // Only keep alignments above cutoff
    while (alignments_finalAlignments_multi[queryNum]->numEntries > 0) {
        // Get last alignment in list
        finalAlignment =
            memSingleBlock_getLastEntry(alignments_finalAlignments_multi[queryNum]);

        // Stop if above cutoff
        if (finalAlignment->highestNominalScore != 0)
            break;

        // Otherwise remove it from the list
        alignments_finalAlignments_multi[queryNum]->numEntries--;
    }
}

// Get the subject sequence descriptions for all of the final alignments
void alignments_getFinalAlignmentDescriptions_multi(int queryNum) {
    struct finalAlignment *finalAlignment;

    // Sort descriptions in order of description location (to speed up disk
    // access)
    qsort(alignments_finalAlignments_multi[queryNum]->block,
            alignments_finalAlignments_multi[queryNum]->numEntries,
            sizeof(struct finalAlignment),
            alignments_compareAlignmentDescriptionLocations);

    // For each alignment, read the description
    memSingleBlock_resetCurrent(alignments_finalAlignments_multi[queryNum]);
    while ((finalAlignment = memSingleBlock_getCurrent(
                    alignments_finalAlignments_multi[queryNum])) != NULL) {
        finalAlignment->description = descriptions_getDescription_multi(
                finalAlignment->alignment->descriptionLocation,
                finalAlignment->alignment->descriptionLength, queryNum);
    }

    // Re-sort the final alignments by score
    qsort(alignments_finalAlignments_multi[queryNum]->block,
            alignments_finalAlignments_multi[queryNum]->numEntries,
            sizeof(struct finalAlignment), alignments_compareFinalAlignments);
}

// Sort the array of final alignments in order of score
void alignments_sortFinalAlignments_multi(int queryNum) {
    qsort(alignments_finalAlignments_multi[queryNum]->block,
            alignments_finalAlignments_multi[queryNum]->numEntries,
            sizeof(struct finalAlignment), alignments_compareFinalAlignments);
}

void alignments_sortFinalAlignments_multi2(int queryNum) {
    qsort(alignments_finalAlignments_multi[queryNum]->block,
            alignments_finalAlignments_multi[queryNum]->numEntries,
            sizeof(struct finalAlignment), alignments_compareFinalAlignments2);
}

// Add the current alignment (which contains at least one gapped extension
// scoring above cutoff) to to-be-sorted list of final alignments
struct finalAlignment * alignments_addFinalAlignment_multi(int4 highestNominalScore,
        struct alignment *alignment,
        int queryNum) {
    struct finalAlignment *finalAlignment;

    // Get a new final alignment entry
    finalAlignment =
        memSingleBlock_newEntry(alignments_finalAlignments_multi[queryNum]);

    // Insert alignment information into the new final alignment
    finalAlignment->highestNominalScore = highestNominalScore;
    finalAlignment->alignment = alignment;
    finalAlignment->description = NULL;
    return finalAlignment;
}

// Free memory used to store alignments
void alignments_free_db_multi(int queryNum) {
    struct alignment *alignment;
    struct finalAlignment *finalAlignment;
    uint4 volumeCount;

    memBlocks_resetCurrent(alignments_alignments_multi[queryNum]);

    // Free unused block of alignments
    memBlocks_free(alignments_alignments_multi[queryNum]);

    // For each volume
    volumeCount = 0;
    while (volumeCount < alignments_numVolumes_multi[queryNum]) {
        alignments_alignments_multi[queryNum] =
            alignments_volumeAlignments_multi[queryNum][volumeCount];

        // For each alignment
        memBlocks_resetCurrent(alignments_alignments_multi[queryNum]);
        while ((alignment = (struct alignment *)memBlocks_getCurrent(
                        alignments_alignments_multi[queryNum])) != NULL) {
            alignment_freeAlignment(alignment);
        }

        // Free blocks of alignments
        memBlocks_free(alignments_alignments_multi[queryNum]);

        volumeCount++;
    }

    // Free memory used by ungapped extensions
    //int ii;
    //for (ii = 0; ii < blast_tgSize; ii++)
    //memBlocks_free(
    //ungappedExtension_extensions_multi[queryNum * blast_tgSize + ii]);

    // Free good alignments
    memSingleBlock_free(alignments_goodAlignments_multi[queryNum]);

    // For each final alignment, free the description
    memSingleBlock_resetCurrent(alignments_finalAlignments_multi[queryNum]);
    while ((finalAlignment = memSingleBlock_getCurrent(
                    alignments_finalAlignments_multi[queryNum])) != NULL) {
        free(finalAlignment->description);
    }

    // Free final alignments
    memSingleBlock_free(alignments_finalAlignments_multi[queryNum]);

    // Free all the unpacked regions
    unpack_free_multi(queryNum);
}

// Free memory used to store alignments
void alignments_free_multi2(int queryNum) {
    struct alignment *alignment;
    struct finalAlignment *finalAlignment;
    uint4 volumeCount;

    // For each volume
    volumeCount = 0;
    while (volumeCount < alignments_numVolumes_multi[queryNum]) {
        alignments_alignments =
            alignments_volumeAlignments_multi2[queryNum][volumeCount];

        // For each alignment
        memBlocks_resetCurrent(alignments_alignments);
        while ((alignment = (struct alignment *)memBlocks_getCurrent(
                        alignments_alignments)) != NULL) {
            alignment_freeAlignment(alignment);
        }

        // Free blocks of alignments
        //memBlocks_free(alignments_alignments);

        volumeCount++;
    }

    //int blockNum;
    //for (blockNum = 0; blockNum < blast_numBlocks; blockNum++) {
    //memBlocks_resetCurrent(alignments_alignments_multi2[queryNum][blockNum]);

    //// Free unused block of alignments
    //memBlocks_free(alignments_alignments_multi2[queryNum][blockNum]);
    //}

    //memBlocks_resetCurrent(alignments_alignments_multi[queryNum]);

    // Free unused block of alignments
    //memBlocks_free(alignments_alignments_multi[queryNum]);

    //for (blockNum = 0; blockNum < blast_numBlocks; blockNum++) {
    //memBlocks_free(ungappedExtension_extensions_multi2[queryNum][blockNum]);
    //}

    // Free memory used by ungapped extensions
    //memBlocks_free(ungappedExtension_extensions_multi[queryNum]);

    // Free good alignments
    memSingleBlock_free(alignments_goodAlignments_multi[queryNum]);

    // For each final alignment, free the description
    memSingleBlock_resetCurrent(alignments_finalAlignments_multi[queryNum]);
    //while ((finalAlignment = memSingleBlock_getCurrent(
                    //alignments_finalAlignments_multi[queryNum])) != NULL) {
        //free(finalAlignment->description);
    //}

    // Free final alignments
    memSingleBlock_free(alignments_finalAlignments_multi[queryNum]);

    // Free all the unpacked regions
    unpack_free_multi(queryNum);
}

void alignments_global_free_multi2() {
    int ii;
    for (ii = 0; ii < blast_numQuery; ii++) {
        free(ungappedExtension_extensions_multi2[ii]);
        free(alignments_alignments_multi2[ii]);
        free(alignments_volumeAlignments_multi2[ii]);
    }

    free(ungappedExtension_extensions_multi2);
    free(alignments_alignments_multi2);
    free(alignments_volumeAlignments_multi2);
}

// Free memory used to store alignments
void alignments_free_multi(int queryNum) {
    struct alignment *alignment;
    struct finalAlignment *finalAlignment;
    uint4 volumeCount;

    //memBlocks_resetCurrent(alignments_alignments_multi[queryNum]);

    // Free unused block of alignments
    //memBlocks_free(alignments_alignments_multi[queryNum]);

    // For each volume
    //volumeCount = 0;
    //while (volumeCount < alignments_numVolumes_multi[queryNum]) {
    //alignments_alignments_multi[queryNum] =
    //alignments_volumeAlignments_multi[queryNum][volumeCount];

    //// For each alignment
    //memBlocks_resetCurrent(alignments_alignments_multi[queryNum]);
    //while ((alignment = (struct alignment *)memBlocks_getCurrent(
    //alignments_alignments_multi[queryNum])) != NULL) {
    //alignment_freeAlignment(alignment);
    //}

    //// Free blocks of alignments
    //memBlocks_free(alignments_alignments_multi[queryNum]);

    //volumeCount++;
    //}

    // Free memory used by ungapped extensions
    //int ii;
    //for (ii = 0; ii < blast_tgSize; ii++)
    //memBlocks_free(
    //ungappedExtension_extensions_multi[queryNum * blast_tgSize + ii]);

    // Free good alignments
    memSingleBlock_free(alignments_goodAlignments_multi[queryNum]);

    // For each final alignment, free the description
    memSingleBlock_resetCurrent(alignments_finalAlignments_multi[queryNum]);
    while ((finalAlignment = memSingleBlock_getCurrent(
                    alignments_finalAlignments_multi[queryNum])) != NULL) {
        free(finalAlignment->description);
    }

    // Free final alignments
    memSingleBlock_free(alignments_finalAlignments_multi[queryNum]);

    // Free all the unpacked regions
    unpack_free_multi(queryNum);
}

void alignments_query_serial(struct PSSMatrix *PSSMatrix_arr, 
        int numQuery, 
        struct scoreMatrix scoreMatrix) {
    int ii, jj;

    struct timeval start, end;

    hit_t *firstBin = (hit_t *)malloc(sizeof(hit_t) * MAX_NUM_HIT_PER_QUERY * numQuery);
    hit_t *secondBin = (hit_t *)malloc(sizeof(hit_t) * MAX_NUM_HIT_PER_QUERY * numQuery);

    int maxNumFirstBin = longestQueryLength + readdb_longestSequenceLength + 1;
    int numSecondBin = numQuery;

    int4 *binOffset = (int4 *)malloc(sizeof(int4) * maxNumFirstBin);
    uint4 *numExtHit = (uint4 *)malloc(sizeof(uint4) * numSecondBin);
    uint4 *numHit = (uint4 *)malloc(sizeof(uint4) * numSecondBin);
    hit_t *lastHits = (hit_t *)malloc(sizeof(hit_t) * numSecondBin);

    //fprintf(stderr, "bin size: %d MB\n", (sizeof(hit_t) * MAX_NUM_HIT_PER_QUERY * numQuery)>>20);

    gettimeofday(&start, NULL);

    int goodExtensionBufSize = MAX_EXTENSIONS_PER_QUERY * numQuery;

    struct ungappedExtension *goodExtensionBuf = (struct ungappedExtension *)malloc(sizeof(struct ungappedExtension) * goodExtensionBufSize);


    long long goodAlignBufSize = MAX_ALIGNMENTS_PER_QUERY * numQuery;

    struct alignment *goodAlignBuf = (struct alignment *)malloc(sizeof(struct alignment) * goodAlignBufSize);


    struct ungappedExtension **ungappedExtension_new = (struct ungappedExtension **)malloc(sizeof(struct ungappedExtension *) * MAX_NUM_TRIGGER_EXT );


    BlastGapDP *dp_mem = (BlastGapDP *)malloc(sizeof(BlastGapDP) * DP_MEM_SIZE );

    BlastIntervalTree *tree = Blast_IntervalTreeInit(0, longestQueryLength + 1,
            0, SEMIGAPPED_ROWSIZE + 1);

    BlastIntervalTree *private_tree = Blast_IntervalTreeInit(0, 
            longestQueryLength + 1, 0, SEMIGAPPED_ROWSIZE + 1);

    BlastHSP *BlastHSPs = (BlastHSP *)malloc(sizeof(BlastHSP) * MAX_NUM_HSP);


    int alphsize = COMPO_NUM_TRUE_AA;
    int n  = alphsize * alphsize;
    int mA = 2 * alphsize - 1;
    int m  = mA + 1;

    Blast_CompositionWorkspace *NRrecord = Blast_CompositionWorkspaceNew();
    Blast_MatrixInfo matrixInfo;
    Blast_MatrixInfoNew2(&matrixInfo, BLASTAA_SIZE, BLASTAA_SIZE, FALSE);

    BlastGapAlignStruct *gap_align;
    BLAST_GapAlignStructNew(readdb_longestSequenceLength, &gap_align);

    int4 **matrix = (void *)_PSIAllocateMatrix(BLASTAA_SIZE, BLASTAA_SIZE, sizeof(Int4)); 

    unsigned char *seq_data;
    seq_data = (unsigned char *)malloc((readdb_longestSequenceLength + 10) * sizeof(unsigned char));

    ReNewtonSystem *newton_system;
    newton_system = ReNewtonSystemNew(alphsize); 

    double **grads, **Scores, *resids_x, *resids_z, *old_scores, *workspace, *z;
    grads = Nlm_DenseMatrixNew(2, n);
    Scores = Nlm_DenseMatrixNew(BLASTAA_SIZE, BLASTAA_SIZE);

    resids_x = (double *) malloc(n * sizeof(double));
    resids_z = (double *) malloc((mA + 1) * sizeof(double));
    z = (double *) calloc( mA + 1,   sizeof(double));
    old_scores = (double *) malloc(n * sizeof(double));
    workspace = (double *) malloc(n * sizeof(double));

    char **querySeq_arr = (char **)malloc(sizeof(char *) * numQuery);
    Blast_AminoAcidComposition *query_composition_arr = (Blast_AminoAcidComposition *)malloc(sizeof(Blast_AminoAcidComposition) * numQuery);

    for(ii = 0; ii < numQuery; ii++)
    {
        querySeq_arr[ii] = (char *)malloc(sizeof(char) * PSSMatrix_arr[ii].length);
        for(jj = 0; jj < PSSMatrix_arr[ii].length; jj++)
        {
            ASSERT(PSSMatrix_arr[ii].queryCodes[jj] < FSA_AA_SIZE);
            querySeq_arr[ii][jj] = to_ncbi[PSSMatrix_arr[ii].queryCodes[jj]];
        }
        Blast_ReadAaComposition_fsa(&query_composition_arr[ii], BLASTAA_SIZE, PSSMatrix_arr[ii].queryCodes, PSSMatrix_arr[ii].length);
    }

    struct gappedExtension *gappedExtension = (struct gappedExtension *)malloc(sizeof(struct gappedExtension) * MAX_NUM_GAPEXT_PER_QUERY * numQuery );

    int goodExtensionCount = 0;
    int goodAlignCount = 0;
    int gappedExtension_cnt = 0;

#if 1
    search_protein2hit_lookup_multi(PSSMatrix_arr, readdb_sequenceData,
            readdb_numberOfSequences, numQuery,
            firstBin, secondBin, binOffset, 
            numExtHit, numHit, lastHits, 
            scoreMatrix, 
            goodExtensionBuf, 
            &goodExtensionCount,
            goodAlignBuf,
            &goodAlignCount,
            ungappedExtension_new,
            dp_mem,
            tree,
            private_tree,
            BlastHSPs
            );
#else
    for(ii = 0; ii < numQuery; ii++)
        search_protein2hit_multi(PSSMatrix_arr, readdb_sequenceData,
                readdb_numberOfSequences, ii);
#endif

#ifndef NO_STAGE4
    for (jj = 0; jj < numQuery; jj++) {
        alignments_getTracebacks_ncbi_multi(PSSMatrix_arr[jj], scoreMatrix, &query_composition_arr[jj], jj, NRrecord, querySeq_arr[jj], &matrixInfo, gap_align, matrix, seq_data, newton_system, grads, Scores, z, resids_x, resids_z, old_scores, workspace, gappedExtension, &gappedExtension_cnt);

        if(gappedExtension_cnt >= MAX_NUM_GAPEXT_PER_QUERY * numQuery)
        {
            fprintf(stderr, "gappedExtension_arr overflow\n");
            exit(1);
        }
    }
#endif


    gettimeofday(&end, NULL);
    long hit_time = ((end.tv_sec * 1000000 + end.tv_usec) -
            (start.tv_sec * 1000000 + start.tv_usec));
    fprintf(stderr, "Hit time: %f\n", (float)hit_time * 1e-6);


    free(numExtHit);
    free(numHit);
    free(lastHits);
    free(binOffset);
    free(firstBin);
    free(secondBin);
    free(BlastHSPs);
    free(dp_mem);
    free(ungappedExtension_new);

    free(workspace);
    free(old_scores);
    free(z);
    free(resids_z);
    free(resids_x);
    Nlm_DenseMatrixFree(&grads);
    Nlm_DenseMatrixFree(&Scores);
    ReNewtonSystemFree(&newton_system);
    free(seq_data);
    _PSIDeallocateMatrix((void **)matrix, BLASTAA_SIZE);
    BLAST_GapAlignStructFree(gap_align);
    Blast_CompositionWorkspaceFree(&NRrecord);
    Blast_MatrixInfoFree2(&matrixInfo);

    free(query_composition_arr);

    Blast_IntervalTreeFree2(private_tree);
    Blast_IntervalTreeFree2(tree);

    free(goodExtensionBuf);
    free(goodAlignBuf);
}

void alignments_query_omp(struct PSSMatrix *PSSMatrix_arr, int numQuery) {
    int jj;
#pragma omp parallel num_threads(parameters_num_threads) default(none) shared( \
        numQuery, PSSMatrix_arr, readdb_sequenceData, readdb_numberOfSequences)
    {
        int4 thread_id = omp_get_thread_num();
        printf("fsablast_omp: affinity - thread_id %d bind to core_id: %d\n",
                thread_id, sched_getcpu());

#pragma omp for schedule(dynamic)
        for (jj = 0; jj < numQuery; jj++) {
            search_protein2hit_multi(PSSMatrix_arr, readdb_sequenceData,
                    readdb_numberOfSequences, jj);
        }

#pragma omp for schedule(dynamic)
        for (jj = 0; jj < numQuery; jj++) {
            alignments_findGoodAlignments_multi(PSSMatrix_arr[jj], jj);
            alignments_findFinalAlignments_multi(PSSMatrix_arr[jj], jj);
            alignments_getFinalAlignmentDescriptions_multi(jj);
            alignments_getTracebacks_multi(PSSMatrix_arr[jj], jj);
        }
    }
}

void alignments_db_serial(struct PSSMatrix *PSSMatrix_arr,
        struct scoreMatrix scoreMatrix, char **query_arr, char **queryDescription_arr, int numQuery) {
    int ii, jj;

    struct timeval start_t, end_t;
    gettimeofday(&start_t, NULL);

    hit_t *hits = (hit_t *)malloc(sizeof(hit_t) * MAX_HITS_PER_SEQ);
    struct ungappedExtension **ungappedExtension_sort, **ungappedExtension_new;
    ungappedExtension_sort = (struct ungappedExtension **)malloc(sizeof(struct ungappedExtension *) * MAX_NUM_TRIGGER_EXT );
    ungappedExtension_new = (struct ungappedExtension **)malloc(sizeof(struct ungappedExtension *) * MAX_NUM_TRIGGER_EXT );
    BlastGapDP *dp_mem = (BlastGapDP *)malloc(sizeof(BlastGapDP) * DP_MEM_SIZE );
    BlastHSP *BlastHSP_arr = (BlastHSP *)malloc(sizeof(BlastHSP) * MAX_NUM_HSP );

    int alphsize = COMPO_NUM_TRUE_AA;
    int n  = alphsize * alphsize;
    int mA = 2 * alphsize - 1;
    int m  = mA + 1;

    Blast_CompositionWorkspace *NRrecord = Blast_CompositionWorkspaceNew();
    BlastIntervalTree *tree = Blast_IntervalTreeInit(0, longestQueryLength + 1,
            0, SEMIGAPPED_ROWSIZE + 1);
    BlastIntervalTree *private_tree = Blast_IntervalTreeInit(0, 
            longestQueryLength + 1, 0, SEMIGAPPED_ROWSIZE + 1);

    Blast_MatrixInfo matrixInfo;
    Blast_MatrixInfoNew2(&matrixInfo, BLASTAA_SIZE, BLASTAA_SIZE, FALSE);

    BlastGapAlignStruct *gap_align;
    BLAST_GapAlignStructNew(readdb_longestSequenceLength, &gap_align);

    int4 **matrix = (void *)_PSIAllocateMatrix(BLASTAA_SIZE, BLASTAA_SIZE, sizeof(Int4)); 

    unsigned char *seq_data;
    seq_data = (unsigned char *)malloc((readdb_longestSequenceLength + 10) * sizeof(unsigned char));

    ReNewtonSystem *newton_system;
    newton_system = ReNewtonSystemNew(alphsize); 

    double **grads, **Scores, *resids_x, *resids_z, *old_scores, *workspace, *z;
    grads = Nlm_DenseMatrixNew(2, n);
    Scores = Nlm_DenseMatrixNew(BLASTAA_SIZE, BLASTAA_SIZE);

    resids_x = (double *) malloc(n * sizeof(double));
    resids_z = (double *) malloc((mA + 1) * sizeof(double));
    z = (double *) calloc( mA + 1,   sizeof(double));
    old_scores = (double *) malloc(n * sizeof(double));
    workspace = (double *) malloc(n * sizeof(double));

    char **querySeq_arr = (char **)malloc(sizeof(char *) * numQuery);
    Blast_AminoAcidComposition *query_composition_arr = (Blast_AminoAcidComposition *)malloc(sizeof(Blast_AminoAcidComposition) * numQuery);

    for(ii = 0; ii < numQuery; ii++)
    {
        querySeq_arr[ii] = (char *)malloc(sizeof(char) * PSSMatrix_arr[ii].length);
        for(jj = 0; jj < PSSMatrix_arr[ii].length; jj++)
        {
            ASSERT(PSSMatrix_arr[ii].queryCodes[jj] < FSA_AA_SIZE);
            querySeq_arr[ii][jj] = to_ncbi[PSSMatrix_arr[ii].queryCodes[jj]];
        }
        Blast_ReadAaComposition_fsa(&query_composition_arr[ii], BLASTAA_SIZE, PSSMatrix_arr[ii].queryCodes, PSSMatrix_arr[ii].length);
    }

    int numQueriesPerThread = numQuery;
    long long goodAlignBufSize = MAX_ALIGNMENTS_PER_QUERY * numQueriesPerThread ;
    int goodExtensionBufSize = MAX_EXTENSIONS_PER_QUERY * numQueriesPerThread ;
#if defined(__ICC) || defined(__INTEL_COMPILER)
    hit_t *binnedHits = (hit_t *)_mm_malloc(
            sizeof(hit_t) * MAX_HITS_PER_SEQ , 64);
    struct alignment *goodAlignBuf = (struct alignment *)_mm_malloc(sizeof(struct alignment) * goodAlignBufSize, 64);
    struct ungappedExtension *goodExtensionBuf = (struct ungappedExtension *)_mm_malloc(sizeof(struct ungappedExtension) * goodExtensionBufSize, 64);
#else
    struct ungappedExtension *goodExtensionBuf = (struct ungappedExtension *)malloc(sizeof(struct ungappedExtension) * goodExtensionBufSize);
    struct alignment *goodAlignBuf = (struct alignment *)malloc(sizeof(struct alignment) * goodAlignBufSize);
    hit_t *binnedHits = (hit_t *)malloc(sizeof(hit_t) * MAX_HITS_PER_SEQ);
#endif
    uint4 numBins = readdb_longestSequenceLength + longestQueryLength + 1;
    numBins = ((numBins + 15) / 16) * 16;
    uint4 *binOffset =
        (uint4 *)malloc(sizeof(uint4) * numBins );
    int maxNumSeqBlkAlign = ((maxNumSeqBlk + 15) / 16) * 16;
    int longestQLenAlign = ((longestQueryLength + 15) / 16) * 16;
    uint4 *numHitsPerQPos = (uint4 *)malloc(sizeof(uint4) * longestQLenAlign);
    hit_t *lastHits = (hit_t *)malloc(sizeof(hit_t) * maxNumSeqBlkAlign);
    uint4 *numExtHit = (uint4 *)malloc(sizeof(uint4) * maxNumSeqBlkAlign);
    //uint4 *allocExtHit = (uint4 *)malloc(sizeof(uint4) * maxNumSeqBlkAlign);
    //memset(allocExtHit, 0,
    //sizeof(uint4) * maxNumSeqBlkAlign );
    hit_t *seqHits = (hit_t *)malloc(sizeof(hit_t) * MAX_NUM_HIT_PER_BLK);
    //hit_t **seqHits = (hit_t **)malloc(sizeof(hit_t *) * maxNumSeqBlkAlign);
    //memset(seqHits, 0,
    //sizeof(hit_t *) * maxNumSeqBlkAlign );
    //int4 *semiGappedScoring_insertQrow_t = (int4 *)malloc(sizeof(int4) * SEMIGAPPED_ROWSIZE * blast_numQuery);
    //int4 *semiGappedScoring_bestRow_t = (int4 *)malloc(sizeof(int4) * SEMIGAPPED_ROWSIZE * blast_numQuery);
    //for(ii = 0; ii < blast_numQuery; ii++)
    //{
        //semiGappedScoring_rowSizes_multi[ii] = SEMIGAPPED_ROWSIZE;
        //semiGappedScoring_insertQrow_multi[ii] = semiGappedScoring_insertQrow_t + ii * SEMIGAPPED_ROWSIZE;
        //semiGappedScoring_bestRow_multi[ii] = semiGappedScoring_bestRow_t + ii * SEMIGAPPED_ROWSIZE;
    //}

    struct gappedExtension *gappedExtension_arr = (struct gappedExtension *)malloc(sizeof(struct gappedExtension) * MAX_NUM_GAPEXT_PER_QUERY * numQueriesPerThread );
    Profile *profiles =
        (Profile *)malloc(sizeof(Profile) );
    memset(profiles, 0, sizeof(Profile) );

    long long start, end;
    start = rdtsc();
    RDTSC_INIT;

    //int4 thread_id = omp_get_thread_num();
    int4 thread_id = 0;
    Profile *profile = profiles + thread_id;
    uint4 *binOffset_p = binOffset + thread_id * numBins;
    hit_t *hits_p = hits + thread_id * MAX_HITS_PER_SEQ;
    uint4 *numHitsPerQPos_p = numHitsPerQPos + thread_id * longestQLenAlign;
    hit_t *binnedHits_p = binnedHits + thread_id * MAX_HITS_PER_SEQ;
    hit_t *lastHits_p = lastHits + thread_id * maxNumSeqBlkAlign;
    uint4 *numExtHit_p = numExtHit + thread_id * maxNumSeqBlkAlign;
    //uint4 *allocExtHit_p = allocExtHit + thread_id * maxNumSeqBlkAlign;
    hit_t *seqHits_p = seqHits + MAX_NUM_HIT_PER_BLK * thread_id;
    struct ungappedExtension **ungappedExtension_sort_t = ungappedExtension_sort + MAX_NUM_TRIGGER_EXT * thread_id; 
    struct ungappedExtension **ungappedExtension_new_t = ungappedExtension_new + MAX_NUM_TRIGGER_EXT * thread_id; 
    BlastGapDP *dp_mem_t = dp_mem + DP_MEM_SIZE * thread_id;
    BlastHSP *BlastHSP_arr_t = BlastHSP_arr + MAX_NUM_HSP * thread_id;

    struct alignment *goodAlignBuf_p = goodAlignBuf + thread_id * MAX_ALIGNMENTS_PER_QUERY * numQueriesPerThread;
    struct ungappedExtension *goodExtensionBuf_p = goodExtensionBuf + thread_id * MAX_EXTENSIONS_PER_QUERY * numQueriesPerThread;
    struct gappedExtension *gappedExtension_t = gappedExtension_arr + MAX_NUM_GAPEXT_PER_QUERY * thread_id * numQueriesPerThread;
    int goodAlignCount = 0;
    int goodExtensionCount = 0;
    int gappedExtension_cnt = 0;

    uint4 blast_numHits_multi_t;
    uint4 blast_numUngappedExtensions_multi_t;
    uint4 blast_numTriggerExtensions_multi_t;
    uint4 blast_numTriggerSequences_multi_t;

    int kk;

    for (ii = 0; ii < proteinLookup_numBlocks; ii++) 
    {
        for (kk = 0; kk < numQuery; kk++) {
            blast_numHits_multi_t = 0;
            blast_numUngappedExtensions_multi_t = 0;
            blast_numTriggerExtensions_multi_t = 0;
            blast_numTriggerSequences_multi_t = 0;

#ifndef COMPRESS_INDEX
            search_protein2hit_lookup_db(
#else
            search_protein2hit_lookup_db_cp(
#endif
                    thread_id, PSSMatrix_arr, scoreMatrix, kk, readdb_sequenceData,
                    ii, binOffset_p, hits_p, numHitsPerQPos_p, binnedHits_p,
                    lastHits_p, seqHits_p, numExtHit_p, 
                    goodAlignBuf_p, &goodAlignCount,
                    goodExtensionBuf_p, &goodExtensionCount,
                    ungappedExtension_new_t, dp_mem_t, tree, private_tree, BlastHSP_arr_t,
                    &blast_numHits_multi_t, &blast_numUngappedExtensions_multi_t,
                    &blast_numTriggerExtensions_multi_t,
                    &blast_numTriggerSequences_multi_t, profile);

            blast_numHits_multi[kk] += blast_numHits_multi_t;
            blast_numUngappedExtensions_multi[kk] += blast_numUngappedExtensions_multi_t;
            blast_numTriggerExtensions_multi[kk] += blast_numTriggerExtensions_multi_t;
            blast_numTriggerSequences_multi[kk] += blast_numTriggerSequences_multi_t;

            if(goodAlignCount >= MAX_ALIGNMENTS_PER_QUERY * numQueriesPerThread)
            {
                printf("ERROR! goodAlignCount = %d >= MAX_ALIGNMENTS_PER_QUERY * numQueriesPerThread\n", goodAlignCount);
                exit(0);
            }
            if(goodExtensionCount >= MAX_EXTENSIONS_PER_QUERY * numQueriesPerThread)
            {
                printf("ERROR! goodExtensionCount = %d >= MAX_EXTENSIONS_PER_QUERY * numQueriesPerThread\n", goodExtensionCount);
                exit(0);
            }

        }

    }


#if 0
    for(ii = 0; ii < goodAlignCount; ii++)
    {

        alignments_findGoodAlignments_ncbi_multi3(
                &goodAlignBuf_p[ii],
                PSSMatrix_arr[goodAlignBuf_p[ii].queryCount], scoreMatrix,
                goodAlignBuf_p[ii].queryCount, ungappedExtension_new, dp_mem,
                tree, private_tree, BlastHSP_arr_t);
    }
#endif



#ifndef NO_STAGE4
    for (jj = 0; jj < numQuery; jj++) {
        RDTSC_START;
        alignments_getTracebacks_ncbi_multi(PSSMatrix_arr[jj], scoreMatrix, &query_composition_arr[jj], jj, NRrecord, querySeq_arr[jj], &matrixInfo, gap_align, matrix, seq_data, newton_system, grads, Scores, z, resids_x, resids_z, old_scores, workspace, gappedExtension_t, &gappedExtension_cnt);
        RDTSC_END;
        GET_RDTSC_ATOMIC(&profile->blast_getTrackCycle);

        if(gappedExtension_cnt >= MAX_NUM_GAPEXT_PER_QUERY * numQueriesPerThread)
        {
            fprintf(stderr, "gappedExtension_arr overflow\n");
            exit(1);
        }
    }
#endif
    // printf("threadid) blast_pseudoCycle, blast_hitDetectCycle,
    // blast_sortCycle, blast_ungappedExtCycle\n");
    end = rdtsc();
#ifdef PROFILE
    fprintf(stderr, "%d) %lld %lld %lld %lld %lld %lld %lld %lld %lld\n", thread_id,
            profile->blast_pseudoCycle, profile->blast_hitDetectCycle,
            profile->blast_sortCycle, profile->blast_ungappedExtCycle,
            profile->blast_findGoodCycle, profile->blast_findFinalCycle,
            profile->blast_getFinalCycle, profile->blast_getTrackCycle, end - start);
#endif

    gettimeofday(&end_t, NULL);
    long hit_time = ((end_t.tv_sec * 1000000 + end_t.tv_usec) -
            (start_t.tv_sec * 1000000 + start_t.tv_usec));
    //#ifdef PROFILE
    fprintf(stderr, "Hit time: %f\n", (float)hit_time * 1e-6);
    //#endif

    free(profiles);

    for (ii = 0; ii < numQuery; ii++) {


#if 0
        printf("Query %d alignment info:\n", ii);
        printf("Query length: %d\n", PSSMatrix_arr[ii].length);
        printf("Number of hits: %d\n", blast_numHits_multi[ii]);
        printf("Number of extensions: %u\n", blast_numUngappedExtensions_multi[ii]);
        printf("Number of successful extensions: %u\n",
                blast_numTriggerExtensions_multi[ii]);
        printf("Number of sequences with successful extensions: %u\n",
                blast_numTriggerSequences_multi[ii]);
        printf("Number of sequences with semi-gapped score above cutoff: %u\n",
                blast_numGoodAlignments_multi[ii]);
        printf("Number of sequences better than %g: %u\n\n", parameters_cutoff,
                alignments_finalAlignments_multi[ii]->numEntries);
        //printf("query: %d numHits: %d numExts: %d numTriggerExts: %d numTriggerSeqs: %d numGapExts: %d numGoodExts: %d numSeqPassed: %d\n", ii, blast_numHits_multi[ii], blast_numUngappedExtensions_multi[ii], blast_numTriggerExtensions_multi[ii], blast_numTriggerSequences_multi[ii], blast_numGappedExtension_multi[ii], blast_numGoodExtensions_multi[ii], blast_numGoodAlignments_multi[ii]);
#else
        printf("Query= %s\n\n", queryDescription_arr[ii]);
        printf("Length=%ld\n\n", strlen(query_arr[ii]));
        print_gappedAlignmentsBrief_multi(ii);
        print_gappedAlignmentsFull_multi(query_arr[ii], PSSMatrix_arr[ii], ii);
#endif



        fasterGappedExtension_free_multi(ii);
        gappedExtension_free_multi(ii);
        alignments_free_multi2(ii);
    }

        free(hits);

#if defined(__ICC) || defined(__INTEL_COMPILER)
        _mm_free(binnedHits);
#else
        free(binnedHits);
        free(binOffset);
        free(numHitsPerQPos);
        free(lastHits);
        free(numExtHit);
        //free(allocExtHit);
        //free(semiGappedScoring_insertQrow_t);
        //free(semiGappedScoring_bestRow_t);
        free(ungappedExtension_new);
        free(ungappedExtension_sort);
        free(dp_mem);
        free(BlastHSP_arr);

        //for(ii = 0; ii < parameters_num_threads; ii++)
        free(workspace);
        free(old_scores);
        free(z);
        free(resids_z);
        free(resids_x);
        Nlm_DenseMatrixFree(&grads);
        Nlm_DenseMatrixFree(&Scores);
        ReNewtonSystemFree(&newton_system);
        free(seq_data);
        _PSIDeallocateMatrix((void **)matrix, BLASTAA_SIZE);
        BLAST_GapAlignStructFree(gap_align);
        Blast_CompositionWorkspaceFree(&NRrecord);
        Blast_MatrixInfoFree2(&matrixInfo);
        Blast_IntervalTreeFree2(private_tree);
        Blast_IntervalTreeFree2(tree);

        for (jj = 0; jj < numQuery; jj++)
        {
            free(querySeq_arr[jj]);
        }
        free(query_composition_arr);

        free(querySeq_arr);

        free(seqHits);
#endif

#if defined(__ICC) || defined(__INTEL_COMPILER)
    _mm_free(goodExtensionBuf);
    _mm_free(goodAlignBuf);
#else
    free(goodExtensionBuf);
    free(goodAlignBuf);
#endif
    free(gappedExtension_arr);

    
}

void alignments_db_omp(struct PSSMatrix *PSSMatrix_arr,
        struct scoreMatrix scoreMatrix, char **query_arr, char **queryDescription_arr, int numQuery) {
    int ii, jj;

    struct timeval start_t, end_t;
    gettimeofday(&start_t, NULL);
    // hit_t *seqHit_arr = (hit_t *)malloc(sizeof(hit_t) * parameters_num_threads
    // * EXT_PER_SEQ * maxNumSeqBlk);
    hit_t *hits = (hit_t *)malloc(sizeof(hit_t) * MAX_HITS_PER_SEQ *
            parameters_num_threads);

    struct ungappedExtension **ungappedExtension_sort, **ungappedExtension_new;
    ungappedExtension_sort = (struct ungappedExtension **)malloc(sizeof(struct ungappedExtension *) * MAX_NUM_TRIGGER_EXT * parameters_num_threads);
    ungappedExtension_new = (struct ungappedExtension **)malloc(sizeof(struct ungappedExtension *) * MAX_NUM_TRIGGER_EXT * parameters_num_threads);
    BlastGapDP *dp_mem = (BlastGapDP *)malloc(sizeof(BlastGapDP) * DP_MEM_SIZE * parameters_num_threads);

    BlastHSP *BlastHSP_arr = (BlastHSP *)malloc(sizeof(BlastHSP) * MAX_NUM_HSP * parameters_num_threads);

    BlastIntervalTree **tree = (BlastIntervalTree **)malloc(sizeof(BlastIntervalTree *) * parameters_num_threads);

    BlastIntervalTree **private_tree = (BlastIntervalTree **)malloc(sizeof(BlastIntervalTree *) * parameters_num_threads);

    Blast_MatrixInfo *matrixInfo_arr = (Blast_MatrixInfo *)malloc(sizeof(Blast_MatrixInfo) * parameters_num_threads);

    Blast_CompositionWorkspace **NRrecord = (Blast_CompositionWorkspace **)malloc(sizeof(Blast_CompositionWorkspace *) * parameters_num_threads); ;

    BlastGapAlignStruct **gap_align_arr = (BlastGapAlignStruct **)malloc(sizeof(BlastGapAlignStruct *) * parameters_num_threads);

    void **matrix_arr = (void **)malloc(sizeof(void *) * parameters_num_threads);

    unsigned char **seq_data_arr = (unsigned char **)malloc(sizeof(unsigned char *) * parameters_num_threads);

    ReNewtonSystem **newton_system_arr = (ReNewtonSystem **)malloc(sizeof(ReNewtonSystem *) * parameters_num_threads);

    double ***grads_arr = (double ***)malloc(sizeof(double **) * parameters_num_threads);
    double ***Scores_arr = (double ***)malloc(sizeof(double **) * parameters_num_threads);

    //ReNewtonSystem *newton_system;
    double **z_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    double **resids_x_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    double **resids_z_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    double **old_scores_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    double **workspace_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    //double **grads;
    //double **Scores;




    int alphsize = COMPO_NUM_TRUE_AA;
    int n  = alphsize * alphsize;
    int mA = 2 * alphsize - 1;
    int m  = mA + 1;

    for(ii = 0; ii < parameters_num_threads; ii++)
    {
        NRrecord[ii] = Blast_CompositionWorkspaceNew();
        tree[ii] = Blast_IntervalTreeInit(0, longestQueryLength + 1,
                0, SEMIGAPPED_ROWSIZE + 1);
        private_tree[ii] = Blast_IntervalTreeInit(0, 
                longestQueryLength + 1, 0, SEMIGAPPED_ROWSIZE + 1);
        Blast_MatrixInfoNew2(&matrixInfo_arr[ii], BLASTAA_SIZE, BLASTAA_SIZE, FALSE);

        BLAST_GapAlignStructNew(readdb_longestSequenceLength, &gap_align_arr[ii]);

        matrix_arr[ii] = (void *)_PSIAllocateMatrix(BLASTAA_SIZE, BLASTAA_SIZE, sizeof(Int4)); 

        seq_data_arr[ii] = (unsigned char *)malloc((readdb_longestSequenceLength + 10) * sizeof(unsigned char));

        newton_system_arr[ii] = ReNewtonSystemNew(alphsize); 

        grads_arr[ii] = Nlm_DenseMatrixNew(2, n);
        Scores_arr[ii] = Nlm_DenseMatrixNew(BLASTAA_SIZE, BLASTAA_SIZE);

        resids_x_arr[ii] = (double *) malloc(n * sizeof(double));
        resids_z_arr[ii] = (double *) malloc((mA + 1) * sizeof(double));
        z_arr[ii] = (double *) calloc( mA + 1,   sizeof(double));
        old_scores_arr[ii] = (double *) malloc(n * sizeof(double));
        workspace_arr[ii] = (double *) malloc(n * sizeof(double));

    }



    char **querySeq_arr = (char **)malloc(sizeof(char *) * numQuery);


    Blast_AminoAcidComposition *query_composition_arr = (Blast_AminoAcidComposition *)malloc(sizeof(Blast_AminoAcidComposition) * numQuery);

    for(ii = 0; ii < numQuery; ii++)
    {
        querySeq_arr[ii] = (char *)malloc(sizeof(char) * PSSMatrix_arr[ii].length);
        for(jj = 0; jj < PSSMatrix_arr[ii].length; jj++)
        {
            ASSERT(PSSMatrix_arr[ii].queryCodes[jj] < FSA_AA_SIZE);
            querySeq_arr[ii][jj] = to_ncbi[PSSMatrix_arr[ii].queryCodes[jj]];
        }
        Blast_ReadAaComposition_fsa(&query_composition_arr[ii], BLASTAA_SIZE, PSSMatrix_arr[ii].queryCodes, PSSMatrix_arr[ii].length);
    }

    int numQueriesPerThread = numQuery / parameters_num_threads;
    long long goodAlignBufSize = MAX_ALIGNMENTS_PER_QUERY * numQueriesPerThread * parameters_num_threads;
    int goodExtensionBufSize = MAX_EXTENSIONS_PER_QUERY * numQueriesPerThread * parameters_num_threads;
#if defined(__ICC) || defined(__INTEL_COMPILER)
    hit_t *binnedHits = (hit_t *)_mm_malloc(
            sizeof(hit_t) * MAX_HITS_PER_SEQ * parameters_num_threads, 64);
    struct alignment *goodAlignBuf = (struct alignment *)_mm_malloc(sizeof(struct alignment) * goodAlignBufSize, 64);
    struct ungappedExtension *goodExtensionBuf = (struct ungappedExtension *)_mm_malloc(sizeof(struct ungappedExtension) * goodExtensionBufSize, 64);
#else
    struct ungappedExtension *goodExtensionBuf = (struct ungappedExtension *)malloc(sizeof(struct ungappedExtension) * goodExtensionBufSize);
    struct alignment *goodAlignBuf = (struct alignment *)malloc(sizeof(struct alignment) * goodAlignBufSize);
    hit_t *binnedHits = (hit_t *)malloc(sizeof(hit_t) * MAX_HITS_PER_SEQ *
            parameters_num_threads);
#endif
    uint4 numBins = readdb_longestSequenceLength + longestQueryLength + 1;
    numBins = ((numBins + 15) / 16) * 16;
    uint4 *binOffset =
        (uint4 *)malloc(sizeof(uint4) * numBins * parameters_num_threads);
    int maxNumSeqBlkAlign = ((maxNumSeqBlk + 15) / 16) * 16;
    int longestQLenAlign = ((longestQueryLength + 15) / 16) * 16;
    uint4 *numHitsPerQPos = (uint4 *)malloc(sizeof(uint4) * longestQLenAlign *
            parameters_num_threads);
    hit_t *lastHits = (hit_t *)malloc(sizeof(hit_t) * maxNumSeqBlkAlign *
            parameters_num_threads);
    uint4 *numExtHit = (uint4 *)malloc(sizeof(uint4) * maxNumSeqBlkAlign *
            parameters_num_threads);
    //uint4 *allocExtHit = (uint4 *)malloc(sizeof(uint4) * maxNumSeqBlkAlign *
    //parameters_num_threads);
    //memset(allocExtHit, 0,
    //sizeof(uint4) * maxNumSeqBlkAlign * parameters_num_threads);
    hit_t *seqHits = (hit_t *)malloc(sizeof(hit_t) * MAX_NUM_HIT_PER_BLK * parameters_num_threads);
    //int4 *semiGappedScoring_insertQrow_t = (int4 *)malloc(sizeof(int4) * SEMIGAPPED_ROWSIZE * blast_numQuery);
    //int4 *semiGappedScoring_bestRow_t = (int4 *)malloc(sizeof(int4) * SEMIGAPPED_ROWSIZE * blast_numQuery);
    //for(ii = 0; ii < blast_numQuery; ii++)
    //{
        //semiGappedScoring_rowSizes_multi[ii] = SEMIGAPPED_ROWSIZE;
        //semiGappedScoring_insertQrow_multi[ii] = semiGappedScoring_insertQrow_t + ii * SEMIGAPPED_ROWSIZE;
        //semiGappedScoring_bestRow_multi[ii] = semiGappedScoring_bestRow_t + ii * SEMIGAPPED_ROWSIZE;
    //}

    struct gappedExtension *gappedExtension_arr = (struct gappedExtension *)malloc(sizeof(struct gappedExtension) * MAX_NUM_GAPEXT_PER_QUERY * numQueriesPerThread * parameters_num_threads);
    Profile *profiles =
        (Profile *)malloc(sizeof(Profile) * parameters_num_threads);
    memset(profiles, 0, sizeof(Profile) * parameters_num_threads);

    //#pragma omp parallel num_threads(parameters_num_threads)
    //{
    // int4 thread_id = omp_get_thread_num();
    // printf("dbblast_omp: affinity - thread_id %d bind to core_id: %d\n",
    // thread_id, sched_getcpu());
    //}
    long long start, end;
    start = rdtsc();

#pragma omp parallel num_threads(parameters_num_threads) default(none) shared( \
        proteinLookup_numBlocks, numQuery, PSSMatrix_arr, readdb_sequenceData,     \
        numBins, maxNumSeqBlkAlign, longestQLenAlign, binOffset, hits,             \
        numHitsPerQPos, binnedHits, lastHits, numExtHit, seqHits,     \
        blast_numHits_multi, blast_numUngappedExtensions_multi,                    \
        blast_numTriggerExtensions_multi, blast_numTriggerSequences_multi,         \
        profiles, scoreMatrix, blast_numGoodAlignments_multi, stderr, \
        ungappedExtension_sort, ungappedExtension_new, dp_mem, tree, \
        private_tree, query_composition_arr, NRrecord, querySeq_arr,\
        BlastHSP_arr, matrixInfo_arr, gap_align_arr, matrix_arr, seq_data_arr,\
        newton_system_arr, grads_arr, Scores_arr, z_arr, resids_x_arr, resids_z_arr,\
        old_scores_arr, workspace_arr, goodAlignBuf, numQueriesPerThread, goodExtensionBuf, gappedExtension_arr, readdb_longestSequenceLength)

    {
        RDTSC_INIT;

        int4 thread_id = omp_get_thread_num();
        Profile *profile = profiles + thread_id;
        uint4 *binOffset_p = binOffset + thread_id * numBins;
        hit_t *hits_p = hits + thread_id * MAX_HITS_PER_SEQ;
        uint4 *numHitsPerQPos_p = numHitsPerQPos + thread_id * longestQLenAlign;
        hit_t *binnedHits_p = binnedHits + thread_id * MAX_HITS_PER_SEQ;
        hit_t *lastHits_p = lastHits + thread_id * maxNumSeqBlkAlign;
        uint4 *numExtHit_p = numExtHit + thread_id * maxNumSeqBlkAlign;
        //uint4 *allocExtHit_p = allocExtHit + thread_id * maxNumSeqBlkAlign;
        hit_t *seqHits_p = seqHits + (long)MAX_NUM_HIT_PER_BLK * thread_id;
        struct ungappedExtension **ungappedExtension_sort_t = ungappedExtension_sort + MAX_NUM_TRIGGER_EXT * thread_id; 
        struct ungappedExtension **ungappedExtension_new_t = ungappedExtension_new + MAX_NUM_TRIGGER_EXT * thread_id; 
        BlastGapDP *dp_mem_t = dp_mem + DP_MEM_SIZE * thread_id;
        BlastHSP *BlastHSP_arr_t = BlastHSP_arr + MAX_NUM_HSP * thread_id;

        struct alignment *goodAlignBuf_p = goodAlignBuf + thread_id * MAX_ALIGNMENTS_PER_QUERY * numQueriesPerThread;
        struct ungappedExtension *goodExtensionBuf_p = goodExtensionBuf + thread_id * MAX_EXTENSIONS_PER_QUERY * numQueriesPerThread;
        struct gappedExtension *gappedExtension_t = gappedExtension_arr + MAX_NUM_GAPEXT_PER_QUERY * thread_id * numQueriesPerThread;
        int goodAlignCount = 0;
        int goodExtensionCount = 0;
        int gappedExtension_cnt = 0;
        uint4 blast_numHits_multi_t;
        uint4 blast_numUngappedExtensions_multi_t;
        uint4 blast_numTriggerExtensions_multi_t;
        uint4 blast_numTriggerSequences_multi_t;

        int ii;
        //#pragma omp for schedule(dynamic)
        for (ii = 0; ii < proteinLookup_numBlocks; ii++) {
            int kk;
#pragma omp for schedule(dynamic) nowait
            //#pragma omp for schedule(dynamic)
            for (kk = 0; kk < numQuery; kk++) {
                blast_numHits_multi_t = 0;
                blast_numUngappedExtensions_multi_t = 0;
                blast_numTriggerExtensions_multi_t = 0;
                blast_numTriggerSequences_multi_t = 0;

                // printf("tid: %d queryId: %d blockId: %d\n", thread_id, kk, ii);
                // for(ii = 0; ii < proteinLookup_numBlocks; ii++)
                {

#ifndef COMPRESS_INDEX
                    search_protein2hit_lookup_db(
#else
                    search_protein2hit_lookup_db_cp(
#endif
                            thread_id, PSSMatrix_arr, scoreMatrix, kk, readdb_sequenceData,
                            ii, binOffset_p, hits_p, numHitsPerQPos_p, binnedHits_p,
                            lastHits_p, seqHits_p, numExtHit_p, 
                            goodAlignBuf_p, &goodAlignCount,
                            goodExtensionBuf_p, &goodExtensionCount,
                            ungappedExtension_new_t, dp_mem_t, tree[thread_id], private_tree[thread_id], BlastHSP_arr_t,
                            &blast_numHits_multi_t, &blast_numUngappedExtensions_multi_t,
                            &blast_numTriggerExtensions_multi_t,
                            &blast_numTriggerSequences_multi_t, profile);
                }
                if(goodAlignCount >= MAX_ALIGNMENTS_PER_QUERY * numQueriesPerThread)
                {
                    printf("ERROR! goodAlignCount = %d >= MAX_ALIGNMENTS_PER_QUERY * numQueriesPerThread\n", goodAlignCount);
                    exit(0);
                }
                if(goodExtensionCount >= MAX_EXTENSIONS_PER_QUERY * numQueriesPerThread)
                {
                    fprintf(stderr, "ERROR! goodExtensionCount = %d >= MAX_EXTENSIONS_PER_QUERY * numQueriesPerThread\n", goodExtensionCount);
                    exit(0);
                }

                __sync_fetch_and_add(blast_numHits_multi + kk, blast_numHits_multi_t);
                __sync_fetch_and_add(blast_numUngappedExtensions_multi + kk,
                        blast_numUngappedExtensions_multi_t);
                __sync_fetch_and_add(blast_numTriggerExtensions_multi + kk,
                        blast_numTriggerExtensions_multi_t);
                __sync_fetch_and_add(blast_numTriggerSequences_multi + kk,
                        blast_numTriggerSequences_multi_t);
            }
        }

#if 0
        for(ii = 0; ii < goodAlignCount; ii++)
        {

            alignments_findGoodAlignments_ncbi_multi3(
                    &goodAlignBuf_p[ii],
                    PSSMatrix_arr[goodAlignBuf_p[ii].queryCount], scoreMatrix,
                    goodAlignBuf_p[ii].queryCount, ungappedExtension_new, dp_mem,
                    tree[thread_id], private_tree[thread_id], BlastHSP_arr_t);
        }
#endif

#ifndef NO_STAGE4
#pragma omp for schedule(dynamic)
        for (jj = 0; jj < numQuery; jj++) {

            RDTSC_START;
            alignments_getTracebacks_ncbi_multi(PSSMatrix_arr[jj], scoreMatrix, &query_composition_arr[jj], jj, NRrecord[thread_id], querySeq_arr[jj], matrixInfo_arr + thread_id, gap_align_arr[thread_id], matrix_arr[thread_id], seq_data_arr[thread_id], newton_system_arr[thread_id], grads_arr[thread_id], Scores_arr[thread_id], z_arr[thread_id], resids_x_arr[thread_id], resids_z_arr[thread_id], old_scores_arr[thread_id], workspace_arr[thread_id], gappedExtension_t, &gappedExtension_cnt);
            RDTSC_END;
            GET_RDTSC_ATOMIC(&profile->blast_getTrackCycle);

            if(gappedExtension_cnt >= MAX_NUM_GAPEXT_PER_QUERY * numQueriesPerThread)
            {
                fprintf(stderr, "gappedExtension_arr overflow\n");
                exit(1);
            }
            ////printf("dbblast_omp: thread %d - query %d finished\n", thread_id, jj);
        }
#endif


#if defined(PROFILE) && defined(THREAD_PROFILE)
        // printf("threadid) blast_pseudoCycle, blast_hitDetectCycle,
        // blast_sortCycle, blast_ungappedExtCycle\n");
        fprintf(stderr, "%d) %lld %lld %lld %lld %lld %lld %lld %lld\n", thread_id,
                profile->blast_pseudoCycle, profile->blast_hitDetectCycle,
                profile->blast_sortCycle, profile->blast_ungappedExtCycle,
                profile->blast_findGoodCycle, profile->blast_findFinalCycle,
                profile->blast_getFinalCycle, profile->blast_getTrackCycle);
#endif
    }

    end = rdtsc();

    long long blast_pseudoCycle = 0, blast_hitDetectCycle = 0,
         blast_sortCycle = 0, blast_ungappedExtCycle = 0,
         blast_findGoodCycle = 0, blast_findFinalCycle = 0,
         blast_getFinalCycle = 0, blast_getTrackCycle = 0;
    for (ii = 0; ii < parameters_num_threads; ii++) {
        blast_pseudoCycle += profiles[ii].blast_pseudoCycle;
        blast_hitDetectCycle += profiles[ii].blast_hitDetectCycle;
        blast_sortCycle += profiles[ii].blast_sortCycle;
        blast_ungappedExtCycle += profiles[ii].blast_ungappedExtCycle;
        blast_findGoodCycle += profiles[ii].blast_findGoodCycle;
        blast_findFinalCycle += profiles[ii].blast_findFinalCycle;
        blast_getFinalCycle += profiles[ii].blast_getFinalCycle;
        blast_getTrackCycle += profiles[ii].blast_getTrackCycle;
    }

#ifdef PROFILE
    fprintf(stderr, "all) %lld %lld %lld %lld %lld %lld %lld %lld %lld\n",
            blast_pseudoCycle, blast_hitDetectCycle, blast_sortCycle,
            blast_ungappedExtCycle, blast_findGoodCycle, blast_findFinalCycle,
            blast_getFinalCycle, blast_getTrackCycle, end - start);
#endif

    gettimeofday(&end_t, NULL);
    long hit_time = ((end_t.tv_sec * 1000000 + end_t.tv_usec) -
            (start_t.tv_sec * 1000000 + start_t.tv_usec));
    fprintf(stderr, "Hit time: %f\n", (float)hit_time * 1e-6);

    free(profiles);

    free(hits);

#if defined(__ICC) || defined(__INTEL_COMPILER)
    _mm_free(binnedHits);
#else
    free(binnedHits);
#endif
    free(binOffset);
    free(numHitsPerQPos);
    free(lastHits);
    free(numExtHit);
    //free(allocExtHit);
    //free(semiGappedScoring_insertQrow_t);
    //free(semiGappedScoring_bestRow_t);
    free(ungappedExtension_new);
    free(ungappedExtension_sort);
    free(dp_mem);
    free(BlastHSP_arr);



    for(ii = 0; ii < parameters_num_threads; ii++)
    {

        free(workspace_arr[ii]);
        free(old_scores_arr[ii]);
        free(z_arr[ii]);
        free(resids_z_arr[ii]);
        free(resids_x_arr[ii]);
        Nlm_DenseMatrixFree(&grads_arr[ii]);
        Nlm_DenseMatrixFree(&Scores_arr[ii]);
        ReNewtonSystemFree(&newton_system_arr[ii]);
        free(seq_data_arr[ii]);
        _PSIDeallocateMatrix((void **)matrix_arr[ii], BLASTAA_SIZE);
        BLAST_GapAlignStructFree(gap_align_arr[ii]);
        Blast_CompositionWorkspaceFree(&NRrecord[ii]);
        Blast_MatrixInfoFree2(&matrixInfo_arr[ii]);
        Blast_IntervalTreeFree2(private_tree[ii]);
        Blast_IntervalTreeFree2(tree[ii]);
    }

    free(workspace_arr);
    free(old_scores_arr);
    free(z_arr);
    free(resids_z_arr);
    free(resids_x_arr);
    free(grads_arr);
    free(Scores_arr);
    free(newton_system_arr);
    free(gap_align_arr);
    free(matrix_arr);
    free(seq_data_arr);

    free(NRrecord);

    free(matrixInfo_arr);

    free(tree);
    free(private_tree);

    for (jj = 0; jj < numQuery; jj++)
    {
        free(querySeq_arr[jj]);
    }
    free(query_composition_arr);

    free(querySeq_arr);

    //for (jj = 0; jj < maxNumSeqBlkAlign * parameters_num_threads; jj++) {
    //free(seqHits[jj]);
    //}

    free(seqHits);

    for (ii = 0; ii < numQuery; ii++) {
        printf("Query= %s\n\n", queryDescription_arr[ii]);
        printf("Length=%ld\n\n", strlen(query_arr[ii]));

        //print_info(ii, PSSMatrix_arr);
        print_gappedAlignmentsBrief_multi(ii);
        print_gappedAlignmentsFull_multi(query_arr[ii], PSSMatrix_arr[ii], ii);

#if 0
        printf("Query %d alignment info:\n", ii);
        printf("Query length: %d\n", PSSMatrix_arr[ii].length);
        printf("Number of hits: %d\n", blast_numHits_multi[ii]);
        printf("Number of extensions: %u\n", blast_numUngappedExtensions_multi[ii]);
        printf("Number of successful extensions: %u\n",
                blast_numTriggerExtensions_multi[ii]);
        printf("Number of sequences with successful extensions: %u\n",
                blast_numTriggerSequences_multi[ii]);
        printf("Number of sequences with semi-gapped score above cutoff: %u\n",
                blast_numGoodAlignments_multi[ii]);
        printf("Number of sequences better than %g: %u\n\n", parameters_cutoff,
                alignments_finalAlignments_multi[ii]->numEntries);
        printf("%d numHits: %d numExts: %d numTriggerExts: %d numGapExts: %d numGoodExts: %d numSeqPassed: %d\n", ii, blast_numHits_multi[ii], blast_numUngappedExtensions_multi[ii], blast_numTriggerExtensions_multi[ii], blast_numGappedExtension_multi[ii], blast_numGoodExtensions_multi[ii], blast_numGoodAlignments_multi[ii]);
#endif

        fasterGappedExtension_free_multi(ii);
        gappedExtension_free_multi(ii);
        //PSSMatrix_free(PSSMatrix_arr[ii]);
        alignments_free_multi2(ii);
    }


#if defined(__ICC) || defined(__INTEL_COMPILER)
    _mm_free(goodExtensionBuf);
    _mm_free(goodAlignBuf);
#else
    free(goodExtensionBuf);
    free(goodAlignBuf);
#endif
    free(gappedExtension_arr);
}

// Perform regular gapped scoring on an extension
void alignments_regularGappedAlignment_multi(
        struct PSSMatrix PSSMatrix, struct ungappedExtension *ungappedExtension,
        struct alignment *alignment, int queryNum) {
    int4 newScore;
    unsigned char *subject;

    subject = alignments_selectRegion(alignment, ungappedExtension);

    //    printf("[%d] subject=%p", subject[ungappedExtension->seed.subjectOffset
    // / 4], subject); fflush(stdout);

    blast_dloc_multi[queryNum] = alignment->descriptionLocation;

    newScore = gappedScoring_score_multi(
            ungappedExtension, PSSMatrix, alignment->subjectLength, subject,
            statistics_gappedNominalDropoff_multi[queryNum], queryNum);

#ifdef VERBOSE
    if (blast_dloc_multi[queryNum] == parameters_verboseDloc)
        printf("Was %d Now %d\n", ungappedExtension->nominalScore, newScore);
#endif

    // If the extension's score has dropped considerably
    if (newScore * parameters_semiGappedR2 < ungappedExtension->nominalScore) {
        // Rescore with larger dropoff
        newScore = gappedScoring_score_multi(
                ungappedExtension, PSSMatrix, alignment->subjectLength, subject,
                statistics_gappedFinalNominalDropoff_multi[queryNum], queryNum);

#ifdef VERBOSE
        if (blast_dloc_multi[queryNum] == parameters_verboseDloc)
            printf("Rescore now %d\n", newScore);
#endif
    }

    ungappedExtension->nominalScore = newScore;

    if (ungappedExtension->nominalScore <
            blast_gappedNominalCutoff_multi[queryNum]) {
        // If extension scores below cutoff, mark it as deleted
        ungappedExtension->status = ungappedExtension_DELETED;
    } else {
        // Else mark it as scored with regular gapped alignment
        ungappedExtension->status = ungappedExtension_GAPPED;
    }

    // Remove any ungapped extensions in this alignment that are in the area
    // covered
    // by the gapped scoring just performed
    alignments_pruneRegion_multi(alignment, ungappedExtension, queryNum);
}

// Check for ungapped extensions that were mistakenly pruned. If one is found,
// add a non-deleted
// copy to the end of list of ungapped extensions
void
alignments_unpruneRegion_multi(struct alignment *alignment,
        struct ungappedExtension *oldUngappedExtension,
        struct ungappedExtension *ungappedExtension,
        struct PSSMatrix PSSMatrix, int queryNum) {
    struct ungappedExtension *currentExtension, *previousExtension,
    *restoreExtension;

#ifdef VERBOSE
    if (alignment->descriptionLocation == parameters_verboseDloc) {
        printf("unpruneRegion dloc=%d oldRegion:", alignment->descriptionLocation);
        printf("Pruning old: ");
        ungappedExtension_print(oldUngappedExtension);
        printf("                      newRegion:");
        printf("Pruning new: ");
        ungappedExtension_print(ungappedExtension);
    }
#endif

    currentExtension = alignment->ungappedExtensions;
    previousExtension = NULL;
    while (currentExtension != NULL) {
        // If old table-driven alignment contains extension but new gapped alignment
        // does not, undelete it
        if (currentExtension->status == ungappedExtension_DELETED &&
                alignments_contains(oldUngappedExtension, currentExtension) &&
                !alignments_contains(ungappedExtension, currentExtension)) {
            restoreExtension = currentExtension;

#ifdef VERBOSE
            if (alignment->descriptionLocation == parameters_verboseDloc) {
                printf("Restoring: ");
                ungappedExtension_print(restoreExtension);
            }
#endif

            // Removing from front of list
            if (previousExtension == NULL) {
                alignment->ungappedExtensions = currentExtension->next;
                currentExtension = alignment->ungappedExtensions;
            }
            // Removing from the middle of the list
            else {
                previousExtension->next = currentExtension->next;
                currentExtension = previousExtension->next;
            }

            // Add the restored extension to end of the list and change status
            alignments_addUngappedExtensionAtEnd(alignment, restoreExtension);
            restoreExtension->status = ungappedExtension_UNGAPPED;

            // Find the seed
            ungappedExtension_findSeed(restoreExtension, PSSMatrix,
                    alignment->subject);
        } else {
            // Advance to next ungapped extension
            previousExtension = currentExtension;
            currentExtension = currentExtension->next;
        }
    }
}

// Check other ungapped extensions that have been scored to see if we can join
// it to this one to get a higher scoring HSP
void alignments_checkForJoin_multi(struct alignment *alignment,
        struct ungappedExtension *extension1,
        struct PSSMatrix PSSMatrix, int queryNum) {
    struct ungappedExtension *extension2;
    int4 queryDistance, subjectDistance, newScore;
    unsigned char *subject;

    subject = alignments_selectRegion(alignment, extension1);

    extension2 = alignment->ungappedExtensions;
    blast_dloc = alignment->descriptionLocation;

    // For each extension that has already been gapped scored
    while (extension2 != NULL) {
        // Check extension2 has been scored, and that the combined scores
        // of extension 1&2 could possibly exceed the cutoff
        if (extension2 != extension1 &&
                extension2->status != ungappedExtension_DELETED) {
            // If extension2 comes after extension1, determine distance between them
            if (extension1->start.queryOffset < extension2->start.queryOffset) {
                queryDistance =
                    extension2->start.queryOffset - extension1->end.queryOffset;
                subjectDistance =
                    extension2->start.subjectOffset - extension1->end.subjectOffset;
            }
            // Else extension1 comes after extension2
            else {
                queryDistance =
                    extension1->start.queryOffset - extension2->end.queryOffset;
                subjectDistance =
                    extension1->start.subjectOffset - extension2->end.subjectOffset;
            }

#ifdef VERBOSE
            if (parameters_verboseDloc == alignment->descriptionLocation) {
                printf("Check for join:\n");
                ungappedExtension_print(extension1);
                ungappedExtension_print(extension2);
            }
#endif

            // Quick check the distance between HSPs is not too great
            if ((queryDistance >= 0 || subjectDistance >= 0) &&
                    minimum(abs(queryDistance), abs(subjectDistance)) +
                    abs(queryDistance - subjectDistance) * parameters_extendGap <=
                    2 * statistics_gappedFinalNominalDropoff_multi[queryNum]) {
                // Perform gappedScoring with higher dropoff to try and bridge the gap
                newScore = gappedScoring_score_multi(
                        extension1, PSSMatrix, alignment->subjectLength, subject,
                        statistics_gappedFinalNominalDropoff_multi[queryNum], queryNum);

#ifdef VERBOSE
                if (parameters_verboseDloc == alignment->descriptionLocation) {
                    printf("Performed join. New score=%d\n", newScore);
                }
#endif

                blast_numAttemptedJoin++;
                // Check if we successfully joined two HSPs
                if (newScore > extension1->nominalScore) {
                    extension1->nominalScore = newScore;
                    extension1->status = ungappedExtension_JOINED;
                    extension2->status = ungappedExtension_JOINED;
                    blast_numSuccessfullyJoined++;
                }

                return;
            }
        }
        extension2 = extension2->next;
    }
}

// Expand a cluster and perform stages 1 and 2 of search on the children
int alignments_expandCluster_multi(struct alignment *alignment,
        struct PSSMatrix PSSMatrix, int queryNum) {
    struct child *children, *child;
    int4 numChildren, childNum, numNewAlignments = 0;
    struct sequenceData sequenceData;

    // Get children
    children = readdb_getChildren(alignment->subject, alignment->subjectLength,
            alignment->encodedLength,
            alignment->descriptionLocation, &numChildren);

    //	printf("Alignment dloc=%d\n", alignment->descriptionLocation);

    if (numChildren == 0)
        return 0;

    // For each child
    childNum = 0;
    while (childNum < numChildren) {
        child = children + childNum;

        sequenceData.descriptionStart = child->descriptionLocation;
        sequenceData.descriptionLength = child->descriptionLength;
        sequenceData.sequenceLength = child->length;
        sequenceData.encodedLength = child->length;
        sequenceData.sequence = child->sequence;

        // Re-initialize the hitMatrix
        hitMatrix_reinitialize(PSSMatrix.length, readdb_longestSequenceLength,
                child->sequence);

        // Perform stage 1 and 2 of search
        if (parameters_oneHitTrigger) {
            search_protein1hit(PSSMatrix, &sequenceData, 1, constants_maxSignedInt);
        } else {
            search_protein2hit(PSSMatrix, &sequenceData, 1, constants_maxSignedInt);
        }

        if (alignments_currentAlignment != NULL) {
            if (numNewAlignments == 0)
                alignments_numClusters++;

            numNewAlignments++;
            alignments_currentAlignment->inMemorySubject = 1;
            alignments_currentAlignment->cluster = alignments_numClusters;
            //            printf("Child %d start=%d length=%d dloc=%d\n", childNum,
            // child->regionStart,
            //                   child->length, child->descriptionLocation);
        }

        blast_numExpandedSequences_multi[queryNum]++;
        childNum++;
    }

    /*    if (numNewAlignments == 0)
          {
          print_singleSequence(alignment->subject, alignment->subjectLength);
          printf("\n");
          }*/

    return 1;
}




