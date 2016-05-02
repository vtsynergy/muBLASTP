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

#include "blast.h"

int get_diag_arr_length(int queryLength)
{
    int diag_arr_length = 1;
    while(diag_arr_length < queryLength + parameters_A)
    {
        diag_arr_length = diag_arr_length << 1;
    }

    return diag_arr_length;
}

#define BLAST_CMP(a,b) ((a)>(b) ? 1 : ((a)<(b) ? -1 : 0))

static int compareUngappedExtQuery(const void *v1, const void *v2)
{
    struct ungappedExtension *h1, *h2;
    int result = 0;

    h1 = ((struct ungappedExtension *)v1);
    h2 = ((struct ungappedExtension *)v2);

    int4 l1 = h1->end.subjectOffset - h1->start.subjectOffset;
    int4 l2 = h2->end.subjectOffset - h2->start.subjectOffset;

    if (    0 == (result = BLAST_CMP(h1->queryCount,
                    h2->queryCount)) &&
            0 == (result = BLAST_CMP(h2->nominalScore,
                    h1->nominalScore)) &&
            0 == (result = BLAST_CMP(h1->start.subjectOffset,
                    h2->start.subjectOffset)) &&
            0 == (result = BLAST_CMP(l2,
                    l1)) &&
            0 == (result = BLAST_CMP(h1->start.queryOffset,
                    h2->start.queryOffset))) {
        result = BLAST_CMP(l2,
                l1);
    }

    return result;
}


void alignments_sortUngapedExtensionQuery(struct ungappedExtension *ungappedExtensions, int numExts)
{
    qsort(ungappedExtensions,
            numExts,
            sizeof(struct ungappedExtension), compareUngappedExtQuery);
}


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
        )
{

    uint2 subjectOffset, wordLengthMinusOne, count = 0;
    uint2 queryOffset;
    struct ungappedExtension *ungappedExtension;
    int4 diagonal;

    struct initialWord_protein_db *initialWord;
    wordLengthMinusOne = parameters_wordSize - 1;
    unsigned char *ungappedExtension_subjectEndReached_t;

#ifdef PROFILE
    RDTSC_INIT;
#endif

    uint4 sequenceCount = subSequenceId;
    uint4 descriptionStart = 0, descriptionLength = 0, encodedLength;
    descriptionLength = sequenceData[sequenceCount].descriptionLength;
    descriptionStart = sequenceData[sequenceCount].descriptionStart;
    uint4 subjectLength = sequenceData[sequenceCount].sequenceLength;
    //address = subject = sequenceData[sequenceCount].sequence;


    int4 length = sequenceData[subSequenceId].sequenceLength;  

    unsigned char *subject = sequenceData[subSequenceId].sequence;

    int numHits = 0, numUngappedExt = 0, numTriggerExt = 0, numTriggerSeq = 0;
    int maxDiag = longestQueryLength + length;
    //int maxDiag2 = get_diag_arr_length(longestQueryLength + parameters_A);
    //int diagMask = maxDiag2 - 1;

#ifdef PROFILE
    RDTSC_START;
#endif

    //int numSecondBins = ((numQuery + 1) * maxDiag2) + 1;
    int numSecondBins = ((numQuery + 1) * maxDiag) + 1;
    int4 numGoodAlign = *goodAlignCount;
    int4 numGoodExtensions = *goodExtensionCount;
    int ii, jj, kk;

    //int bitCount = 0, bitValue = 1;
    //while(bitValue < numSecondBins)
    //{
        //bitValue *= 2;
        //bitCount++;
    //}

    //bitCount = ceil((float)bitCount / LOG_NUM_BINS) * LOG_NUM_BINS;

    int4 maxNumHitPerSeq = MAX_NUM_HIT_PER_QUERY / parameters_num_threads;  

    if(numSecondBins >= MAX_NUM_SECONDARY_BINS)
    {
        fprintf(stderr, "ERROR! numSecondBins = %d\n", numSecondBins);
        exit(0);
    }

    for(ii = 0; ii < numSecondBins; ii++)
    {
        lastHits[ii] = 0xFFFF;
    }

    memset(numExtHit, 0, sizeof(uint4) * numQuery);

    struct initialWord_protein_query *proteinLookup_query_currBlk =
        proteinLookup_query_blk[block_id];
    PV_ARRAY_TYPE *pv_currBlk = pv_arr[block_id];

    int numExt = 0;
    int seqLen = length - wordLengthMinusOne;
    int seqPos;
    int numUngappedExtSeq = 0;
    for(seqPos = 0; seqPos < seqLen; seqPos++) {

        int codeWord = wordLookupDFA_getCodeword(subject + seqPos, parameters_wordSize);

        //if (PV_TEST(pv_currBlk, codeWord, PV_ARRAY_BTS)) 
        {

            uint32_t *pos_arr = NULL;

            if(proteinLookup_query_currBlk[codeWord].numQueryPositions > 3)
                pos_arr = proteinLookup_query_currBlk[codeWord].querySequencePositions; 
            else
                pos_arr = proteinLookup_query_currBlk[codeWord].embeddedPositions; 

            //if(proteinLookup_query[codeWord].numQueryPositions > 3)
                //pos_arr = proteinLookup_query[codeWord].querySequencePositions; 
            //else
                //pos_arr = proteinLookup_query[codeWord].embeddedPositions; 

            for (jj = 0; jj < proteinLookup_query_currBlk[codeWord].numQueryPositions; jj++) 
            {
                hit_t hit = pos_arr[jj];
                queryOffset = hit >> 16;
                diagonal =
                    queryOffset - seqPos + (length - wordLengthMinusOne);

                int queryNum = hit & 0xFFFF;
                //blast_numHits_multi_t[seqId]++;

                uint4 secondBinId = queryNum * maxDiag + diagonal;

                uint4 currHit = seqPos;
                uint4 lastHit = lastHits[secondBinId];

                int distance = currHit - lastHit;

                if(distance >= parameters_A || lastHit == 0xFFFF)
                {
                    lastHits[secondBinId] = currHit;
                }
                else if(distance >= parameters_overlap)
                {
                    lastHits[secondBinId] = currHit;
                    int subjectOffset = currHit;
                    subjectOffset += wordLengthMinusOne;
                    queryOffset += wordLengthMinusOne;

                    int distance = currHit - lastHit;
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
                            lastHits[secondBinId] = subjectOffset;
                        }
                    }
                }
            }
        }
    }


    alignments_sortUngapedExtensionQuery(goodExtensionBuf + numGoodExtensions, numUngappedExtSeq);

    numGoodExtensions += numUngappedExtSeq;

    //for(ii = 1; ii < numUngappedExtSeq; ii++)
    //{
        //ASSERT(goodExtensionBuf[numGoodExtensions + ii].queryCount >= goodExtensionBuf[numGoodExtensions + ii - 1].queryCount); 
        //ASSERT(goodExtensionBuf[numGoodExtensions + ii].sequenceCount == sequenceCount); 
    //}

    //numGoodExtensions += numUngappedExtSeq;
#if 1
    int prevSeqId = -1;
    struct alignment *alignment = NULL;
    int new_numUngappedExtSeq = 0;
    int numAlignBlk = 0;
    int queryNum = -1;
    for(ii = 0; ii < numUngappedExtSeq; ii++)
    {
        queryNum = goodExtensionBuf[numGoodExtensions + ii].queryCount;

        //fprintf(stderr, "sequenceCount: %d start: %d %d\n",
                //sequenceCount, goodExtensionBuf[numGoodExtensions+ ii].start.queryOffset,
                //goodExtensionBuf[numGoodExtensions+ ii].start.subjectOffset);

        if (queryNum != prevSeqId) {


            if(alignment != NULL)
            {
                if(alignments_findGoodAlignments_ncbi(
                            alignment,
                            PSSMatrix_arr[queryNum], scoreMatrix,
                            queryNum, ungappedExtension_new, dp_mem,
                            tree, private_tree, BlastHSP))
                {
                    int kk;
                    for(kk = 0; kk < alignment->numExtensions; kk++)
                    {
                        goodExtensionBuf[numGoodExtensions + new_numUngappedExtSeq + kk] =
                            alignment->ungappedExtensions[kk];
                    }

                    alignment->ungappedExtensions = goodExtensionBuf + numGoodExtensions + new_numUngappedExtSeq;

                    numGoodAlign++;
                    new_numUngappedExtSeq += alignment->numExtensions;
                }
            }

            uint4 descriptionStart =
                sequenceData[sequenceCount].descriptionStart;
            uint4 encodedLength = sequenceData[sequenceCount].encodedLength;
            uint4 descriptionLength =
                sequenceData[sequenceCount].descriptionLength;
            unsigned char *subject = sequenceData[sequenceCount].sequence;
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
            alignment->ungappedExtensions = goodExtensionBuf + numGoodExtensions + ii;
            alignment->gappedExtensions = NULL;
            alignment->unpackRegions = NULL;
            alignment->edits = NULL;
            alignment->numExtensions = 0;
            alignment->queryCount = queryNum;
            numAlignBlk = 0;
        }

        numAlignBlk++;
        alignment->numExtensions++;

        prevSeqId = queryNum;
    }

    if(numAlignBlk > 0)
    {

        //fprintf(stderr, "seq: %d numExt: %d\n", sequenceCount, alignment->numExtensions);
        if(alignments_findGoodAlignments_ncbi(
                    alignment,
                    PSSMatrix_arr[queryNum], scoreMatrix,
                    queryNum, ungappedExtension_new, dp_mem,
                    tree, private_tree, BlastHSP))
        {
            int kk;
            for(kk = 0; kk < alignment->numExtensions; kk++)
            {
                goodExtensionBuf[numGoodExtensions + new_numUngappedExtSeq + kk] =
                    alignment->ungappedExtensions[kk];
            }

            alignment->ungappedExtensions = goodExtensionBuf + numGoodExtensions + new_numUngappedExtSeq;
            numGoodAlign++;
            new_numUngappedExtSeq += alignment->numExtensions;
        }
    }

    numGoodExtensions += new_numUngappedExtSeq;

#endif


#ifdef PROFILE
    RDTSC_END;
    GET_RDTSC_ATOMIC(&profile->blast_ungappedExtCycle);
#endif

    (*goodAlignCount) = numGoodAlign;
    (*goodExtensionCount) = numGoodExtensions;

}

