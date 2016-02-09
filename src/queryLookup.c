#include "blast.h"

struct initialWord_protein_query *proteinLookup_query;

PV_ARRAY_TYPE *pv;


void proteinLookup_query_initial(int4 numCodes, int wordLength) {
    struct initialWord_protein_query *initialLookup, *initialWord;
    struct initialWord_neighborLookup *initialLookup_n;
    uint4 codeword, numEntries;

    wordLookupDFA_numCodes = numCodes;
    wordLookupDFA_wordLength = wordLength;
    numEntries = ceil(pow(numCodes, wordLength));
    int4 proteinLookup_wordLength = wordLength;

    //uint4 numBlocks = ceil((float)readdb_numberOfLetters / dbIdx_block_size) + 1;
    //proteinLookup_numBlocks = numBlocks;

    // Declare memory for initial DB index blocks
    //proteinLookup_db_b = (struct proteinLookup_db_blk *)global_malloc(
    //sizeof(struct proteinLookup_db_blk) * numBlocks);

    // Declare memory for initial lookup table
    proteinLookup_query = initialLookup =
        (struct initialWord_protein_query *)global_malloc(
                sizeof(struct initialWord_protein_query) * numEntries);

    proteinLookup_numWords = numEntries;

    pv = (PV_ARRAY_TYPE *)calloc( (proteinLookup_numWords >> PV_ARRAY_BTS) + 1, sizeof(PV_ARRAY_TYPE));

    // Iterate through every possible codeword
    codeword = 0;
    //while (codeword < numBlocks) {
        //// Initialize list of query positions as empty
        ////proteinLookup_db_b[codeword].proteinLookup_db =
        ////initialLookup + codeword * numEntries;
        //proteinLookup_db_b[codeword].subPositionOffset = (uint4 *)global_malloc(sizeof(uint4) * (numEntries + 1));
        //proteinLookup_db_b[codeword].dbIdxblk_longestSeq = 0;
        //proteinLookup_db_b[codeword].numSeqBlk = 0;
        //codeword++;
    //}

    codeword = 0;
    while (codeword < numEntries) {
        // Initialize list of query positions as empty
        initialLookup[codeword].numQueryPositions = 0;
        initialLookup[codeword].allocQueryPositions = 0;
        initialLookup[codeword].querySequencePositions = NULL;
        codeword++;
    }
}

void queryIndex_free() {
    int codeword = 0;
    while (codeword < proteinLookup_numWords) {
        // Initialize list of query positions as empty
        proteinLookup_query[codeword].numQueryPositions = 0;
        proteinLookup_query[codeword].allocQueryPositions = 0;
        free(proteinLookup_query[codeword].querySequencePositions);
        codeword++;
    }

    free(proteinLookup_query);
}

void proteinLookup_query_add(uint4 sequenceNum, unsigned char *sequence,
        int4 subjectLength, int4 wordLength,
        struct scoreMatrix scoreMatrix) {

    uint4 codeword, numEntries, byteNumber, queryPosition;
    struct initialWord_protein_query *initialLookup, *initialWord;

    // Number of entries in the table
    numEntries = proteinLookup_numWords;

    initialLookup = proteinLookup_query;

    // Slide a word-sized window across the query
    queryPosition = 0;

    int4 numNeighbours;
    struct neighbour *neighbours = (struct neighbour *)global_malloc(
            sizeof(struct neighbour) * proteinLookup_numWords);

    while (queryPosition < subjectLength - wordLength + 1) {

        numNeighbours = 0;
        wordLookupSM_getNeighbours(sequence, scoreMatrix, queryPosition,
                &numNeighbours, neighbours);

        while (numNeighbours > 0) {
            numNeighbours--;
            int codeword = neighbours[numNeighbours].codeword;
            initialWord = initialLookup + codeword;
            initialWord->numQueryPositions++;

            if (initialWord->numQueryPositions > initialWord->allocQueryPositions) {
                initialWord->allocQueryPositions += 20; 
                    initialWord->querySequencePositions =
                    (uint32_t *)global_realloc(
                            initialWord->querySequencePositions,
                            sizeof(uint32_t) * initialWord->allocQueryPositions);
            }

            // initialWord->subSequencePositions[initialWord->numSubPositions - 1] =
            //((sequenceNum - blockNum * dbIdx_block_size) << 16) + queryPosition;

            initialWord->querySequencePositions[initialWord->numQueryPositions - 1] =
                (queryPosition << 16) + (sequenceNum);
        }
        queryPosition++;
    }
    free(neighbours);
}

int compareQueryPos(const void *a, const void *b) {
    int pos_a = *((uint32_t *)a);
    int pos_b = *((uint32_t *)b);
    return pos_a - pos_b;
}

void proteinLookup_sort()
{
    int ii;
    for(ii = 0; ii < proteinLookup_numWords; ii++)
    {
        if(proteinLookup_query[ii].numQueryPositions > 0)
        {
            PV_SET(pv, ii, PV_ARRAY_BTS);
            qsort(proteinLookup_query[ii].querySequencePositions, proteinLookup_query[ii].numQueryPositions, sizeof(uint32_t), compareQueryPos);
        }
    }
}
