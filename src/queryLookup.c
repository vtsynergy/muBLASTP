#include "blast.h"

struct initialWord_protein_query *proteinLookup_query;

int numQueryBlk = 0;
struct initialWord_protein_query *proteinLookup_query_blk[NUM_QUERY_BLOCK] = {NULL};

PV_ARRAY_TYPE *pv;
PV_ARRAY_TYPE *pv_arr[NUM_QUERY_BLOCK] = {NULL};

uint32_t *external_positions_blk[NUM_QUERY_BLOCK] = {NULL};
int total_positions_blk[NUM_QUERY_BLOCK] = {0};

uint32_t *external_positions = NULL;
int total_positions = 0;

void proteinLookup_query_initial_blk(int4 numCodes, int wordLength, int block_id) {
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
    proteinLookup_query_blk[block_id] = initialLookup =
        (struct initialWord_protein_query *)global_malloc(
                sizeof(struct initialWord_protein_query) * numEntries);

    proteinLookup_numWords = numEntries;

    pv_arr[block_id] = (PV_ARRAY_TYPE *)calloc(
            (proteinLookup_numWords >> PV_ARRAY_BTS) + 1, sizeof(PV_ARRAY_TYPE));

    // Iterate through every possible codeword
    codeword = 0;

    codeword = 0;
    while (codeword < numEntries) {
        // Initialize list of query positions as empty
        initialLookup[codeword].numQueryPositions = 0;
        initialLookup[codeword].allocQueryPositions = 0;
        initialLookup[codeword].querySequencePositions = NULL;
        codeword++;
    }
}

void proteinLookup_query_initial(int4 numCodes, int wordLength) {
    struct initialWord_protein_query *initialLookup, *initialWord;
    struct initialWord_neighborLookup *initialLookup_n;
    uint4 codeword, numEntries;

    wordLookupDFA_numCodes = numCodes;
    wordLookupDFA_wordLength = wordLength;
    numEntries = ceil(pow(numCodes, wordLength));
    int4 proteinLookup_wordLength = wordLength;

    // Declare memory for initial lookup table
    proteinLookup_query = initialLookup =
        (struct initialWord_protein_query *)global_malloc(
                sizeof(struct initialWord_protein_query) * numEntries);

    proteinLookup_numWords = numEntries;

    pv = (PV_ARRAY_TYPE *)calloc( (proteinLookup_numWords >> PV_ARRAY_BTS) + 1, sizeof(PV_ARRAY_TYPE));

    // Iterate through every possible codeword
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
        //free(proteinLookup_query[codeword].querySequencePositions);
        codeword++;
    }

    free(proteinLookup_query);
    free(external_positions);
    total_positions = 0;
    free(pv);
}

void queryIndex_free_blk() {
    int t;
    for(t = 0; t < numQueryBlk; t++)
    {
        if(proteinLookup_query_blk[t] == NULL)
            break;

        int codeword = 0;
        while (codeword < proteinLookup_numWords) {
            // Initialize list of query positions as empty
            proteinLookup_query_blk[t][codeword].numQueryPositions = 0;
            proteinLookup_query_blk[t][codeword].allocQueryPositions = 0;
            //free(proteinLookup_query_blk[t][codeword].querySequencePositions);
            codeword++;
        }

        free(external_positions_blk[t]);
        total_positions_blk[t] = 0;
        free(proteinLookup_query_blk[t]);
        free(pv_arr[t]);
    }

    //free(proteinLookup_query_blk);
}

void proteinLookup_query_add_blk(uint4 sequenceNum, unsigned char *sequence,
        int4 subjectLength, int4 wordLength,
        struct scoreMatrix scoreMatrix, int block_id) {

    uint4 codeword, numEntries, byteNumber, queryPosition;
    struct initialWord_protein_query *initialLookup, *initialWord;

    // Number of entries in the table
    numEntries = proteinLookup_numWords;

    initialLookup = proteinLookup_query_blk[block_id];

    // Slide a word-sized window across the query
    queryPosition = 0;

    int4 numNeighbours;
    struct neighbour *neighbours = (struct neighbour *)global_malloc(
            sizeof(struct neighbour) * proteinLookup_numWords);

    while (queryPosition < subjectLength - wordLength + 1) {

        numNeighbours = 0;
        wordLookupSM_getNeighbours(sequence, scoreMatrix, queryPosition,
                &numNeighbours, neighbours);

        total_positions_blk[block_id] += numNeighbours;

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

        total_positions += numNeighbours;

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


int compareQuerySeqId(const void *a, const void *b) {
    int pos_a = *((uint32_t *)a) & 0xffff;
    int pos_b = *((uint32_t *)b) & 0xffff;
    return pos_b - pos_a;
}

void proteinLookup_sort()
{
    int ii, numEmptyCell = 0, numSmallCell = 0;
    int external_positions_cnt = 0;
    external_positions = (uint32_t *)malloc(sizeof(uint32_t) * total_positions);
    for(ii = 0; ii < proteinLookup_numWords; ii++)
    {
        if(proteinLookup_query[ii].numQueryPositions > 0)
        {
            PV_SET(pv, ii, PV_ARRAY_BTS);
            //qsort(proteinLookup_query[ii].querySequencePositions, 
            //proteinLookup_query[ii].numQueryPositions, sizeof(uint32_t), compareQueryPos);

            if(proteinLookup_query[ii].numQueryPositions <= 3)
            {
                memcpy(proteinLookup_query[ii].embeddedPositions, 
                        proteinLookup_query[ii].querySequencePositions, 
                        sizeof(uint32_t) * proteinLookup_query[ii].numQueryPositions);
                numSmallCell++;
                //free(proteinLookup_query[ii].querySequencePositions);
            }

            //else
            {
                memcpy(external_positions + external_positions_cnt, 
                        proteinLookup_query[ii].querySequencePositions, 
                        sizeof(uint32_t) * proteinLookup_query[ii].numQueryPositions);
                free(proteinLookup_query[ii].querySequencePositions);
                proteinLookup_query[ii].querySequencePositions = external_positions + external_positions_cnt;
                external_positions_cnt+=proteinLookup_query[ii].numQueryPositions;
            }
        }
        else
        {
            numEmptyCell++;
        }
    }

    //fprintf(stderr, "numEmptyCell: %d numSmallCell: %d totalNumCell: %d\n", numEmptyCell, numSmallCell, proteinLookup_numWords);
}

void proteinLookup_sort_blk(int block_id)
{
    int ii, numEmptyCell = 0;
    int numSmallCell = 0;
    external_positions_blk[block_id] = (uint32_t *)malloc(sizeof(uint32_t) * total_positions_blk[block_id]);
    int external_positions_cnt = 0;
    for(ii = 0; ii < proteinLookup_numWords; ii++)
    {
        if(proteinLookup_query_blk[block_id][ii].numQueryPositions > 0)
        {
            PV_SET(pv_arr[block_id], ii, PV_ARRAY_BTS);
            //qsort(proteinLookup_query_blk[block_id][ii].querySequencePositions, 
            //proteinLookup_query_blk[block_id][ii].numQueryPositions, sizeof(uint32_t), compareQueryPos);

            if(proteinLookup_query_blk[block_id][ii].numQueryPositions <= 3)
            {
                memcpy(proteinLookup_query_blk[block_id][ii].embeddedPositions, 
                        proteinLookup_query_blk[block_id][ii].querySequencePositions, 
                        sizeof(uint32_t) * proteinLookup_query_blk[block_id][ii].numQueryPositions);
                numSmallCell++;
                //free(proteinLookup_query_blk[block_id][ii].querySequencePositions);
            }
            //else
            {
                memcpy(external_positions_blk[block_id] + external_positions_cnt, 
                        proteinLookup_query_blk[block_id][ii].querySequencePositions, 
                        sizeof(uint32_t) * proteinLookup_query_blk[block_id][ii].numQueryPositions);
                free(proteinLookup_query_blk[block_id][ii].querySequencePositions);
                proteinLookup_query_blk[block_id][ii].querySequencePositions = external_positions_blk[block_id] + external_positions_cnt;
                external_positions_cnt+=proteinLookup_query_blk[block_id][ii].numQueryPositions;
            }

        }
        else
        {
            numEmptyCell++;
        }
    }
}
