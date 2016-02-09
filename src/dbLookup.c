#include "blast.h"

struct initialWord_protein_db *proteinLookup_db;
struct proteinLookup_db_blk *proteinLookup_db_b;
struct initialWord_neighborLookup *neighborLookup;

uint4 proteinLookup_numBlocks = 0;
uint4 proteinLookup_numWords = 0;
uint4 maxNumSeqBlk = 0;

uint8 total_allocHits = 0;
uint8 numHits_blk = 0;
uint8 numHits_db = 0;

#define min(a, b) ((a) <= (b) ? (a) : (b))

void write_dbIdxBlock(char *write_dbLookupFilename, int blockNum);
void write_dbLookupAux(char *write_dbLookupFilename);
//uint32_t **subSequencePositions_t;
//uint4 *allocSubPositions;

inline int4 getCodeword(unsigned char *codes, int4 wordLength) {
    int4 codeword;
    uint4 codeCount;

    codeword = 0;
    codeCount = 0;
    while (codeCount < wordLength) {
        codeword *= wordLookupDFA_numCodes;
        if (codes[codeCount] < wordLookupDFA_numCodes)
            codeword += codes[codeCount];
        codeCount++;
    }

    return codeword;
}



void proteinLookup_db_initial(int4 numCodes, int wordLength) {
    struct initialWord_protein_db *initialLookup, *initialWord;
    struct initialWord_neighborLookup *initialLookup_n;
    uint4 codeword, numEntries;

    wordLookupDFA_numCodes = numCodes;
    wordLookupDFA_wordLength = wordLength;
    numEntries = ceil(pow(numCodes, wordLength));
    int4 proteinLookup_wordLength = wordLength;

    uint4 numBlocks = ceil((float)readdb_numberOfLetters / dbIdx_block_size) + 1;
    proteinLookup_numBlocks = numBlocks;

    // Declare memory for initial DB index blocks
    proteinLookup_db_b = (struct proteinLookup_db_blk *)global_malloc(
            sizeof(struct proteinLookup_db_blk) * numBlocks);

    // Declare memory for initial lookup table
    proteinLookup_db = initialLookup =
        (struct initialWord_protein_db *)global_malloc(
                sizeof(struct initialWord_protein_db) * numEntries);

    proteinLookup_numWords = numEntries;

    // Iterate through every possible codeword
    codeword = 0;
    while (codeword < numBlocks) {
        // Initialize list of query positions as empty
        //proteinLookup_db_b[codeword].proteinLookup_db =
        //initialLookup + codeword * numEntries;
        proteinLookup_db_b[codeword].subPositionOffset = (uint4 *)global_malloc(sizeof(uint4) * (numEntries + 1));
        proteinLookup_db_b[codeword].dbIdxblk_longestSeq = 0;
        proteinLookup_db_b[codeword].numSeqBlk = 0;
        codeword++;
    }

    codeword = 0;
    while (codeword < numEntries) {
        // Initialize list of query positions as empty
        initialLookup[codeword].numSubPositions = 0;
        initialLookup[codeword].allocSubPositions = 0;
        initialLookup[codeword].subSequencePositions = NULL;
        //initialLookup[codeword].subPositionOffset = 0;
        codeword++;
    }
}



void free_dbindex() {
    uint4 codeword;
    codeword = 0;
    while (codeword < proteinLookup_numBlocks) {
        free(proteinLookup_db_b[codeword].subSequencePositions);
        free(proteinLookup_db_b[codeword].subPositionOffset);
        codeword++;
    }
    free(proteinLookup_db_b);
    free(proteinLookup_db);
}


void proteinLookup_db_sub(uint4 sequenceNum, unsigned char *sequence,
        int4 subjectLength, int4 wordLength, int blockNum,
        int startSeqNum) {

    uint4 codeword, numEntries, byteNumber, queryPosition;
    struct initialWord_protein_db *initialLookup, *initialWord;

    // Number of entries in the table
    numEntries = proteinLookup_numWords;

    initialLookup = proteinLookup_db;

    // Slide a word-sized window across the query
    queryPosition = 0;

    while (queryPosition < subjectLength - wordLength + 1) {

        codeword = getCodeword(sequence + queryPosition, wordLength);

        initialWord = initialLookup + codeword;

        initialWord->numSubPositions++;
        if (initialWord->numSubPositions > initialWord->allocSubPositions) {
            // total_allocHits += allocSubPositions[codeword];
            initialWord->allocSubPositions = initialWord->allocSubPositions == 0
                ? 50
                : 2 * initialWord->allocSubPositions;
            initialWord->subSequencePositions = (subPos_t *)global_realloc(
                    initialWord->subSequencePositions,
                    sizeof(subPos_t) * initialWord->allocSubPositions);
        }

        // initialWord->subSequencePositions[initialWord->numSubPositions - 1] =
        //((sequenceNum - blockNum * dbIdx_block_size) << 16) + queryPosition;

        initialWord->subSequencePositions[initialWord->numSubPositions - 1].subOff = queryPosition;
        initialWord->subSequencePositions[initialWord->numSubPositions - 1].seqId = sequenceNum - startSeqNum;

        //initialWord->subSequencePositions[initialWord->numSubPositions - 1] =
        //(queryPosition << 16) + (sequenceNum - startSeqNum);

        numHits_blk++;
        queryPosition++;
    }
}

void neighbourLookup_init() {
    int4 numEntries = proteinLookup_numWords;

    neighborLookup = (struct initialWord_neighborLookup *)global_malloc(
            sizeof(struct initialWord_neighborLookup) * numEntries);

    int4 codeword = 0;
    while (codeword < numEntries) {
        neighborLookup[codeword].numNeighbours = 0;
        neighborLookup[codeword].neighbours = NULL;
        codeword++;
    }
}


void neighbourLookup_free() {
    int4 numEntries = proteinLookup_numWords;
    int4 codeword = 0;
    while (codeword < numEntries) {
        free(neighborLookup[codeword].neighbours);
        codeword++;
    }
    free(neighborLookup);
}

void neighbourLookup_build(struct PSSMatrix PSSMatrix,
        struct scoreMatrix scoreMatrix, int4 wordLength) {
    int4 queryPosition = 0;
    int4 numNeighbours;
    int4 codeword;
    int4 numWords = proteinLookup_numWords;
    struct neighbour *neighbours =
        (struct neighbour *)global_malloc(sizeof(struct neighbour) * numWords);

    while (queryPosition < PSSMatrix.length - wordLength + 1) {
        codeword =
            getCodeword(PSSMatrix.bestMatchCodes + queryPosition, wordLength);

        if (neighborLookup[codeword].numNeighbours == 0) {
            numNeighbours = 0;
            // wordLookupDFA_getNeighbours(PSSMatrix, queryPosition, &numNeighbours,
            // neighbours);

            wordLookupSM_getNeighbours(PSSMatrix.bestMatchCodes, scoreMatrix,
                    queryPosition, &numNeighbours, neighbours);

            neighborLookup[codeword].numNeighbours = numNeighbours;
            neighborLookup[codeword].neighbours =
                (int4 *)global_malloc(sizeof(int4) * numNeighbours);

            while (numNeighbours > 0) {
                numNeighbours--;
                neighborLookup[codeword].neighbours[numNeighbours] =
                    neighbours[numNeighbours].codeword;
            }
        }

        // printf("%d %d\n", codeword, neighborLookup[codeword].numNeighbours);
        queryPosition++;
    }
    free(neighbours);
}


void proteinLookup_db_build(int4 numCodes, int wordLength,
        struct scoreMatrix scoreMatrix,
        char *write_dbLookupFilename) {
    proteinLookup_db_initial(numCodes, wordLength);
    int ii, jj;

    printf("numBlocks: %d numWords: %d\n", proteinLookup_numBlocks,
            proteinLookup_numWords);

    int seqStartBlk = 0;
    for (jj = 0; jj < proteinLookup_numBlocks; jj++) {
        int numLetterBlk = 0;
        numHits_blk = 0;
        printf("buiding index block %d of %d ... ", jj + 1,
                proteinLookup_numBlocks);
        proteinLookup_db_b[jj].seqOffset = seqStartBlk;
        for (ii = seqStartBlk;
                ii < readdb_numberOfSequences && numLetterBlk < dbIdx_block_size;
                ii++) {

            numLetterBlk += readdb_sequenceData[ii].sequenceLength;
            proteinLookup_db_sub(ii, readdb_sequenceData[ii].sequence,
                    readdb_sequenceData[ii].sequenceLength, wordLength,
                    jj, seqStartBlk);
            proteinLookup_db_b[jj].dbIdxblk_longestSeq =
                MAX(proteinLookup_db_b[jj].dbIdxblk_longestSeq,
                        readdb_sequenceData[ii].sequenceLength);
            if((ii - seqStartBlk) >= ((1 << 15) - 1))
            {
                break;
            }
        }

        if (numLetterBlk >= dbIdx_block_size) {
            numLetterBlk -= readdb_sequenceData[ii - 1].sequenceLength;
        }

        proteinLookup_db_b[jj].numSeqBlk = ii - seqStartBlk;

        //if(proteinLookup_db_b[jj].numSeqBlk >= (1 << 15))
        //{
            //fprintf(stderr, "Err: too many sequences (%d) in a block, please use a smaller block size.\n", proteinLookup_db_b[jj].numSeqBlk);
            //exit(1);

            //}
        seqStartBlk = ii;
        numHits_db += numHits_blk;

#if 0 
        uint4 numSubPositions = 0;
        for (ii = 0; ii < proteinLookup_numWords; ii++) {
            // proteinLookup_db_b[jj].proteinLookup_db[ii].wordCode = ii;
            std_sort(subSequencePositions_t[ii],
                    proteinLookup_db[jj * proteinLookup_numWords + ii].numSubPositions);
            proteinLookup_db_b[jj].subPositionOffset[ii] = numSubPositions;
            numSubPositions += proteinLookup_db[jj * proteinLookup_numWords + ii].numSubPositions;
        }
        proteinLookup_db_b[jj].subPositionOffset[proteinLookup_numWords] = numSubPositions;
#else
        uint4 numSubPositions = 0;
        for (ii = 0; ii < proteinLookup_numWords; ii++) {
            proteinLookup_db_b[jj].subPositionOffset[ii] = numSubPositions;
            numSubPositions += proteinLookup_db[ii].numSubPositions;
            qsort(proteinLookup_db[ii].subSequencePositions,
                    proteinLookup_db[ii].numSubPositions, sizeof(subPos_t), comparePos_subjectOffset);
        }
        proteinLookup_db_b[jj].subPositionOffset[proteinLookup_numWords] = numSubPositions;
#endif

        maxNumSeqBlk = MAX(proteinLookup_db_b[jj].numSeqBlk, maxNumSeqBlk);

        printf("numSeq: %u numLetterBlk: %u total_size: %lu = %lu + %lu "
                "(KB)\n",
                proteinLookup_db_b[jj].numSeqBlk, numLetterBlk,
                (sizeof(struct initialWord_protein_db) * proteinLookup_numWords +
                 sizeof(subPos_t) * numHits_blk) >>
                10,
                sizeof(struct initialWord_protein_db) * proteinLookup_numWords >> 10,
                sizeof(subPos_t) * numHits_blk >> 10);

        write_dbIdxBlock(write_dbLookupFilename, jj);

        for (ii = 0; ii < proteinLookup_numWords; ii++) {
            free(proteinLookup_db[ii].subSequencePositions);
            proteinLookup_db[ii].subSequencePositions = NULL;
            proteinLookup_db[ii].allocSubPositions = 0;
            proteinLookup_db[ii].numSubPositions = 0;
        }

    }

    // detect empty blocks and skip them
    int numBlocks = proteinLookup_numBlocks;
    while (proteinLookup_db_b[numBlocks - 1].numSeqBlk == 0) {
        numBlocks--;
    }

    proteinLookup_numBlocks = numBlocks;

    write_dbLookupAux(write_dbLookupFilename);

    printf("totalHit: %ld numBlocks: %d\n", numHits_db, proteinLookup_numBlocks);

    //free(proteinLookup_db_b);
    //free(allocSubPositions);
    // proteinLookup_db_print();
}

void write_dbIdxBlock(char *write_dbLookupFilename, int blockNum) {
    char *write_dbIdxBlockFilename;
    FILE *write_dbIdxBlockFile;

    write_dbIdxBlockFilename =
        (char *)global_malloc(strlen(write_dbLookupFilename) + 20);
    sprintf(write_dbIdxBlockFilename, "%s.dbLookup.%d", write_dbLookupFilename,
            blockNum);
    if ((write_dbIdxBlockFile = fopen(write_dbIdxBlockFilename, "w")) == NULL) {
        fprintf(stderr, "Error opening file %s for writing\n",
                write_dbIdxBlockFilename);
        exit(-1);
    }

    int ii;
    if (fwrite(proteinLookup_db_b[blockNum].subPositionOffset,
                sizeof(uint4), proteinLookup_numWords + 1,
                write_dbIdxBlockFile) != proteinLookup_numWords + 1) {
        fprintf(stderr, "Error writing header to dbIdxBlock file %s\n",
                write_dbIdxBlockFilename);
        exit(-1);
    }

    for (ii = 0; ii < proteinLookup_numWords; ii++) {
        int4 numSubPositions =
            proteinLookup_db[ii].numSubPositions;
        //uint32_t *subSequencePositions = proteinLookup_db[ii].subSequencePositions;
        if (numSubPositions == 0)
            continue;
        if (fwrite(proteinLookup_db[ii].subSequencePositions, sizeof(subPos_t), numSubPositions,
                    write_dbIdxBlockFile) != numSubPositions) {
            fprintf(stderr, "Error writing numSubPositions to dbIdxBlock file %s\n",
                    write_dbIdxBlockFilename);
            exit(-1);
        }
    }

    fclose(write_dbIdxBlockFile);
    free(write_dbIdxBlockFilename);
}

void read_dbIdxBlock(char *read_dbLookupFilename, int blockNum) {
    char *read_dbIdxBlockFilename;
    FILE *read_dbIdxBlockFile;

    read_dbIdxBlockFilename =
        (char *)global_malloc(strlen(read_dbLookupFilename) + 20);
    sprintf(read_dbIdxBlockFilename, "%s.dbLookup.%d", read_dbLookupFilename,
            blockNum);
    if ((read_dbIdxBlockFile = fopen(read_dbIdxBlockFilename, "r")) == NULL) {
        fprintf(stderr, "Error opening file %s for writing\n",
                read_dbIdxBlockFilename);
        exit(-1);
    }

    proteinLookup_db_b[blockNum].subPositionOffset = (uint4 *)global_malloc(sizeof(uint4) * (proteinLookup_numWords + 1));
    if (fread(proteinLookup_db_b[blockNum].subPositionOffset, sizeof(uint4),
                proteinLookup_numWords + 1,
                read_dbIdxBlockFile) != proteinLookup_numWords + 1) {
        fprintf(stderr, "Error reading LoopkupHeader to dbIdxBlock file %s\n",
                read_dbIdxBlockFilename);
        exit(-1);
    }

    int totalHits = proteinLookup_db_b[blockNum].subPositionOffset[proteinLookup_numWords];

    //uint32_t *subSequencePositions =
    //(uint32_t *)global_malloc(sizeof(uint32_t) * totalHits);
    proteinLookup_db_b[blockNum].subSequencePositions = (subPos_t *)global_malloc(sizeof(subPos_t) * totalHits);
    if (fread(proteinLookup_db_b[blockNum].subSequencePositions, sizeof(subPos_t), totalHits,
                read_dbIdxBlockFile) != totalHits) {
        fprintf(stderr, "Error reading numSubPositions to dbIdxBlock file %s\n",
                read_dbIdxBlockFilename);
        exit(-1);
    }

    fclose(read_dbIdxBlockFile);
    free(read_dbIdxBlockFilename);
}

struct dbLookupAux {
    uint4 proteinLookup_numBlocks;
    uint4 proteinLookup_numWords;
    uint4 proteinLookup_wordLength;
    int4 proteinLookup_numCodes;
    int4 dbIdx_block_size;
    uint8 totalNumPositions;
    int4 maxNumSeqBlk;
};

void write_dbLookupAux(char *write_dbLookupFilename) {
    char *write_dbLookupAuxFilename;
    FILE *write_dbLookupAuxFile;

    write_dbLookupAuxFilename =
        (char *)global_malloc(strlen(write_dbLookupFilename) + 10);
    sprintf(write_dbLookupAuxFilename, "%s.dbLookup", write_dbLookupFilename);
    // Open dbLookup file for writing
    if ((write_dbLookupAuxFile = fopen(write_dbLookupAuxFilename, "w")) == NULL) {
        fprintf(stderr, "Error opening file %s for writing\n",
                write_dbLookupAuxFilename);
        exit(-1);
    }

    struct dbLookupAux dbLookupAux;
    dbLookupAux.proteinLookup_numBlocks = proteinLookup_numBlocks;
    dbLookupAux.proteinLookup_numWords = proteinLookup_numWords;
    dbLookupAux.proteinLookup_wordLength = wordLookupDFA_wordLength;
    dbLookupAux.proteinLookup_numCodes = wordLookupDFA_numCodes;
    dbLookupAux.dbIdx_block_size = dbIdx_block_size;
    dbLookupAux.totalNumPositions = numHits_db;
    dbLookupAux.maxNumSeqBlk = maxNumSeqBlk;

    if (fwrite(&dbLookupAux, sizeof(struct dbLookupAux), 1,
                write_dbLookupAuxFile) != 1) {
        fprintf(stderr, "Error writing data to dbLookup aux file %s\n",
                write_dbLookupAuxFilename);
        exit(-1);
    }

    if (fwrite(proteinLookup_db_b, sizeof(struct proteinLookup_db_blk),
                proteinLookup_numBlocks,
                write_dbLookupAuxFile) != proteinLookup_numBlocks) {
        fprintf(stderr, "Error writing data to dbLookup aux file %s\n",
                write_dbLookupAuxFilename);
        exit(-1);
    }

    // if (fwrite(numSeqBlk, sizeof(int), proteinLookup_numBlocks,
    // write_dbLookupAuxFile) != proteinLookup_numBlocks)
    //{
    // fprintf(stderr, "Error writing data to dbLookup aux file %s\n",
    // write_dbLookupAuxFilename);
    // exit(-1);
    //}

    free(write_dbLookupAuxFilename);
    fclose(write_dbLookupAuxFile);
}

void read_dbLookupAux(char *read_dbLookupFilename) {
    char *read_dbLookupAuxFilename;
    FILE *read_dbLookupAuxFile;

    read_dbLookupAuxFilename =
        (char *)global_malloc(strlen(read_dbLookupFilename) + 10);
    sprintf(read_dbLookupAuxFilename, "%s.dbLookup", read_dbLookupFilename);

    // Open dbLookup file for writing
    if ((read_dbLookupAuxFile = fopen(read_dbLookupAuxFilename, "r")) == NULL) {
        fprintf(stderr, "Error opening file %s for reading\n",
                read_dbLookupAuxFilename);
        exit(-1);
    }

    struct dbLookupAux dbLookupAux;

    fread(&dbLookupAux, sizeof(struct dbLookupAux), 1, read_dbLookupAuxFile);

    blast_numBlocks = proteinLookup_numBlocks =
        dbLookupAux.proteinLookup_numBlocks;
    proteinLookup_numWords = dbLookupAux.proteinLookup_numWords;
    wordLookupDFA_wordLength = dbLookupAux.proteinLookup_wordLength;
    wordLookupDFA_numCodes = dbLookupAux.proteinLookup_numCodes;
    dbIdx_block_size = dbLookupAux.dbIdx_block_size;
    numHits_db = dbLookupAux.totalNumPositions;
    maxNumSeqBlk = dbLookupAux.maxNumSeqBlk;

    proteinLookup_db_b = (struct proteinLookup_db_blk *)malloc(
            sizeof(struct proteinLookup_db_blk) * proteinLookup_numBlocks);
    fread(proteinLookup_db_b, sizeof(struct proteinLookup_db_blk),
            proteinLookup_numBlocks, read_dbLookupAuxFile);

#if 0
    fprintf(stderr,
            "dbIndex info:\nproteinLookup_numBlocks: %d\nproteinLookup_numWords: %d\n"
            "proteinLookup_wordLength: %d\nproteinLookup_numCodes: "
            "%d\nproteinLookup_dbIdxBlockSize: %d (K letters)\ntotal_numPositions: "
            "%ld\nmaxNumSeqBlk: %d\nlongestSequenceLength: %d\n\n",
            proteinLookup_numBlocks, proteinLookup_numWords, wordLookupDFA_wordLength,
            wordLookupDFA_numCodes, dbIdx_block_size / 1024, numHits_db,
            maxNumSeqBlk, readdb_longestSequenceLength);
#endif

    proteinLookup_db = (struct initialWord_protein_db *)malloc(
            sizeof(struct initialWord_protein_db) * proteinLookup_numWords *
            proteinLookup_numBlocks);

    free(read_dbLookupAuxFilename);
    fclose(read_dbLookupAuxFile);
}

void read_dbLookup(char *read_dbLookupFilename) {

#ifdef PROFILE
    struct timeval start, end;
    gettimeofday(&start, NULL);
#endif
    read_dbLookupAux(read_dbLookupFilename);

    int ii;
    for (ii = 0; ii < proteinLookup_numBlocks; ii++) {
        // printf("reading DB_BLOCK %d ", ii);
        read_dbIdxBlock(read_dbLookupFilename, ii);
    }

#ifdef PROFILE
    gettimeofday(&end, NULL);
    long read_time = ((end.tv_sec * 1000000 + end.tv_usec) -
            (start.tv_sec * 1000000 + start.tv_usec));
    //fprintf(stderr, "Index read time: %f\n", (float)read_time * 1e-6);
#endif

}

void write_dbLookup(char *write_dbLookupFilename) {
    write_dbLookupAux(write_dbLookupFilename);

    int ii;
    for (ii = 0; ii < proteinLookup_numBlocks; ii++) {
        // printf("writing DB_BLOCK %d\n", ii);
        write_dbIdxBlock(write_dbLookupFilename, ii);
    }
}
