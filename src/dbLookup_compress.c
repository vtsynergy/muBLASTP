#include "blast.h"

struct initialWord_protein_db_cp *proteinLookup_db_cp;
struct proteinLookup_db_blk_cp  *proteinLookup_db_blk_cp;

uint4 numPosBlkCP;
uint8 numPosDBCP;

#define CHAR_MAX 255

typedef struct{
    uint4 subOff;
    uint4 seqId;
}initialSubPos;

struct initialWord_protein_db_t
{
    uint4  numSubPositions;
    uint4  allocSubPositions;
    initialSubPos *subSequencePositions;
};

int compareSubPos(const void *a, const void *b) {
    initialSubPos *pos_a = (initialSubPos *)a;
    initialSubPos *pos_b = (initialSubPos *)b;

    return pos_a->subOff - pos_b->subOff; 
}

struct initialWord_protein_db_t *proteinLookup_db_temp;
//struct initialWord_protein_db_t **proteinLookup_db_blk_temp;

void write_dbLookupAux_cp(char *write_dbLookupFilename);
void write_dbIdxBlock_cp(char *write_dbLookupFilename, int blockNum);

//uint4 proteinLookup_numBlocks = 0;

void proteinLookup_db_cp_initial(int4 numCodes, int wordLength) {
    struct initialWord_protein_db_cp *initialLookup, *initialWord;
    struct initialWord_protein_db_t *initialLookup_t;
    uint4 codeword, numEntries;

    wordLookupDFA_numCodes = numCodes;
    wordLookupDFA_wordLength = wordLength;
    numEntries = ceil(pow(numCodes, wordLength));
    int4 proteinLookup_wordLength = wordLength;

    uint4 numBlocks = ceil((float)readdb_numberOfLetters / dbIdx_block_size) + 1;
    proteinLookup_numBlocks = numBlocks;

    // Declare memory for initial DB index blocks
    proteinLookup_db_blk_cp = (struct proteinLookup_db_blk_cp *)global_malloc(
            sizeof(struct proteinLookup_db_blk_cp) * numBlocks);

    //proteinLookup_db_blk_temp = (struct initialWord_protein_db_t **)global_malloc(
    //sizeof(struct initialWord_protein_db_t *) * numBlocks);

    // Declare memory for compressed lookup table
    proteinLookup_db_cp = initialLookup =
        (struct initialWord_protein_db_cp *)global_malloc(
                sizeof(struct initialWord_protein_db_cp) * numEntries);

    //Declare memory for temp lookup table
    proteinLookup_db_temp = initialLookup_t =
        (struct initialWord_protein_db_t *)global_malloc(
                sizeof(struct initialWord_protein_db_t) * numEntries);

    proteinLookup_numWords = numEntries;

    // Iterate through every possible codeword
    codeword = 0;
    while (codeword < numBlocks) {
        // Initialize list of lookup table block 
        proteinLookup_db_blk_cp[codeword].proteinLookup_db_cp = initialLookup + codeword * numEntries;
        proteinLookup_db_blk_cp[codeword].dbIdxblk_longestSeq = 0;
        proteinLookup_db_blk_cp[codeword].numSeqBlk = 0;
        //proteinLookup_db_blk_temp[codeword] = initialLookup_t + codeword * numEntries;
        codeword++;
    }

    codeword = 0;
    while (codeword < numEntries) {
        // Initialize list of temp lookup table 
        initialLookup_t[codeword].numSubPositions = 0;
        initialLookup_t[codeword].allocSubPositions = 0;
        initialLookup_t[codeword].subSequencePositions = NULL;
        codeword++;
    }

    codeword = 0;
    while (codeword < numEntries) {
        // Initialize list of compressed lookup table 
        initialLookup[codeword].numPos = 0;
        initialLookup[codeword].numSubOff = 0;
        initialLookup[codeword].subOff = NULL;
        initialLookup[codeword].seqId = NULL;
        codeword++;
    }

}



void reinit_indexdb_temp() {
    uint4 codeword;
    codeword = 0;
    while (codeword < proteinLookup_numWords) {
        if (proteinLookup_db_temp[codeword].numSubPositions != 0) {
            free(proteinLookup_db_temp[codeword].subSequencePositions);
        }
        proteinLookup_db_temp[codeword].numSubPositions = 0;
        proteinLookup_db_temp[codeword].allocSubPositions = 0;
        proteinLookup_db_temp[codeword].subSequencePositions = NULL;
        codeword++;
    }
    //free(proteinLookup_db_blk_temp);
    //free(proteinLookup_db_temp);
}

void free_indexdb_cp() {
    uint4 codeword;
    codeword = 0;
    while (codeword < proteinLookup_numBlocks) {
        free(proteinLookup_db_blk_cp[codeword].proteinLookup_db_cp[0].subOff);
        free(proteinLookup_db_blk_cp[codeword].proteinLookup_db_cp[0].seqId);
        codeword++;
    }
    free(proteinLookup_db_blk_cp);
    free(proteinLookup_db_cp);
}

void reinit_indexdb_cp() {
    uint4 codeword;
    codeword = 0;
    while (codeword < proteinLookup_numWords) {
        free(proteinLookup_db_cp[codeword].subOff);
        free(proteinLookup_db_cp[codeword].seqId);
        // Initialize list of compressed lookup table 
        proteinLookup_db_cp[codeword].numPos = 0;
        proteinLookup_db_cp[codeword].numSubOff = 0;
        proteinLookup_db_cp[codeword].subOff = NULL;
        proteinLookup_db_cp[codeword].seqId = NULL;
        codeword++;
    }
}

void proteinLookup_db_cp_sub(uint4 sequenceNum, unsigned char *sequence,
        int4 subjectLength, int4 wordLength,
        struct scoreMatrix scoreMatrix, int startSeqNum) {

    uint4 codeword, numEntries, byteNumber, queryPosition;
    struct initialWord_protein_db_t *initialLookup, *initialWord;

    // Number of entries in the table
    numEntries = proteinLookup_numWords;

    initialLookup = proteinLookup_db_temp;

    // Slide a word-sized window across the query
    queryPosition = 0;

    while (queryPosition < subjectLength - wordLength + 1) {

        codeword = getCodeword(sequence + queryPosition,
                wordLength);

        initialWord = initialLookup + codeword;

        initialWord->numSubPositions++;
        if (initialWord->numSubPositions > initialWord->allocSubPositions) {
            initialWord->allocSubPositions += 50;
            //total_allocHits += 50;
            initialWord->subSequencePositions = (initialSubPos *)global_realloc(
                    initialWord->subSequencePositions,
                    sizeof(initialSubPos) * initialWord->allocSubPositions);
        }

        //initialWord->subSequencePositions[initialWord->numSubPositions - 1] =
        //((sequenceNum - blockNum * dbIdx_block_size) << 16) + queryPosition;

        initialWord->subSequencePositions[initialWord->numSubPositions - 1].subOff = queryPosition;
        initialWord->subSequencePositions[initialWord->numSubPositions - 1].seqId = sequenceNum - startSeqNum;
            //(queryPosition << 16) + (sequenceNum - startSeqNum);

        numPosBlkCP++;
        queryPosition++;
    }

}



void proteinLookup_db_cp_build(int4 numCodes, int wordLength,
        struct scoreMatrix scoreMatrix, char *filename_prefix) {
    proteinLookup_db_cp_initial(numCodes, wordLength);
    int ii, jj;

    printf("numBlocks: %d numWords: %d\n", proteinLookup_numBlocks,
            proteinLookup_numWords);

    numPosDBCP = 0;
    int seqStartBlk = 0;
    int numSplit = 0;
    long totalNumSubOff = 0;

    for (jj = 0; jj < proteinLookup_numBlocks; jj++) {
        int numLetterBlk = 0;
        numPosBlkCP = 0;
        int numSubOffBlk = 0;
        printf("buiding index block %d of %d ... ", jj, proteinLookup_numBlocks);
        proteinLookup_db_blk_cp[jj].seqOffset = seqStartBlk;
        for (ii = seqStartBlk; ii < readdb_numberOfSequences && numLetterBlk < dbIdx_block_size; ii++) {
            numLetterBlk += readdb_sequenceData[ii].sequenceLength;
            proteinLookup_db_cp_sub(ii, readdb_sequenceData[ii].sequence,
                    readdb_sequenceData[ii].sequenceLength, wordLength,
                    scoreMatrix, seqStartBlk);
            proteinLookup_db_blk_cp[jj].dbIdxblk_longestSeq = 
                max(proteinLookup_db_blk_cp[jj].dbIdxblk_longestSeq, readdb_sequenceData[ii].sequenceLength);

            if((ii - seqStartBlk) >= (INT4_MAX - 1))
            {
                break;
            }
        }

        if(numLetterBlk >= dbIdx_block_size)
        {
            numLetterBlk -= readdb_sequenceData[ii - 1].sequenceLength;
        }

        //printf(" ( %d %d %d ) ", seqStartBlk, ii, numLetterBlk);

        proteinLookup_db_blk_cp[jj].numSeqBlk = ii - seqStartBlk;

        //if(proteinLookup_db_blk_cp[jj].numSeqBlk >= INT2_MAX)
        //{
            //fprintf(stderr, "Err: too many sequences (%d) in a block, please use a smaller block size.\n", proteinLookup_db_blk_cp[jj].numSeqBlk);
            //exit(1);
        //}

        seqStartBlk = ii;
        numPosDBCP += numPosBlkCP;
        maxNumSeqBlk = max(proteinLookup_db_blk_cp[jj].numSeqBlk, maxNumSeqBlk);

        for (ii = 0; ii < proteinLookup_numWords; ii++)
        {
            if(proteinLookup_db_temp[ii].numSubPositions == 0)
                continue;

            qsort(proteinLookup_db_temp[ii].subSequencePositions, proteinLookup_db_temp[ii].numSubPositions, sizeof(initialSubPos), compareSubPos);
            int kk;
            int prev_pos = proteinLookup_db_temp[ii].subSequencePositions[0].subOff;
            int prev_none_zero_pos = 0;
            int same_pos_cnt = 0;
            int none_zero_pos_cnt = 0;
            ASSERT(proteinLookup_db_blk_cp[jj].dbIdxblk_longestSeq < INT4_MAX);
            proteinLookup_db_cp[ii].seqId = (uint2 *)global_malloc(sizeof(uint2) * 
                    (proteinLookup_db_temp[ii].numSubPositions));
            int numPosEst = proteinLookup_db_blk_cp[jj].dbIdxblk_longestSeq + 200;
            proteinLookup_db_cp[ii].subOff = (subOff_t *)global_malloc(sizeof(subOff_t) * numPosEst);
            //proteinLookup_db_blk_cp[jj].proteinLookup_db_cp[ii].posNum = (byte *)global_malloc(sizeof(byte) * proteinLookup_db_blk_cp[jj].dbIdxblk_longestSeq);

            for(kk = 0; kk < proteinLookup_db_temp[ii].numSubPositions + 1; kk++)
            {
                int currPos = kk < proteinLookup_db_temp[ii].numSubPositions 
                    ? (proteinLookup_db_temp[ii].subSequencePositions[kk].subOff) : -1;
                if(prev_pos == currPos)
                {
                    same_pos_cnt++;

                    if(same_pos_cnt > CHAR_MAX)
                    {
                        //printf("warning: same_pos_cnt: %d > CHAR_MAX, split group\n", same_pos_cnt);
                        proteinLookup_db_cp[ii].subOff[none_zero_pos_cnt].subOff = prev_pos - prev_none_zero_pos;
                        proteinLookup_db_cp[ii].subOff[none_zero_pos_cnt].numPos = CHAR_MAX;
                        prev_none_zero_pos = prev_pos;
                        same_pos_cnt = 1;
                        numSplit++;
                        none_zero_pos_cnt++;
                    }
                }
                else
                {
                    while((prev_pos - prev_none_zero_pos) > CHAR_MAX)
                    {
                        //printf("warning: posNum diff: %d > CHAR_MAX, split group\n", (prev_pos - prev_none_zero_pos));
                        proteinLookup_db_cp[ii].subOff[none_zero_pos_cnt].subOff = CHAR_MAX;
                        proteinLookup_db_cp[ii].subOff[none_zero_pos_cnt].numPos = 0;
                        prev_none_zero_pos = prev_none_zero_pos + CHAR_MAX;
                        none_zero_pos_cnt++;
                        numSplit++;
                    }
                    proteinLookup_db_cp[ii].subOff[none_zero_pos_cnt].subOff = prev_pos - prev_none_zero_pos;
                    proteinLookup_db_cp[ii].subOff[none_zero_pos_cnt].numPos = same_pos_cnt;
                    prev_none_zero_pos = prev_pos;
                    prev_pos = proteinLookup_db_temp[ii].subSequencePositions[kk].subOff;
                    same_pos_cnt = 1;
                    none_zero_pos_cnt++;
                }

                if(kk < proteinLookup_db_temp[ii].numSubPositions)
                    proteinLookup_db_cp[ii].seqId[kk] = proteinLookup_db_temp[ii].subSequencePositions[kk].seqId;

                ASSERT(none_zero_pos_cnt < numPosEst);
            }


            proteinLookup_db_cp[ii].numSubOff = none_zero_pos_cnt;
            proteinLookup_db_cp[ii].numPos = proteinLookup_db_temp[ii].numSubPositions;
            numSubOffBlk += none_zero_pos_cnt;

        }

        //Verify the correctness of compressed index
#if 1
        for (ii = 0; ii < proteinLookup_numWords; ii++)
        {
            uint4 trueSubOff = 0;
            uint4 seqIdCnt = 0;
            int kk;
            for(kk = 0; kk < proteinLookup_db_cp[ii].numSubOff; kk++)
            {
                uint2 numPos =  proteinLookup_db_cp[ii].subOff[kk].numPos;
                uint2 subOff = proteinLookup_db_cp[ii].subOff[kk].subOff;
                trueSubOff += subOff;
                uint4 tt;
                for(tt = 0; tt < numPos; tt++)
                {
                    ASSERT(proteinLookup_db_cp[ii].seqId[seqIdCnt] < proteinLookup_db_blk_cp[jj].numSeqBlk);

                    if(trueSubOff != proteinLookup_db_temp[ii].subSequencePositions[seqIdCnt].subOff || 
                            proteinLookup_db_cp[ii].seqId[seqIdCnt] != proteinLookup_db_temp[ii].subSequencePositions[seqIdCnt].seqId)
                    {
                        fprintf(stderr, "Err: word: %d subOff: %d pos: %d compressed index: %d %d <> original index: %d %d\n", 
                                ii, subOff, seqIdCnt,
                                trueSubOff, proteinLookup_db_cp[ii].seqId[seqIdCnt], 
                                proteinLookup_db_temp[ii].subSequencePositions[seqIdCnt].subOff, 
                                proteinLookup_db_temp[ii].subSequencePositions[seqIdCnt].seqId);
                        fflush(stderr);
                        exit(1);
                    }

                    seqIdCnt++;
                }

                ASSERT(trueSubOff <= readdb_longestSequenceLength);
            }

            if(seqIdCnt != proteinLookup_db_temp[ii].numSubPositions)
            {
                fprintf(stderr, "Err: blk: %d word: %d numSubOff: %d seqIdCnt %d <> numPos %d\n", jj, ii, proteinLookup_db_cp[ii].numSubOff, seqIdCnt, proteinLookup_db_temp[ii].numSubPositions);
                fflush(stderr);
            }
        }
#endif

        if(numPosBlkCP > 0)
        {

            printf("numSeqs: %d longestSeq: %u numPos: %d blkSize: %lu (MB)\n",
                    proteinLookup_db_blk_cp[jj].numSeqBlk,
                    proteinLookup_db_blk_cp[jj].dbIdxblk_longestSeq, numPosBlkCP,
                    (sizeof(struct initialWord_protein_db_cp) * proteinLookup_numWords +
                     sizeof(uint2) * numPosBlkCP + numSubOffBlk * sizeof(subOff_t)) >> 20);
            write_dbIdxBlock_cp(filename_prefix, jj);
        }
        else
        {
            printf("empty blk, skip\n");
            proteinLookup_numBlocks = jj;
            reinit_indexdb_temp();
            reinit_indexdb_cp();
            break;
        }

        totalNumSubOff += numSubOffBlk;

        reinit_indexdb_temp();
        reinit_indexdb_cp();
    }

    long total_index_size = sizeof(struct initialWord_protein_db_cp) * proteinLookup_numWords * proteinLookup_numBlocks
        + sizeof(uint2) * numPosDBCP 
        + sizeof(subOff_t) * totalNumSubOff;

    printf("totalNumBlk: %d totalNumPos: %ld totalIndexSize: %ld (MB) numSubOff: %ld numSplit: %d\n", proteinLookup_numBlocks, numPosDBCP, total_index_size >> 20, totalNumSubOff, numSplit);

    write_dbLookupAux_cp(filename_prefix);

    free(proteinLookup_db_blk_cp);
    free(proteinLookup_db_cp);
    free(proteinLookup_db_temp);

    //proteinLookup_db_print();
}

void write_dbIdxBlock_cp(char *write_dbLookupFilename, int blockNum)
{
    char* write_dbIdxBlockFilename;
    FILE* write_dbIdxBlockFile;

    write_dbIdxBlockFilename = (char *)global_malloc(strlen(write_dbLookupFilename) + 20);
    sprintf(write_dbIdxBlockFilename, "%s.dbLookupCP.%d", write_dbLookupFilename, blockNum);
    if ((write_dbIdxBlockFile = fopen(write_dbIdxBlockFilename, "w")) == NULL)
    {
        fprintf(stderr, "Error opening file %s for writing\n", write_dbIdxBlockFilename);
        exit(-1);
    }

    int ii;
    if (fwrite(proteinLookup_db_cp, sizeof(struct initialWord_protein_db_cp), proteinLookup_numWords, write_dbIdxBlockFile) != proteinLookup_numWords)
    {
        fprintf(stderr, "Error writing header to dbIdxBlock file %s\n", write_dbIdxBlockFilename);
        exit(-1);
    }   

    for(ii = 0; ii < proteinLookup_numWords; ii++)
    {
        //byte *numSeqIds  = proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp[ii].numSeqIds;
        subOff_t *subOffs = proteinLookup_db_cp[ii].subOff;
        if (fwrite(subOffs, sizeof(subOff_t), proteinLookup_db_cp[ii].numSubOff, write_dbIdxBlockFile) != proteinLookup_db_cp[ii].numSubOff)
        {
            fprintf(stderr, "Error writing data to dbIdxBlock file %s\n", write_dbIdxBlockFilename);
            exit(-1);
        }   
    }

    //for(ii = 0; ii < proteinLookup_numWords; ii++)
    //{
        //byte *posNum  = proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp[ii].posNum;
        //if (fwrite(posNum, sizeof(byte), proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp[ii].numPosGroups, write_dbIdxBlockFile) != proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp[ii].numPosGroups)
        //{
            //fprintf(stderr, "Error writing data to dbIdxBlock file %s\n", write_dbIdxBlockFilename);
            //exit(-1);
        //}   
    //}


    for(ii = 0; ii < proteinLookup_numWords; ii++)
    {
        uint2* seqId = proteinLookup_db_cp[ii].seqId;
        if (fwrite(seqId, sizeof(uint2), proteinLookup_db_cp[ii].numPos, write_dbIdxBlockFile) != proteinLookup_db_cp[ii].numPos)
        {
            fprintf(stderr, "Error writing numSubPositions to dbIdxBlock file %s\n", write_dbIdxBlockFilename);
            exit(-1);
        }
    }

    fclose(write_dbIdxBlockFile);
    free(write_dbIdxBlockFilename);

}

void read_dbIdxBlock_cp(char *read_dbLookupFilename, int blockNum)
{
    char* read_dbIdxBlockFilename;
    FILE* read_dbIdxBlockFile;

    read_dbIdxBlockFilename = (char *)global_malloc(strlen(read_dbLookupFilename) + 20);
    sprintf(read_dbIdxBlockFilename, "%s.dbLookupCP.%d", read_dbLookupFilename, blockNum);
    if ((read_dbIdxBlockFile = fopen(read_dbIdxBlockFilename, "r")) == NULL)
    {
        fprintf(stderr, "Error opening file %s for reading\n", read_dbIdxBlockFilename);
        exit(-1);
    }

    if (fread(proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp, sizeof(struct initialWord_protein_db_cp), proteinLookup_numWords, read_dbIdxBlockFile) != proteinLookup_numWords)
    {
        fprintf(stderr, "Error reading LoopkupHeader to dbIdxBlock file %s\n", read_dbIdxBlockFilename);
        exit(-1);
    }

    int ii;
    int totalPos = 0;
    int totalSubOff = 0;
    for(ii = 0; ii < proteinLookup_numWords; ii++)
    {
        totalPos += proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp[ii].numPos;
        totalSubOff += proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp[ii].numSubOff;
        //printf("%d : %d\n", ii, proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp[ii].numPos);
    }

    subOff_t *subOffs = (subOff_t *)global_malloc(sizeof(subOff_t) * totalSubOff); 
    if (fread(subOffs, sizeof(subOff_t), totalSubOff, read_dbIdxBlockFile) != totalSubOff)
    {
        fprintf(stderr, "Error reading data from dbIdxBlock file %s\n", read_dbIdxBlockFilename);
        exit(-1);
    }


    //byte *posNum = (byte *)global_malloc(sizeof(byte) * totalNumSubOff); 
    //if (fread(posNum, sizeof(byte), totalNumSubOff, read_dbIdxBlockFile) != totalNumSubOff)
    //{
        //fprintf(stderr, "Error reading data from dbIdxBlock file %s\n", read_dbIdxBlockFilename);
        //exit(-1);
    //}

    uint2* seqId  = (uint2*)global_malloc(sizeof(uint2) * totalPos);
    if (fread(seqId, sizeof(uint2), totalPos, read_dbIdxBlockFile) != totalPos)
    {
        fprintf(stderr, "Error reading data from dbIdxBlock file %s\n", read_dbIdxBlockFilename);
        exit(-1);
    }


    totalPos = 0;
    totalSubOff = 0;
    for(ii = 0; ii < proteinLookup_numWords; ii++)
    {
        proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp[ii].seqId = seqId + totalPos;
        proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp[ii].subOff = subOffs + totalSubOff;
        //proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp[ii].posNum = posNum + totalNumSubOff;
        totalPos += proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp[ii].numPos;
        totalSubOff += proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp[ii].numSubOff;
    }

    //printf("totalHits: %d - %d (MB)\n", totalHits, sizeof(uint4) * totalHits >> 20);

    fclose(read_dbIdxBlockFile);
    free(read_dbIdxBlockFilename);
}

struct dbLookupAuxCP{
    uint4 proteinLookup_numBlocks;
    uint4 proteinLookup_numWords;
    uint4 proteinLookup_wordLength;
    int4 proteinLookup_numCodes;
    int4 dbIdx_block_size;
    int4 totalNumSubOffitions;
    int4 maxNumSeqBlk;
};

void write_dbLookupAux_cp(char *write_dbLookupFilename)
{
    char *write_dbLookupAuxCPFilename;
    FILE *write_dbLookupAuxCPFile;

    write_dbLookupAuxCPFilename = (char *)global_malloc(strlen(write_dbLookupFilename) + 20);
    sprintf(write_dbLookupAuxCPFilename, "%s.dbLookupAuxCP", write_dbLookupFilename);
    // Open dbLookup file for writing
    if ((write_dbLookupAuxCPFile = fopen(write_dbLookupAuxCPFilename, "w")) == NULL)
    {
        fprintf(stderr, "Error opening file %s for writing\n", write_dbLookupAuxCPFilename);
        exit(-1);
    }


    struct dbLookupAuxCP dbLookupAuxCP;
    dbLookupAuxCP.proteinLookup_numBlocks = proteinLookup_numBlocks;
    dbLookupAuxCP.proteinLookup_numWords = proteinLookup_numWords;
    dbLookupAuxCP.proteinLookup_wordLength = wordLookupDFA_wordLength;
    dbLookupAuxCP.proteinLookup_numCodes = wordLookupDFA_numCodes;
    dbLookupAuxCP.dbIdx_block_size = dbIdx_block_size;
    dbLookupAuxCP.totalNumSubOffitions = numPosDBCP;
    dbLookupAuxCP.maxNumSeqBlk = maxNumSeqBlk;

    if (fwrite(&dbLookupAuxCP, sizeof(struct dbLookupAuxCP), 1, write_dbLookupAuxCPFile) != 1)
    {
        fprintf(stderr, "Error writing data to dbLookup aux file %s\n", write_dbLookupAuxCPFilename);
        exit(-1);
    }

    if (fwrite(proteinLookup_db_blk_cp, sizeof(struct proteinLookup_db_blk_cp), proteinLookup_numBlocks, write_dbLookupAuxCPFile) != proteinLookup_numBlocks)
    {
        fprintf(stderr, "Error writing data to dbLookup aux file %s\n", write_dbLookupAuxCPFilename);
        exit(-1);
    }

    //if (fwrite(numSeqBlk, sizeof(int), proteinLookup_numBlocks, write_dbLookupAuxCPFile) != proteinLookup_numBlocks)
    //{
    //fprintf(stderr, "Error writing data to dbLookup aux file %s\n", write_dbLookupAuxCPFilename);
    //exit(-1);
    //}

    free(write_dbLookupAuxCPFilename);
    fclose(write_dbLookupAuxCPFile);
}

void read_dbLookupAux_cp(char *read_dbLookupFilename)
{
    char *read_dbLookupAuxCPFilename;
    FILE *read_dbLookupAuxCPFile;

    read_dbLookupAuxCPFilename = (char *)global_malloc(strlen(read_dbLookupFilename) + 20);
    sprintf(read_dbLookupAuxCPFilename, "%s.dbLookupAuxCP", read_dbLookupFilename);

    // Open dbLookup file for writing
    if ((read_dbLookupAuxCPFile = fopen(read_dbLookupAuxCPFilename, "r")) == NULL)
    {
        fprintf(stderr, "Error opening file %s for reading\n", read_dbLookupAuxCPFilename);
        exit(-1);
    }

    struct dbLookupAuxCP dbLookupAuxCP;

    fread(&dbLookupAuxCP, sizeof(struct dbLookupAuxCP), 1, read_dbLookupAuxCPFile);

    proteinLookup_numBlocks = dbLookupAuxCP.proteinLookup_numBlocks;
    proteinLookup_numWords = dbLookupAuxCP.proteinLookup_numWords;
    wordLookupDFA_wordLength = dbLookupAuxCP.proteinLookup_wordLength;
    wordLookupDFA_numCodes = dbLookupAuxCP.proteinLookup_numCodes;
    dbIdx_block_size = dbLookupAuxCP.dbIdx_block_size;
    numPosDBCP = dbLookupAuxCP.totalNumSubOffitions;
    maxNumSeqBlk = dbLookupAuxCP.maxNumSeqBlk;

    proteinLookup_db_blk_cp = (struct proteinLookup_db_blk_cp *)malloc(sizeof(struct proteinLookup_db_blk_cp) * proteinLookup_numBlocks);
    fread(proteinLookup_db_blk_cp, sizeof(struct proteinLookup_db_blk_cp), proteinLookup_numBlocks, read_dbLookupAuxCPFile);

    printf("dbIndex info:\nproteinLookup_numBlocks: %d\nproteinLookup_numWords: %d\n"
            "proteinLookup_wordLength: %d\nproteinLookup_numCodes: %d\nproteinLookup_dbIdxBlockSize: %d (K letters)\ntotalNumSubOffitions: %ld\nmaxNumSeqBlk: %d\n\n",
            proteinLookup_numBlocks, proteinLookup_numWords,
            wordLookupDFA_wordLength, wordLookupDFA_numCodes, dbIdx_block_size/1024, numPosDBCP, maxNumSeqBlk);
    fflush(stdout);

    proteinLookup_db_cp = (struct initialWord_protein_db_cp *)malloc(
            sizeof(struct initialWord_protein_db_cp) * proteinLookup_numWords *
            proteinLookup_numBlocks);

    int ii;
    for(ii = 0; ii < proteinLookup_numBlocks; ii++)
    {
        proteinLookup_db_blk_cp[ii].proteinLookup_db_cp = proteinLookup_db_cp + ii * proteinLookup_numWords;
        //printf("seqOffset: %d numSeq: %d\n", proteinLookup_db_blk_cp[ii].seqOffset, proteinLookup_db_blk_cp[ii].numSeqBlk);
        //printf("dbIdx %d longSeqLength: %d\n", ii, dbIdxblk_longestSeq[ii]);
    }

    free(read_dbLookupAuxCPFilename);
    fclose(read_dbLookupAuxCPFile);
}

void read_dbLookup_cp(char *read_dbLookupFilename)
{

#ifdef PROFILE
    struct timeval start, end;
    gettimeofday(&start, NULL);
#endif
    read_dbLookupAux_cp(read_dbLookupFilename); 

    int ii;
    for(ii = 0; ii < proteinLookup_numBlocks; ii++)
    {
        //printf("reading DB_BLOCK %d ", ii);
        read_dbIdxBlock_cp(read_dbLookupFilename, ii);
    }

#ifdef PROFILE
    gettimeofday(&end, NULL);
    long read_time = ((end.tv_sec * 1000000 + end.tv_usec)
            - (start.tv_sec * 1000000 + start.tv_usec));
    printf("Index read time: %f\n", (float)read_time * 1e-6);
#endif


}

