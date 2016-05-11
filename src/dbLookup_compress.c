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

void write_dbLookupAux_cp(char *write_dbLookupFilename);
void write_dbIdxBlock_cp(FILE *dbLookupFile_cp, int blockNum);

void proteinLookup_db_cp_initial(int4 numCodes, int wordLength) {
    struct initialWord_protein_db_cp *initialLookup, *initialWord;

    uint4 codeword, numEntries;

    wordLookupDFA_numCodes = numCodes;
    wordLookupDFA_wordLength = wordLength;
    numEntries = ceil(pow(numCodes, wordLength));
    int4 proteinLookup_wordLength = wordLength;

    uint4 numBlocks = ceil((float)readdb_numVolumeLetters / dbIdx_block_size) + 1;
    proteinLookup_numBlocks = numBlocks;

    // Declare memory for initial DB index blocks
    proteinLookup_db_blk_cp = (struct proteinLookup_db_blk_cp *)global_malloc(
            sizeof(struct proteinLookup_db_blk_cp) * numBlocks);

    proteinLookup_numWords = numEntries;

    // Declare memory for compressed lookup table
    proteinLookup_db_cp = (struct initialWord_protein_db_cp *)global_malloc(
            sizeof(struct initialWord_protein_db_cp) * proteinLookup_numWords);

    //Declare memory for temp lookup table
    proteinLookup_db_temp = 
        (struct initialWord_protein_db_t *)global_malloc(
                sizeof(struct initialWord_protein_db_t) * proteinLookup_numWords);


    // Iterate through every possible codeword
    codeword = 0;
    while (codeword < numBlocks) {
        // Initialize list of lookup table block 
        proteinLookup_db_blk_cp[codeword].dbIdxblk_longestSeq = 0;
        proteinLookup_db_blk_cp[codeword].numSeqBlk = 0;
        codeword++;
    }

    codeword = 0;
    while (codeword < numEntries) {
        // Initialize list of temp lookup table 
        proteinLookup_db_temp[codeword].numSubPositions = 0;
        proteinLookup_db_temp[codeword].allocSubPositions = 0;
        proteinLookup_db_temp[codeword].subSequencePositions = NULL;
        codeword++;
    }

    codeword = 0;
    while (codeword < numEntries) {
        // Initialize list of compressed lookup table 
        proteinLookup_db_cp[codeword].numPos = 0;
        proteinLookup_db_cp[codeword].numSubOff = 0;
        proteinLookup_db_cp[codeword].subOff = NULL;
        proteinLookup_db_cp[codeword].seqId = NULL;
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
        struct scoreMatrix scoreMatrix) {

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

        initialWord->subSequencePositions[initialWord->numSubPositions - 1].subOff = queryPosition;
        initialWord->subSequencePositions[initialWord->numSubPositions - 1].seqId = sequenceNum;

        numPosBlkCP++;
        queryPosition++;
    }

}



void proteinLookup_db_cp_build(int4 numCodes, int wordLength,
        struct scoreMatrix scoreMatrix, char *write_dbLookupFilename) {

    int ii, jj;
    
    FILE *dbLookupFile_cp;
    char *dbLookupFileName;

    dbLookupFileName = (char *)
        malloc(sizeof(char *) * strlen(write_dbLookupFilename) + 40);

    sprintf(dbLookupFileName, 
            "%s.sequence%d.dbLookup", write_dbLookupFilename, readdb_volume);

    if ((dbLookupFile_cp = fopen(dbLookupFileName, "w")) == NULL) {
        fprintf(stderr, "Error opening file %s for writing\n",
                dbLookupFileName);
        exit(-1);
    }

    free(dbLookupFileName);

    proteinLookup_db_cp_initial(numCodes, wordLength);

    fprintf(stderr, "volumn: %d seqOffset: %d volumNumSeq: %d numBlocks: %d numWords: %d\nindexing...", 
            readdb_volume,
            readdb_volumeOffset,
            readdb_numVolumeSequences,
            proteinLookup_numBlocks,
            proteinLookup_numWords);

    numPosDBCP = 0;
    int seqStartBlk = 0;
    int numSplit = 0;
    long totalNumSubOff = 0;

    for (jj = 0; jj < proteinLookup_numBlocks; jj++) {
        if(!(jj % 10))
        {
            fprintf(stderr, ".");
        }
        int numLetterBlk = 0;
        numPosBlkCP = 0;
        int numSubOffBlk = 0;
        proteinLookup_db_blk_cp[jj].seqOffset = seqStartBlk + readdb_volumeOffset;
        for (ii = seqStartBlk; ii < readdb_numVolumeSequences && numLetterBlk < dbIdx_block_size; ii++) {
            numLetterBlk += readdb_sequenceData[ii + readdb_volumeOffset].sequenceLength;
            proteinLookup_db_cp_sub(ii - seqStartBlk, readdb_sequenceData[ii + readdb_volumeOffset].sequence,
                    readdb_sequenceData[ii + readdb_volumeOffset].sequenceLength, wordLength,
                    scoreMatrix);
            proteinLookup_db_blk_cp[jj].dbIdxblk_longestSeq = 
                MAX(proteinLookup_db_blk_cp[jj].dbIdxblk_longestSeq, readdb_sequenceData[ii + readdb_volumeOffset].sequenceLength);

            if((ii - seqStartBlk) >= (INT4_MAX - 1))
            {
                break;
            }
        }

        if(numLetterBlk >= dbIdx_block_size)
        {
            numLetterBlk -= readdb_sequenceData[ii + readdb_volumeOffset - 1].sequenceLength;
            ii--;
        }

        proteinLookup_db_blk_cp[jj].numSeqBlk = ii - seqStartBlk + 1;

        seqStartBlk = ii;
        numPosDBCP += numPosBlkCP;
        maxNumSeqBlk = MAX(proteinLookup_db_blk_cp[jj].numSeqBlk, maxNumSeqBlk);

        for (ii = 0; ii < proteinLookup_numWords; ii++)
        {
            if(proteinLookup_db_temp[ii].numSubPositions == 0)
                continue;

            qsort(proteinLookup_db_temp[ii].subSequencePositions, 
                    proteinLookup_db_temp[ii].numSubPositions, 
                    sizeof(initialSubPos), compareSubPos);
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

            //if(prev_pos > CHAR_MAX)
            //{
                //fprintf(stderr, "prev_pos: %d\n", prev_pos);
            //}

            for(kk = 0; kk < proteinLookup_db_temp[ii].numSubPositions + 1; kk++)
            {
                int currPos = kk < proteinLookup_db_temp[ii].numSubPositions 
                    ? (proteinLookup_db_temp[ii].subSequencePositions[kk].subOff) : -1;

                if(prev_pos == currPos && (prev_pos - prev_none_zero_pos) <= CHAR_MAX)
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

                    if(trueSubOff != 
                            proteinLookup_db_temp[ii].subSequencePositions[seqIdCnt].subOff || 
                            proteinLookup_db_cp[ii].seqId[seqIdCnt] != 
                            proteinLookup_db_temp[ii].subSequencePositions[seqIdCnt].seqId)
                    {
                        fprintf(stderr, "\nErr (block: %d word: %d subOff: %d pos: %d): compressed index: %d %d <> original index: %d %d\n", 
                                jj,
                                ii, subOff, seqIdCnt,
                                trueSubOff, proteinLookup_db_cp[ii].seqId[seqIdCnt], 
                                proteinLookup_db_temp[ii].subSequencePositions[seqIdCnt].subOff, 
                                proteinLookup_db_temp[ii].subSequencePositions[seqIdCnt].seqId);
                        fflush(stderr);
                        //exit(1);
                    }

                    seqIdCnt++;
                }

                ASSERT(trueSubOff <= readdb_longestSequenceLength);
            }

            if(seqIdCnt != proteinLookup_db_temp[ii].numSubPositions)
            {
                fprintf(stderr, "\nErr: blk: %d word: %d numSubOff: %d seqIdCnt %d <> numPos %d\n", 
                        jj, ii, proteinLookup_db_cp[ii].numSubOff, seqIdCnt, proteinLookup_db_temp[ii].numSubPositions);
                fflush(stderr);
            }
        }
#endif

        totalNumSubOff += numSubOffBlk;


        if(numPosBlkCP > 0)
        {
            write_dbIdxBlock_cp(dbLookupFile_cp, jj);
            reinit_indexdb_temp();
            reinit_indexdb_cp();
        }
        else
        {
            //printf("skip empty blk %d\n", jj);
            proteinLookup_numBlocks = jj;
            reinit_indexdb_temp();
            reinit_indexdb_cp();
            break;
        }
    }

    fprintf(stderr, "\n");

    long total_index_size = sizeof(struct initialWord_protein_db_cp) * proteinLookup_numWords * proteinLookup_numBlocks
        + sizeof(uint2) * numPosDBCP 
        + sizeof(subOff_t) * totalNumSubOff;

    write_dbLookupAux_cp(write_dbLookupFilename);

    free(proteinLookup_db_blk_cp);
    free(proteinLookup_db_cp);
    free(proteinLookup_db_temp);

    fclose(dbLookupFile_cp);
}

void write_dbIdxBlock_cp(FILE *write_dbIdxBlockFile, int blockNum)
{
    int ii;
    if (fwrite(proteinLookup_db_cp, sizeof(struct initialWord_protein_db_cp), proteinLookup_numWords, write_dbIdxBlockFile) != proteinLookup_numWords)
    {
        fprintf(stderr, "Error writing header to dbIdxBlock file\n");
        exit(-1);
    }   

    for(ii = 0; ii < proteinLookup_numWords; ii++)
    {
        //byte *numSeqIds  = proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp[ii].numSeqIds;
        subOff_t *subOffs = proteinLookup_db_cp[ii].subOff;
        if (fwrite(subOffs, sizeof(subOff_t), proteinLookup_db_cp[ii].numSubOff, write_dbIdxBlockFile) != proteinLookup_db_cp[ii].numSubOff)
        {
            fprintf(stderr, "Error writing data to dbIdxBlock file\n");
            exit(-1);
        }   
    }


    for(ii = 0; ii < proteinLookup_numWords; ii++)
    {
        uint2* seqId = proteinLookup_db_cp[ii].seqId;
        if (fwrite(seqId, sizeof(uint2), proteinLookup_db_cp[ii].numPos, write_dbIdxBlockFile) != proteinLookup_db_cp[ii].numPos)
        {
            fprintf(stderr, "Error writing numSubPositions to dbIdxBlock file\n");
            exit(-1);
        }
    }

}

void read_dbIdxBlock_cp(FILE *read_dbIdxBlockFile, int blockNum)
{
    //char* read_dbIdxBlockFilename;
    //FILE* read_dbIdxBlockFile;

    //read_dbIdxBlockFilename = (char *)global_malloc(strlen(read_dbLookupFilename) + 20);
    //sprintf(read_dbIdxBlockFilename, "%s.dbLookupCP.%d", read_dbLookupFilename, blockNum);
    //if ((read_dbIdxBlockFile = fopen(read_dbIdxBlockFilename, "r")) == NULL)
    //{
        //fprintf(stderr, "Error opening file %s for reading\n", read_dbIdxBlockFilename);
        //exit(-1);
    //}

    if (fread(proteinLookup_db_blk_cp[blockNum].proteinLookup_db_cp, 
                sizeof(struct initialWord_protein_db_cp), proteinLookup_numWords, 
                read_dbIdxBlockFile) != proteinLookup_numWords)
    {
        fprintf(stderr, "Error reading LoopkupHeader to dbIdxBlock file: %d\n", blockNum);
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
        fprintf(stderr, "Error reading data from dbIdxBlock file\n");
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
        fprintf(stderr, "Error reading data from dbIdxBlock file\n");
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

    //fclose(read_dbIdxBlockFile);
    //free(read_dbIdxBlockFilename);
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

    write_dbLookupAuxCPFilename = (char *)global_malloc(strlen(write_dbLookupFilename) + 40);
    sprintf(write_dbLookupAuxCPFilename, "%s.sequence%d.dbLookupAuxCP", write_dbLookupFilename, readdb_volume);
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


#if 1
    fprintf(stderr,
            "dbIdxBlockSize:%d(KB) "
            "maxNumSeqPerBlk:%d longestSeqLength:%d\n",
            dbIdx_block_size / 1024, 
            maxNumSeqBlk, 
            readdb_longestSequenceLength);
#endif

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

    read_dbLookupAuxCPFilename = (char *)global_malloc(strlen(read_dbLookupFilename) + 40);
    sprintf(read_dbLookupAuxCPFilename, "%s.sequence%d.dbLookupAuxCP", read_dbLookupFilename, readdb_volume);

    // Open dbLookup file for writing
    if ((read_dbLookupAuxCPFile = fopen(read_dbLookupAuxCPFilename, "r")) == NULL)
    {
        fprintf(stderr, "Error opening file %s for reading\n", read_dbLookupAuxCPFilename);
        exit(-1);
    }

    struct dbLookupAuxCP dbLookupAuxCP;

    if(fread(&dbLookupAuxCP, sizeof(struct dbLookupAuxCP), 1, read_dbLookupAuxCPFile) != 1 )
    {
        fprintf(stderr, "fread read_dbLookupAuxCPFile failed!\n");
        exit(0);
    }

    proteinLookup_numBlocks = dbLookupAuxCP.proteinLookup_numBlocks;
    proteinLookup_numWords = dbLookupAuxCP.proteinLookup_numWords;
    wordLookupDFA_wordLength = dbLookupAuxCP.proteinLookup_wordLength;
    wordLookupDFA_numCodes = dbLookupAuxCP.proteinLookup_numCodes;
    dbIdx_block_size = dbLookupAuxCP.dbIdx_block_size;
    numPosDBCP = dbLookupAuxCP.totalNumSubOffitions;
    maxNumSeqBlk = dbLookupAuxCP.maxNumSeqBlk;

    proteinLookup_db_blk_cp = (struct proteinLookup_db_blk_cp *)malloc(sizeof(struct proteinLookup_db_blk_cp) * proteinLookup_numBlocks);
    if(fread(proteinLookup_db_blk_cp, sizeof(struct proteinLookup_db_blk_cp), proteinLookup_numBlocks, read_dbLookupAuxCPFile) != proteinLookup_numBlocks)
    {
        fprintf(stderr, "fread read_dbLookupAuxCPFile failed\n");
        exit(0);
    }

    //fprintf(stderr, "dbIndex info:\nproteinLookup_numBlocks: %d\nproteinLookup_numWords: %d\n"
            //"proteinLookup_wordLength: %d\nproteinLookup_numCodes: %d\nproteinLookup_dbIdxBlockSize: %d (K letters)\ntotalNumSubOffitions: %ld\nmaxNumSeqBlk: %d\n\n",
            //proteinLookup_numBlocks, proteinLookup_numWords,
            //wordLookupDFA_wordLength, wordLookupDFA_numCodes, dbIdx_block_size/1024, numPosDBCP, maxNumSeqBlk);
    //fflush(stdout);

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

    fprintf(stderr, "loading index(%d/%d)...", readdb_volume + 1, readdb_numberOfVolumes);

    struct timeval start, end;
    gettimeofday(&start, NULL);

    read_dbLookupAux_cp(read_dbLookupFilename); 

    char *dbLookupFileName;

    dbLookupFileName = (char *)
        malloc(sizeof(char *) * strlen(read_dbLookupFilename) + 40);

    sprintf(dbLookupFileName, 
            "%s.sequence%d.dbLookup", read_dbLookupFilename, readdb_volume);

    FILE *dbLookupFile_cp;

    if ((dbLookupFile_cp = fopen(dbLookupFileName, "r")) == NULL) {
        fprintf(stderr, "Error opening file %s for reading\n",
                dbLookupFileName);
        exit(-1);
    }

    free(dbLookupFileName);

    int ii;
    for(ii = 0; ii < proteinLookup_numBlocks; ii++)
    {
        //printf("reading DB_BLOCK %d ", ii);
        read_dbIdxBlock_cp(dbLookupFile_cp, ii);
    }

    gettimeofday(&end, NULL);
    long read_time = ((end.tv_sec * 1000000 + end.tv_usec)
            - (start.tv_sec * 1000000 + start.tv_usec));
    fprintf(stderr, "time: %f\n", (float)read_time * 1e-6);

    fclose(dbLookupFile_cp);
}

