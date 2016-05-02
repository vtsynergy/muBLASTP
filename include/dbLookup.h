#ifndef _dbLookup_
#define _dbLookup_

typedef struct{
    uint4 subOff:16;
    uint2 seqId:16;
} subPos_t;

struct initialWord_protein_db
{
    uint4  numSubPositions;
    uint4  allocSubPositions;
    /*uint4 subPositionOffset;*/
    subPos_t* subSequencePositions;
};


struct proteinLookup_db_blk
{
    int4 dbIdxblk_longestSeq;
    int4 numSeqBlk;
    int4 seqOffset;
    subPos_t* subSequencePositions;
    uint4 *subPositionOffset;
    /*struct initialWord_protein_db* proteinLookup_db;*/
};

   
struct initialWord_neighborLookup
{
	int4  numNeighbours;
    int4 *neighbours;
};

void read_dbLookupAux(char *read_dbLookupFilename);
extern struct initialWord_protein_db *proteinLookup_db;
extern struct proteinLookup_db_blk *proteinLookup_db_b;
extern struct initialWord_neighborLookup* neighborLookup;

extern uint4 proteinLookup_numWords;
extern uint4 proteinLookup_numBlocks;
extern uint4 maxNumSeqBlk;



void free_dbIdxBlock(int bid);

void free_dbIdxAux();

void proteinLookup_db_build(int4 numCodes, int wordLength,
        struct scoreMatrix scoreMatrix, char *filename);

/*void write_dbLookup(char *write_dbLookupFilename);*/
void read_dbLookup(char *write_dbLookupFilename);
/*void read_dbIdxBlock(char *read_dbLookupFilename, int blockNum);*/

void neighbourLookup_init();

void neighbourLookup_free();

void neighbourLookup_build(struct PSSMatrix PSSMatrix, struct scoreMatrix, int4 wordLength);

void free_dbindex();
#endif
