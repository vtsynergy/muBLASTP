#ifndef _dbLookup_compress_
#define _dbLookup_compress_

typedef 
struct subOff_s{
    byte numPos;
    byte subOff;
}subOff_t;

struct initialWord_protein_db_cp
{
	uint4 numPos;
    uint2 numSubOff;
    subOff_t *subOff;
    uint2* seqId;
};

struct proteinLookup_db_blk_cp
{
    int dbIdxblk_longestSeq;
    int numSeqBlk;
    int seqOffset;
    struct initialWord_protein_db_cp* proteinLookup_db_cp;
};

extern struct initialWord_protein_db_cp *proteinLookup_db_cp;
extern struct proteinLookup_db_blk_cp *proteinLookup_db_blk_cp;

void proteinLookup_db_cp_build(int4 numCodes, int wordLength,
        struct scoreMatrix scoreMatrix, char *filename);

/*void write_dbLookup_cp(char *write_dbLookupFilename);*/
void read_dbLookup_cp(char *write_dbLookupFilename);
/*void read_dbIdxBlock_cp(char *read_dbLookupFilename, int blockNum);*/

void free_indexdb_cp();
#endif
