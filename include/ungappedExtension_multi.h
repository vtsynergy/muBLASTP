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



#ifndef _ungappedExtension_multi_
#define _ungappedExtension_multi_

/*extern unsigned char* ungappedExtension_subjectEndReached_multi[];*/
/*extern int4 ungappedExtension_bestScore_multi[];*/
extern struct memBlocks *ungappedExtension_extensions_multi[];
extern struct memBlocks ***ungappedExtension_extensions_multi2;

void ungappedExtension_initialize_multi();
void ungappedExtension_initialize_multi2();

struct ungappedExtension* ungappedExtension_extend_multi(int2** queryHit, unsigned char* subjectHit,
                        unsigned char* lastHit, struct PSSMatrix PSSMatrix, unsigned char* subject, int queryNum, unsigned char **ungappedExtension_subjectEndReached);


struct ungappedExtension* ungappedExtension_extend_multi2(int2** queryHit, unsigned char* subjectHit,
                        unsigned char* lastHit, struct PSSMatrix PSSMatrix, unsigned char* subject, int queryNum, int blockNum, unsigned char **ungappedExtension_subjectEndReached);

struct ungappedExtension* ungappedExtension_extend_ncbi_multi2(
        struct PSSMatrix PSSMatrix, struct scoreMatrix scoreMatrix, 
        char *subject, int4 lastHitOffset, 
        int4 subjectOffset, int4 queryOffset, 
        int subjectLength, int queryLength, 
        int4 sequenceCount, unsigned char **ungappedExtension_subjectEndReached, int queryNum, 
    struct ungappedExtension *goodExtensionBuf, int *goodExtensionCount,
 char *rightExtend);


int ungappedExtension_extend_multi3(
                        struct ungappedExtension *newUngappedExtension,
                        int2** queryHit, unsigned char* subjectHit,
                        unsigned char* lastHit, struct PSSMatrix PSSMatrix, unsigned char* subject, int queryNum, int blockNum, unsigned char **ungappedExtension_subjectEndReached);

int ungappedExtension_extend_ncbi_multi3(
                        struct ungappedExtension *newUngappedExtension,
                        struct PSSMatrix PSSMatrix, struct scoreMatrix scoreMatrix, 
                        char *subject, int4 lastHitOffset, 
                        int4 subjectOffset, int4 queryOffset, 
                        int subjectLength, int queryLength, 
                        int4 sequenceCount, unsigned char **ungappedExtension_subjectEndReached, int queryNum, int blockNum, char *rightExtend);


#endif
