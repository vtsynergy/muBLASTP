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



#ifndef _alignments_multi_
#define _alignments_multi_

// All alignments with an ungapped extension above gapping threshold
extern struct memBlocks* alignments_alignments_multi[];

// Good alignments which contain a semi-gapped extension above cutoff
extern struct memSingleBlock* alignments_goodAlignments_multi[];

// Final alignments which contain a gapped extension above cutoff
extern struct memSingleBlock* alignments_finalAlignments_multi[];
extern struct alignment* alignments_currentAlignment_multi[];

void alignments_initialize_multi();
void alignments_initialize_multi2();

void alignments_createNew_multi(uint4 descriptionLocation, uint4 descriptionLength,
                          unsigned char* subject, int4 subjectLength,
                          int4 encodedLength, int queryNum);

struct alignment *alignments_createNew_multi_db(uint4 descriptionLocation, uint4 descriptionLength,
                          unsigned char* subject, int4 subjectLength,
                          int4 encodedLength, int queryNum, uint4 sequenceCount);

struct alignment *alignments_createNew_multi_db2(
    struct alignment *goodAlignBuf, int *goodAlignCount,
        uint4 descriptionLocation, uint4 descriptionLength,
                          unsigned char* subject, int4 subjectLength,
                          int4 encodedLength, uint4 sequenceCount);

void alignments_addUngappedExtension_multi(struct ungappedExtension* newExtension, int queryNum);

void alignments_addUngappedExtension_multi_db(struct ungappedExtension* newExtension, int queryNum,
        struct alignment *alignments_currentAlignment_multi_t);

void alignments_findGoodAlignments_multi(struct PSSMatrix PSSMatrix, int queryNum);

int alignments_findGoodAlignments_multi3(struct alignment *alignment, struct PSSMatrix PSSMatrix, int queryNum);

void alignments_findFinalAlignments_multi(struct PSSMatrix PSSMatrix, int queryNum);

void alignments_getFinalAlignmentDescriptions_multi(int queryNum);


void alignments_getTracebacks_multi(struct PSSMatrix PSSMatrix, int queryNum);

void alignments_free_multi(int queryNum);


void alignments_db_serial(struct PSSMatrix *PSSMatrix_arr,
    struct scoreMatrix scoreMatrix, char **query_arr, char **queryDescription_arr, int numQuery);

void alignments_db_omp(struct PSSMatrix *PSSMatrix_arr, struct scoreMatrix scoreMatrix, char **query_arr, char **queryDescription_arr, int numQuery);

void alignments_sortFinalAlignments_multi2(int queryNum);

void alignments_db_omp_neigbour(struct PSSMatrix *PSSMatrix_arr,
    struct scoreMatrix scoreMatrix, char **query_arr, char **queryDescription_arr, int numQuery);

void alignments_all_multi(struct PSSMatrix *PSSMatrix_arr, int numSequences, int numQuery);
void alignments_multi(struct PSSMatrix *PSSMatrix_arr, int numQuery);

void alignments_regularGappedAlignment_multi(struct PSSMatrix PSSMatrix,
	struct ungappedExtension* ungappedExtension, struct alignment* alignment, int queryNum);


void alignments_pruneRegion_multi(struct alignment* alignment, struct ungappedExtension* ungappedExtension, int queryNum);


void alignments_unpruneRegion_multi(struct alignment* alignment, struct ungappedExtension* oldUngappedExtension,
                              struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix, int queryNum);

void alignments_sortFinalAlignments_multi(int queryNum);

void alignments_global_free_multi2();


void alignments_free_multi2(int queryNum);

int alignments_findGoodAlignments_ncbi_multi3(struct alignment *alignment,
        struct PSSMatrix PSSMatrix,
        struct scoreMatrix scoreMatrix,
        int queryNum,
        struct ungappedExtension **ungappedExtension_new,
        BlastGapDP *dp_mem,
        BlastIntervalTree *tree,
        BlastIntervalTree *private_tree, 
        BlastHSP *BlastHSP_arr);


void alignments_query_serial(struct PSSMatrix *PSSMatrix_arr, 
        int numQuery, 
        struct scoreMatrix scoreMatrix);

struct finalAlignment *alignments_addFinalAlignment_multi(int4 highestNominalScore,
        struct alignment *alignment,
        int queryNum);

#endif
