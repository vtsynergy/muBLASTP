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



#ifndef _semiGappedScoring_multi_
#define _semiGappedScoring_multi_
// Perform semi-gapped alignment with restricted insertion
int4 semiGappedScoring_score_multi(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                        int4 subjectSize, unsigned char* subject, int4 dropoff, int queryNum);

void semiGappedScoring_free_multi(int queryNum);


extern int4 *semiGappedScoring_bestRow_multi[];
extern int4 *semiGappedScoring_insertQrow_multi[];
extern int4 semiGappedScoring_rowSizes_multi[];
#endif
