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



#ifndef _gappedExtension_multi_
#define _gappedExtension_multi_

// Build a gapped extension with a trace and nominal score from the seed point4 of an ungapped
// extension using dynamic programming
struct gappedExtension* gappedExtension_build_multi(struct ungappedExtension* ungappedExtension,
                        struct PSSMatrix PSSMatrix, int4 subjectSize, unsigned char* subject,
                        struct unpackRegion* unpackRegion, int4 dropoff, int queryNum);

void gappedExtension_free_multi(int queryNum);

void gappedExtension_score_multi(struct gappedExtension* gappedExtension, int queryNum);


void gappedExtension_score_ncbi_multi(struct gappedExtension *gappedExtension,
                                 int queryNum);

#endif

