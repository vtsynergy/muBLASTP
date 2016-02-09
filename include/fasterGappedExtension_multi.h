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


#ifndef _fasterGappedExtension_multi_
#define _fasterGappedExtension_multi_
struct gappedExtension* fasterGappedExtension_build_multi(struct ungappedExtension* ungappedExtension,
                                  struct PSSMatrix PSSMatrix, int4 subjectSize,
                                  unsigned char* subject, int4 dropoff, int queryNum);

void fasterGappedExtension_free_multi(int queryNum);

#endif
