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



#ifndef _unpack_multi_
#define _unpack_multi_

void unpack_initialize_multi(int queryNum);

void unpack_unpackSubject_multi(struct PSSMatrix PSSMatrix, struct alignment* alignment, int queryNum);

void unpack_free_multi(int queryNum);

int4 unpack_loadSubject_multi(struct PSSMatrix PSSMatrix, struct alignment* alignment, int queryNum);

#endif
