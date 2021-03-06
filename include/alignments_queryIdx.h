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


#ifndef _alignments_queryIdx_
#define _alignments_queryIdx_


void alignments_queryIdx(
        struct PSSMatrix *PSSMatrix_arr,
        struct scoreMatrix scoreMatrix, 
        char *query_arr[],
        char *queryDescription_arr[],
        int numQuery
        ); 


void alignments_queryIdx2(
        struct PSSMatrix *PSSMatrix_arr, 
        struct scoreMatrix scoreMatrix,
        char **query_arr,
        char **queryDescription_arr,
        int numQuery 
        ); 
#endif
