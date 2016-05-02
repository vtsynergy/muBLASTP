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


#ifndef _alignments_ncbi_
#define _alignments_ncbi_


int alignments_findGoodAlignments_ncbi(
        align_t *alignment,
        struct ungappedExtension *goodExtensionBuf,
        struct PSSMatrix PSSMatrix,
        struct scoreMatrix scoreMatrix,
        int queryNum,
        struct ungappedExtension **ungappedExtension_new,
        BlastGapDP *dp_mem,
        BlastIntervalTree *tree,
        BlastIntervalTree *private_tree, 
        BlastHSP *BlastHSP_arr); 

int alignments_findGoodAlignments_ncbi(
        align_t *alignment,
        struct ungappedExtension *goodExtensionBuf,
        struct PSSMatrix PSSMatrix,
        struct scoreMatrix scoreMatrix,
        int queryNum,
        struct ungappedExtension **ungappedExtension_new,
        BlastGapDP *dp_mem,
        BlastIntervalTree *tree,
        BlastIntervalTree *private_tree, 
        BlastHSP *BlastHSP_arr); 

#endif
