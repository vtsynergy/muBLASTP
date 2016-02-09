#ifndef _NCBI_STAT_
#define _NCBI_STAT_

#include <inttypes.h>

#define Int1 char 
#define Int4 int32_t
#define Uint4 uint32_t
#define Int8 int64_t
#define Uint8 uint64_t
#define Int2 int16_t
#define Uint2 uint16_t

/*typedef unsigned char Uint1;*/
/*typedef char Boolean;*/

typedef struct BlastHSPList {
   Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
   Int4 query_index; /**< Index of the query which this HSPList corresponds to.
                        Set to 0 if not applicable */
   BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
   Int4 hspcnt; /**< Number of HSPs saved */
   Int4 allocated; /**< The allocated size of the hsp_array */
   Int4 hsp_max; /**< The maximal number of HSPs allowed to be saved */
   Boolean do_not_reallocate; /**< Is reallocation of the hsp_array allowed? */
   double best_evalue; /**< Smallest e-value for HSPs in this list. Filled after 
                          e-values are calculated. Necessary because HSPs are
                          sorted by score, but highest scoring HSP may not have
                          the lowest e-value if sum statistics is used. */
} BlastHSPList;



Int2
BLAST_ScoreSetAmbigRes(BlastScoreBlk* sbp, char ambiguous_res);

Int2
Blast_ScoreBlkKbpUngappedCalc(
        BlastScoreBlk* sbp, 
        //int num_query,
        struct PSSMatrix PSSMatrix,
        Int4 *blast_ungappedNominalTrigger,
        Blast_KarlinBlk *kbp
        );

Int4 
BLAST_SpougeEtoS(double e0,
                 Blast_KarlinBlk* kbp,
                 Blast_GumbelBlk* gbp,
                 Int4 m, Int4 n); 

Int2 Blast_HSPListGetEvalues(
        //EBlastProgramType program_number,
        //const BlastQueryInfo* query_info,
        int query_length,
        Int4 subject_length,
        Blast_KarlinBlk* kbp,
        Blast_GumbelBlk* gbp,
        struct ungappedExtension **ungappedExtensions,
        /*struct alignment *alignment,*/
        Int4 hsp_cnt,
        //BlastHSPList* hsp_list, 
        Boolean gapped_calculation, 
        Boolean RPS_prelim,
        double gap_decay_rate,
        double scaling_factor,
        double *best_eValue);

Blast_KarlinBlk *Blast_KarlinBlkNew(void);

Blast_KarlinBlk *Blast_KarlinBlkFree(Blast_KarlinBlk *kbp);


Int2 
s_Blast_HSPListReapByPrelimEvalue(struct ungappedExtension **ungappedExtensions, int4 numExts, double prelim_evalue);

Int4
Blast_HSPListPurgeHSPsWithCommonEndpoints( 
        struct ungappedExtension **hsp_array,
        Int4 hsp_count,
        Boolean purge);


void**
_PSIDeallocateMatrix(void** matrix, unsigned int ncols);


double Blast_KarlinLambdaNR(Blast_ScoreFreq *sfp, double initialLambdaGuess);

void**
_PSIAllocateMatrix(unsigned int ncols, unsigned int nrows, 
                   unsigned int data_type_sz);


Int4 get_ncbi_dropoff_score(struct PSSMatrix PSSMatrix);


Int4 get_ncbi_ungappedNominalTrigger(struct PSSMatrix PSSMatrix, struct scoreMatrix scoreMatrix);

#endif
