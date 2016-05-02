/** Parameters used to compute gapped alignments */
typedef struct BlastCompo_GappingParams {
    int gap_open;        /**< penalty for opening a gap */
    int gap_extend;      /**< penalty for extending a gapped alignment by
                           one residue */
    int decline_align;   /**< penalty for declining to align two characters */
    int x_dropoff;       /**< for x-drop algorithms, once a path falls below
                           the best score by this (positive) amount, the
                           path is no longer searched */
    void * context;      /**< a pointer to any additional gapping parameters
                           that may be needed by the calling routine. */
} BlastCompo_GappingParams;

typedef struct BlastScoringOptions {
   char* matrix;   /**< Name of the matrix containing all scores: needed for
                        finding neighboring words */
   char* matrix_path; /**< Directory path to where matrices are stored. */
   Int2 reward;      /**< Reward for a match */
   Int2 penalty;     /**< Penalty for a mismatch */
   Boolean gapped_calculation; /**< gap-free search if FALSE */
   Boolean complexity_adjusted_scoring; /**< Use cross_match-like complexity
                                           adjustment on raw scores. -RMH- */
   Int4 gap_open;    /**< Extra penalty for starting a gap */
   Int4 gap_extend;  /**< Penalty for each gap residue */

   /* only blastx and tblastn (When query & subj are diff) */
   Boolean is_ooframe; /**< Should out-of-frame gapping be used in a translated
                          search? */
   Int4 shift_pen;   /**< Penalty for shifting a frame in out-of-frame 
                        gapping */
} BlastScoringOptions;
/** Auxiliary structure for dynamic programming gapped extension */
/*typedef struct {*/
/*  Int4 best;            *//**< score of best path that ends in a match*/
/*                             at this position */
/*  Int4 best_gap;        *//**< score of best path that ends in a gap*/
/*                             at this position */
/*} BlastGapDP;*/

/* and data-structures needed to perform a gapped alignment */
typedef struct BlastKappa_GappingParamsContext {
    const BlastScoringParameters*  scoringParams;                /**< scoring parameters for a
                                        gapped alignment */
    BlastGapAlignStruct * gap_align;  /**< additional parameters for a
                                        gapped alignment */
    BlastScoreBlk* sbp;               /**< the score block for this search */
    double localScalingFactor;        /**< the amount by which this
                                        search has been scaled */
    /*    EBlastProgramType prog_number;    *//**< the type of search being*/
    /*                                        performed */
} BlastKappa_GappingParamsContext;

typedef enum EBlastPrelimGapExt {
    eDynProgScoreOnly,          /**< standard affine gapping */
    eGreedyScoreOnly,           /**< Greedy extension (megaBlast) */
    eSmithWatermanScoreOnly     /**< Score-only smith-waterman */
} EBlastPrelimGapExt;

typedef enum EBlastTbackExt {
    eDynProgTbck,          /**< standard affine gapping */
    eGreedyTbck,           /**< Greedy extension (megaBlast) */
    eSmithWatermanTbck,    /**< Smith-waterman finds optimal scores, then 
                                ALIGN_EX to find alignment. */
    eSmithWatermanTbckFull /**< Smith-waterman to find all alignments */
} EBlastTbackExt;

typedef struct BlastExtensionOptions {
   double gap_x_dropoff; /**< X-dropoff value for gapped extension (in bits) */
   double gap_x_dropoff_final;/**< X-dropoff value for the final gapped 
                                  extension (in bits) */
   EBlastPrelimGapExt ePrelimGapExt; /**< type of preliminary gapped extension (normally) for calculating
                              score. */
   EBlastTbackExt eTbackExt; /**< type of traceback extension. */
   Int4 compositionBasedStats; /**< mode of compositional adjustment to use;
                                   if zero then compositional adjustment is
                                   not used */
   Int4 unifiedP; /**< Indicates unified P values to be used in blastp or tblastn */
   /*EBlastProgramType program_number; *//**< indicates blastn, blastp, etc. */
} BlastExtensionOptions;

typedef struct BlastExtensionParameters {
   BlastExtensionOptions* options; /**< The original (unparsed) options. */
   Int4 gap_x_dropoff; /**< X-dropoff value for gapped extension (raw) */
   Int4 gap_x_dropoff_final;/**< X-dropoff value for the final gapped 
                               extension (raw) */
} BlastExtensionParameters;

Int2
BLAST_GapAlignStructNew(
        //const BlastExtensionParameters* ext_params, 
   Uint4 max_subject_length,
   //BlastScoreBlk* sbp, 
   BlastGapAlignStruct** gap_align_ptr);


BlastGapAlignStruct* 
BLAST_GapAlignStructFree(BlastGapAlignStruct* gap_align);


Int2 BLAST_GappedAlignmentWithTraceback( 
        const Uint1* query, const Uint1* subject, BlastGapAlignStruct* gap_align, 
        const BlastScoringParameters* score_params,
        Int4 q_start, Int4 s_start, Int4 query_length, Int4 subject_length,
        Boolean * fence_hit, Int4 **matrix );
