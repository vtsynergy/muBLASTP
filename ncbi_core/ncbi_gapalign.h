#ifndef _NCBI_GAPALIGN_
#define _NCBI_GAPALIGN_

#define TRUE 1
#define FALSE 0

#define Int4 int32_t
#define Uint4 uint32_t
#define Int8 int64_t
#define Uint8 uint64_t
#define Int2 int16_t
#define Uint2 uint16_t



typedef int4 Int4;    
typedef uint4 Uint4;    
typedef int2 Int2;    
typedef unsigned char Uint1;
typedef char Boolean;
/** Scoring parameters block
 *  Contains scoring-related information that is actually used
 *  for the blast search
 */
typedef struct BlastScoringParameters {
   Int2 reward;      /**< Reward for a match */
   Int2 penalty;     /**< Penalty for a mismatch */
   Int4 gap_open;    /**< Extra penalty for starting a gap (scaled version) */
   Int4 gap_extend;  /**< Penalty for each gap residue  (scaled version) */
   Int4 shift_pen;   /**< Penalty for shifting a frame in out-of-frame 
                        gapping (scaled version) */
   double scale_factor; /**< multiplier for all cutoff scores */
} BlastScoringParameters;

/** Auxiliary structure for dynamic programming gapped extension */
typedef struct {
  Int4 best;            /**< score of best path that ends in a match
                             at this position */
  Int4 best_gap;        /**< score of best path that ends in a gap
                             at this position */
} BlastGapDP;


typedef struct BlastSeg {
   Int2 frame;  /**< Translation frame */
   Int4 offset; /**< Start of hsp */
   Int4 end;    /**< End of hsp */
   Int4 gapped_start;/**< Where the gapped extension started. */
} BlastSeg;




/** Structure holding all information about an HSP */
typedef struct BlastHSP {
   Int4 score;           /**< This HSP's raw score */
   Int4 num_ident;       /**< Number of identical base pairs in this HSP */
   double bit_score;     /**< Bit score, calculated from score */
   double evalue;        /**< This HSP's e-value */
   BlastSeg query;       /**< Query sequence info. */
   BlastSeg subject;     /**< Subject sequence info. */
   Int4     context;     /**< Context number of query */
   /*GapEditScript* gap_info;*//**< ALL gapped alignment is here */
   Int4 num;             /**< How many HSP's are linked together for sum 
                              statistics evaluation? If unset (0), this HSP is
                              not part of a linked set, i.e. value 0 is treated
                              the same way as 1. */
   Int2     comp_adjustment_method;  /**< which mode of composition
                                              adjustment was used; relevant
                                              only for blastp and tblastn */
   /*   SPHIHspInfo* pat_info; *//**< In PHI BLAST, information about this pattern*/
   /*                                 match. */
   Int4 num_positives;
} BlastHSP;


/*typedef struct BlastScoringParameters {*/
/*BlastScoringOptions *options; *//**< User-provided values for these params */
/*Int2 reward;      *//**< Reward for a match */
/*Int2 penalty;     *//**< Penalty for a mismatch */
/*Int4 gap_open;    *//**< Extra penalty for starting a gap (scaled version) */
/*Int4 gap_extend;  *//**< Penalty for each gap residue  (scaled version) */
/*   Int4 shift_pen;   *//**< Penalty for shifting a frame in out-of-frame */
/*                        gapping (scaled version) */
/*double scale_factor; *//**< multiplier for all cutoff scores */
/*} BlastScoringParameters;*/
typedef enum EGapAlignOpType { 
   eGapAlignDel = 0, /**< Deletion: a gap in query */
   eGapAlignDel2 = 1,/**< Frame shift deletion of two nucleotides */
   eGapAlignDel1 = 2,/**< Frame shift deletion of one nucleotide */
   eGapAlignSub = 3, /**< Substitution */
   eGapAlignIns1 = 4,/**< Frame shift insertion of one nucleotide */
   eGapAlignIns2 = 5,/**< Frame shift insertion of two nucleotides */
   eGapAlignIns = 6, /**< Insertion: a gap in subject */
   eGapAlignDecline = 7, /**< Non-aligned region */
   eGapAlignInvalid = 8 /**< Invalid operation */
} EGapAlignOpType;

/** A version of GapEditScript used to store initial results
    from the gapped alignment routines */
typedef struct GapPrelimEditScript {
   EGapAlignOpType op_type;    /**< Type of operation */
   Int4 num;                   /**< Number of operations */
} GapPrelimEditScript;

typedef struct Blast_ScoreFreq {
    Int4         score_min; /**< lowest allowed scores */
    Int4         score_max; /**< highest allowed scores */
    Int4         obs_min;   /**< lowest observed (actual) scores */
    Int4         obs_max;   /**< highest observed (actual) scores */
    double       score_avg; /**< average score, must be negative for local alignment. */
    double*      sprob0;    /**< arrays for frequency of given score */
    double*      sprob;     /**< arrays for frequency of given score, shifted down by score_min. */
} Blast_ScoreFreq;



typedef struct Blast_KarlinBlk {
      double   Lambda; /**< Lambda value used in statistics */
      double   K; /**< K value used in statistics */
      double   logK; /**< natural log of K value used in statistics */
      double   H; /**< H value used in statistics */
      double   paramC;  /**< for use in seed. */
} Blast_KarlinBlk;

/**
  Structure to hold the Gumbel parameters (for FSC).
  */
typedef struct Blast_GumbelBlk {
    double  Lambda;    /**< the unscaled Lambda value */
    double  C;
    double  G;         /**< G is the total penalty for extension */
    double  a;         /**< avg(L) = a     y + b    */
    double  Alpha;     /**< var(L) = alpha y + beta */
    double  Sigma;     /**< cov(L) = sigma y + tau  */
    double  a_un;      /**< Ungapped a */
    double  Alpha_un;  /**< Ungapped alpha */

    double  b;         /**< 2*G*(a_un - a) */
    double  Beta;      /**< 2*G*(alpha_un - alpha) */
    double  Tau;       /**< 2*G*(alpha_un - Sigma) */

    Int8 db_length;    /**< total length of database */

    Boolean filled;    /**< flag indicate the values of gbp are prepared */

} Blast_GumbelBlk;



/** Structure to keep memory for state structure. */ 
typedef struct GapStateArrayStruct {
    Int4    length,     /**< length of the state_array. */
        used;       /**< how much of length is used. */
    Uint1* state_array; /**< array to be used. */
    struct GapStateArrayStruct* next; /**< Next link in the list. */
} GapStateArrayStruct;

/** Edit script: linked list of correspondencies between two sequences */
typedef struct GapEditScript {
   EGapAlignOpType* op_type;    /**< Array of type of operation */
   Int4* num;                   /**< Array of number of operations */
   Int4 size;                   /**< Size of above arrays. */
} GapEditScript;

typedef struct GapPrelimEditBlock {
    GapPrelimEditScript *edit_ops;  /**< array of edit operations */
    Int4 num_ops_allocated;        /**< size of allocated array */
    Int4 num_ops;                  /**< number of edit ops presently in use */
    EGapAlignOpType last_op;        /**< most recent operation added */
    Int4 total_num_ops;
} GapPrelimEditBlock;

typedef struct BlastScoreBlk {
   Boolean     protein_alphabet; /**< TRUE if alphabet_code is for a 
protein alphabet (e.g., ncbistdaa etc.), FALSE for nt. alphabets. */
   Uint1    alphabet_code; /**< NCBI alphabet code. */
   Int2     alphabet_size;  /**< size of alphabet. */
   Int2     alphabet_start;  /**< numerical value of 1st letter. */
   char*    name;           /**< name of scoring matrix. */
   Int2 **matrix;
   /*ListNode*   comments;    *//**< Comments about scoring matrix. */
   /*SBlastScoreMatrix* matrix;   *//**< scoring matrix data */
   /*SPsiBlastScoreMatrix* psi_matrix;    *//**< PSSM and associated data. If this*/
   /*is not NULL, then the BLAST search is*/
   /*                                         position specific (i.e.: PSI-BLAST) */
   Boolean  matrix_only_scoring;  /**< Score ungapped/gapped alignment only
                                       using the matrix parameters and
                                       with raw scores. Ignore 
                                       penalty/reward and do not report 
                                       Karlin-Altschul stats.  This is used
                                       by the rmblastn program. -RMH- */
   Boolean complexity_adjusted_scoring; /**< Use cross_match-like complexity
                                           adjustment on raw scores. -RMH- */
   Int4  loscore;   /**< Min.  substitution scores */
   Int4  hiscore;   /**< Max. substitution scores */
   Int4  penalty;   /**< penalty for mismatch in blastn. */
   Int4  reward;    /**< reward for match in blastn. */
        double  scale_factor; /**< multiplier for all cutoff and dropoff scores */
   Boolean     read_in_matrix; /**< If TRUE, matrix is read in, otherwise
               produce one from penalty and reward above. @todo should this be
                an allowed way of specifying the matrix to use? */
   Blast_ScoreFreq** sfp;  /**< score frequencies for scoring matrix. */
   /* kbp & kbp_gap are ptrs that should be set to kbp_std, kbp_psi, etc. */
   Blast_KarlinBlk** kbp;  /**< Karlin-Altschul parameters. Actually just a placeholder. */
   Blast_KarlinBlk** kbp_gap; /**< K-A parameters for gapped alignments.  Actually just a placeholder. */
   /*Blast_GumbelBlk* gbp;  *//**< Gumbel parameters for FSC. */
   /* Below are the Karlin-Altschul parameters for non-position based ('std')
   and position based ('psi') searches. */
   Blast_KarlinBlk **kbp_std,  /**< K-A parameters for ungapped alignments */
                    **kbp_psi,       /**< K-A parameters for position-based alignments. */
                    **kbp_gap_std,  /**< K-A parameters for std (not position-based) alignments */
                    **kbp_gap_psi;  /**< K-A parameters for psi alignments. */
   Blast_KarlinBlk*  kbp_ideal;  /**< Ideal values (for query with average database composition). */
   Int4 number_of_contexts;   /**< Used by sfp and kbp, how large are these*/
   Uint1*   ambiguous_res; /**< Array of ambiguous res. (e.g, 'X', 'N')*/
   Int2     ambig_size, /**< size of array above. FIXME: not needed here? */
         ambig_occupy;  /**< How many occupied? */
   Boolean  round_down; /**< Score must be rounded down to nearest even score if odd. */
} BlastScoreBlk;


/* supporting the gapped alignment */
typedef struct BlastGapAlignStruct {
    Boolean positionBased; /**< Is this PSI-BLAST? */
    GapStateArrayStruct* state_struct; /**< Structure to keep extension */
    GapEditScript* edit_script; /**< The traceback (gap) information */
    GapPrelimEditBlock *fwd_prelim_tback; /**< traceback from right extensions */
    GapPrelimEditBlock *rev_prelim_tback; /**< traceback from left extensions */
    /*   SGreedyAlignMem* greedy_align_mem;*//**< Preallocated memory for the greedy */
    /*                                         gapped extension */
    BlastGapDP* dp_mem; /**< scratch structures for dynamic programming */
    Int4 dp_mem_alloc;  /**< current number of structures allocated */
    BlastScoreBlk* sbp; /**< Pointer to the scoring information block */
    Int4 gap_x_dropoff; /**< X-dropoff parameter to use */
    Int4 query_start; /**< query start offset of current alignment */
    Int4 query_stop; /**< query end offseet of current alignment */
    Int4 subject_start;  /**< subject start offset current alignment */
    Int4 subject_stop; /**< subject end offset of current alignment */
    Int4 greedy_query_seed_start;  /**< for greedy alignments, the query */
    Int4 greedy_subject_seed_start;  /**< for greedy alignments, the subject*/
    Int4 score;   /**< Return value: alignment score */
} BlastGapAlignStruct;

/** Structure to hold a sequence. */
typedef struct BLAST_SequenceBlk {
   Uint1* sequence; /**< Sequence used for search (could be translation). */
   Uint1* sequence_start; /**< Start of sequence, usually one byte before 
                               sequence as that byte is a NULL sentinel byte.*/
   Int4 length;         /**< Length of sequence. */
   Int2 frame; /**< Frame of the query, needed for translated searches */
   Int2 subject_strand; /**< Strand of the subject sequence for translated searches. 
                          Uses the same values as ENa_strand. */
   Int4 oid; /**< The ordinal id of the current sequence */
   Boolean sequence_allocated; /**< TRUE if memory has been allocated for 
                                  sequence */
   Boolean sequence_start_allocated; /**< TRUE if memory has been allocated 
                                        for sequence_start */
   Uint1* sequence_start_nomask; /**< Query sequence without masking. */
   Uint1* sequence_nomask; /**< Start of query sequence without masking. */
   Boolean nomask_allocated; /**< If false the two above are just pointers to
                                   sequence and sequence_start. */
   Uint1* oof_sequence; /**< Mixed-frame protein representation of a
                             nucleotide sequence for out-of-frame alignment */
   Boolean oof_sequence_allocated; /**< TRUE if memory has been allocated 
                                        for oof_sequence */
   Uint1* compressed_nuc_seq; /**< 4-to-1 compressed version of sequence */
   Uint1* compressed_nuc_seq_start; /**< start of compressed_nuc_seq */
   Boolean lcase_mask_allocated; /**< TRUE if memory has been allocated for 
                                    lcase_mask */
   Int4 chunk;  /**< Used for indexing only: the chunk number within the 
                     subject sequence. */
   Uint1 *gen_code_string;  /**< for nucleotide subject sequences (tblast[nx]),
                              the genetic code used to create a translated
                              protein sequence (NULL if not applicable). This
                              field is NOT owned by this data structure, it's
                              owned by the genetic code singleton. 
                              @sa gencode_singleton.h
                              */
   /* BEGIN: Data members needed for masking subjects from a BLAST database */
   Uint4 num_seq_ranges;    /**< Number of elements in seq_ranges */
   Boolean seq_ranges_allocated;   /**< TRUE if memory has been allocated for
                                      seq_ranges */
   Uint1 bases_offset; /* Bases offset in first byte for SRA seq */
} BLAST_SequenceBlk;

/** Structure holding a pair of offsets. Used for storing offsets for the
 * initial seeds. In most programs the offsets are query offset and subject 
 * offset of an initial word match. For PHI BLAST, the offsets are start and 
 * end of the pattern occurrence in subject, with no query information, 
 * because all pattern occurrences in subjects are aligned to all pattern 
 * occurrences in query.
 */
typedef union BlastOffsetPair {
    struct {
        Uint4 q_off;  /**< Query offset */
        Uint4 s_off;  /**< Subject offset */
    } qs_offsets;     /**< Query/subject offset pair */
    struct {
        Uint4 s_start;/**< Start offset of pattern in subject */
        Uint4 s_end;  /**< End offset of pattern in subject */
    } phi_offsets;    /**< Pattern offsets in subject (PHI BLAST only) */
} BlastOffsetPair;

/** Structure to hold ungapped alignment information */
typedef struct BlastUngappedData {
   Int4 q_start; /**< Start of the ungapped alignment in query */
   Int4 s_start; /**< Start of the ungapped alignment in subject */ 
   Int4 length;  /**< Length of the ungapped alignment */
   Int4 score;   /**< Score of the ungapped alignment */
} BlastUngappedData;

/** Structure to hold the initial HSP information */
typedef struct BlastInitHSP {
    BlastOffsetPair offsets; /**< Offsets in query and subject, or, in PHI
                                BLAST, start and end of pattern in subject. */
    BlastUngappedData* ungapped_data; /**< Pointer to a structure holding
                                         ungapped alignment information */
} BlastInitHSP;


Int2 s_BlastProtGappedAlignment(
   BLAST_SequenceBlk* query_blk, 
   Int4 q_off,
   BLAST_SequenceBlk* subject_blk, 
   Int4 s_off,
   BlastGapAlignStruct* gap_align,
   const BlastScoringParameters* score_params,
   Boolean restricted_alignment,
   struct ungappedExtension *ungappedExtension,
   struct scoreMatrix scoreMatrix);


Int4 
BlastGetStartForGappedAlignment (const Uint1* query, const Uint1* subject,
   int2 **matrix, Uint4 q_start, Uint4 q_length, 
   Uint4 s_start, Uint4 s_length);


BlastHSP* Blast_HSPNew(void);

BlastHSP* Blast_HSPInit(BlastHSP *Blast_HSP_arr, int *Blast_HSP_cnt);
#endif
