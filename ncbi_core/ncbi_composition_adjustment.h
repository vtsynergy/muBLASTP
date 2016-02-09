#ifndef _ncbi_composition_adjustment_
#define _ncbi_composition_adjustment_


#define COMPO_LARGEST_ALPHABET 28

/**
 * Represents the composition of an amino-acid sequence, in the ncbistdaa
 * alphabet. */
typedef struct Blast_AminoAcidComposition {
    /** probabilities of each amino acid */
    double prob[COMPO_LARGEST_ALPHABET];
    int numTrueAminoAcids;   /**< number of true amino acids in the sequence,
                               omitting nonstandard amino acids */
} Blast_AminoAcidComposition;

typedef enum ECompoAdjustModes {
    /** Don't use composition based statistics */
    eNoCompositionBasedStats       = 0, 
    /** Composition-based statistics as in NAR 29:2994-3005, 2001 */
    eCompositionBasedStats         = 1, 
    /** Composition-based score adjustment as in Bioinformatics 21:902-911,
     * 2005, conditioned on sequence properties. Cannot be applied to PSSMs. */
    eCompositionMatrixAdjust       = 2, 
    /** Composition-based score adjustment as in Bioinformatics 21:902-911,
     * 2005, unconditionally. Cannot be applied to PSSMs. */
    eCompoForceFullMatrixAdjust    = 3,
    eNumCompoAdjustModes
} ECompoAdjustModes;

typedef struct Blast_CompositionWorkspace {
    double ** mat_b;       /**< joint probabilities for the matrix in
                             standard context */
    double ** mat_final;   /**< optimized target frequencies */

    double * first_standard_freq;     /**< background frequency vector
                                        of the first sequence */
    double * second_standard_freq;    /**< background frequency vector of
                                        the second sequence */
} Blast_CompositionWorkspace;

typedef enum EMatrixAdjustRule {
    eDontAdjustMatrix              = (-1),
    eCompoScaleOldMatrix           = 0,
    eUnconstrainedRelEntropy       = 1,
    eRelEntropyOldMatrixNewContext = 2,
    eRelEntropyOldMatrixOldContext = 3,
    eUserSpecifiedRelEntropy       = 4
} EMatrixAdjustRule;

#define COMPO_NUM_TRUE_AA 20
#define COMPO_LARGEST_ALPHABET 28

#define FSA_AA_SIZE 32

static int to_ncbi[FSA_AA_SIZE] = {
    11,
    1,
    7,
    17,
    19,
    5,
    18,
    10,
    4,
    14,
    9,
    16,
    13,
    15,
    6,
    22,
    12,
    8,
    3,
    20,
    2,
    23,
    21,
    24,
    25,
    26,
    27,
    0,
    0,
    0,
    0,
    0
};

static int to_fsa[COMPO_LARGEST_ALPHABET] = {0, 2, 21, 19, 9, 6, 15, 3, 18, 11, 8, 1, 17, 13, 10, 14, 12, 4, 7, 5, 20, 23, 16, 22, 24, 25, 26, 27};


/** The context related information */
typedef struct BlastContextInfo {
    Int4 query_offset;      /**< Offset of this query, strand or frame in the
                               concatenated super-query. */
    Int4 query_length;      /**< Length of this query, strand or frame */
    Int8 eff_searchsp;      /**< Effective search space for this context. */
    Int4 length_adjustment; /**< Length adjustment for boundary conditions */
    Int4 query_index;       /**< Index of query (same for all frames) */
    /*Int1 frame;             *//**< Frame number (-1, -2, -3, 0, 1, 2, or 3) */
    Boolean is_valid;       /**< Determine if this context is valid or not.
                              This field should be set only by the setup code
                              and read by subsequent stages of the BLAST search
                              */
} BlastContextInfo;

void
s_GetComposition(Blast_AminoAcidComposition * composition,
        int alphsize,
        Uint1 *data,
        int length
        );

int
Blast_AdjustScores(int ** matrix,
        const Blast_AminoAcidComposition * query_composition,
        int queryLength,
        const Blast_AminoAcidComposition * subject_composition,
        int subjectLength,
        const Blast_MatrixInfo * matrixInfo,
        ECompoAdjustModes composition_adjust_mode,
        int RE_pseudocounts,
        Blast_CompositionWorkspace *NRrecord,
        EMatrixAdjustRule *matrix_adjust_rule,
        /*double calc_lambda(double *,int,int,double),*/
        double *pvalueForThisPair,
        int compositionTestIndex,
        double *ratioToPassBack);

void
Blast_ReadAaComposition(Blast_AminoAcidComposition * composition,
        int alphsize,
        const Uint1 * sequence, int length);
void 
Blast_ReadAaComposition_fsa(Blast_AminoAcidComposition * composition,
        int alphsize,
        const Uint1 * sequence, int length);

double
s_CalcLambda(double probs[], int min_score, int max_score, double lambda0);

Blast_CompositionWorkspace * Blast_CompositionWorkspaceNew(void);

int
Blast_CompositionWorkspaceInit(Blast_CompositionWorkspace * NRrecord,
                               const char *matrixName);


void
Blast_CompositionWorkspaceFree(Blast_CompositionWorkspace ** pNRrecord);


void alignments_get_eValue(struct ungappedExtension **ungappedExtensions, int4 numExts, Int4 queryLength, Int4 subjectLength, double *best_eValue);


void finalAlignments_get_eValue(struct ungappedExtension **ungappedExtensions, int4 numExts, Int4 queryLength, Int4 subjectLength, double *best_eValue);

int 
s_HitlistEvaluateAndPurge(int * pbestScore, double *pbestEvalue,
                          struct ungappedExtension ** hsp_list,
                          int hspcnt,
                          int query_length,
                          int subject_length);

typedef struct ReNewtonSystem {
    int alphsize;              /**< the size of the alphabet */
    int constrain_rel_entropy; /**< if true, use the relative entropy
                                    constraint for this optimization
                                    problem */
    double ** W;               /**< A lower-triangular matrix
                                    representing a factorization of
                                    the (2,2) block, -J D^{-1} J^T, of
                                    the condensed linear system */
    double * Dinv;             /**< The diagonal elements of the
                                    inverse of the necessarily
                                    diagonal (1,1) block of the linear
                                    system */
    double * grad_re;          /**< the gradient of the
                                    relative-entropy constraint, if
                                    this constraint is used. */
} ReNewtonSystem;



int
Blast_AdjustScores2(Int4 ** matrix,
                   const Blast_AminoAcidComposition * query_composition,
                   int queryLength,
                   const Blast_AminoAcidComposition * subject_composition,
                   int subjectLength,
                   const Blast_MatrixInfo * matrixInfo,
                   ECompoAdjustModes composition_adjust_mode,
                   int RE_pseudocounts,
                   Blast_CompositionWorkspace *NRrecord,
                   EMatrixAdjustRule *matrix_adjust_rule,
                   //double calc_lambda(double *,int,int,double),
                   double *pvalueForThisPair,
                   int compositionTestIndex,
                   double *ratioToPassBack,
                   ReNewtonSystem *newton_system,
                   double *z, double *resids_x, double *resids_z, 
                   double *old_scores, double *workspace, double **grads, double **Scores);

void
Nlm_DenseMatrixFree(double *** mat);


void
ReNewtonSystemFree(ReNewtonSystem ** newton_system);


ReNewtonSystem * ReNewtonSystemNew(int alphsize);
#endif
