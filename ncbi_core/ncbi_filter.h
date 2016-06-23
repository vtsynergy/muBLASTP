/** This bit is on if the query is protein */
#define PROTEIN_QUERY_MASK      (0x1<<0)
/** This bit is on if the subject is protein */
#define PROTEIN_SUBJECT_MASK    (0x1<<1)
/** This bit is on if the query is nucleotide */
#define NUCLEOTIDE_QUERY_MASK   (0x1<<2)
/** This bit is on if the subject is nucleotide */
#define NUCLEOTIDE_SUBJECT_MASK (0x1<<3)
/** This bit is on if the query is translated */
#define TRANSLATED_QUERY_MASK   (0x1<<4)
/** This bit is on if the subject is translated */
#define TRANSLATED_SUBJECT_MASK (0x1<<5)
/** This bit is on if the query is a PSSM (PSI-BLAST) */
#define PSSM_QUERY_MASK         (0x1<<6)
/** This bit is on if the subject is a PSSM (RPS-BLAST) */
#define PSSM_SUBJECT_MASK       (0x1<<7)
/** This bit is on if the query includes a pattern (PHI-BLAST) */
#define PATTERN_QUERY_MASK      (0x1<<8)

typedef enum {
    eBlastTypeBlastp        = (PROTEIN_QUERY_MASK | PROTEIN_SUBJECT_MASK),
    eBlastTypeBlastn        = (NUCLEOTIDE_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK),
    eBlastTypeBlastx        = (NUCLEOTIDE_QUERY_MASK | PROTEIN_SUBJECT_MASK | 
                               TRANSLATED_QUERY_MASK),
    eBlastTypeTblastn       = (PROTEIN_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK | 
                               TRANSLATED_SUBJECT_MASK),
    eBlastTypeTblastx       = (NUCLEOTIDE_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK
                               | TRANSLATED_QUERY_MASK 
                               | TRANSLATED_SUBJECT_MASK),
    eBlastTypePsiBlast      = (PSSM_QUERY_MASK | eBlastTypeBlastp),
    eBlastTypePsiTblastn    = (PSSM_QUERY_MASK | eBlastTypeTblastn),
    eBlastTypeRpsBlast      = (PSSM_SUBJECT_MASK | eBlastTypeBlastp),
    eBlastTypeRpsTblastn    = (PSSM_SUBJECT_MASK | eBlastTypeBlastx),
    eBlastTypePhiBlastp     = (PATTERN_QUERY_MASK | eBlastTypeBlastp),
    eBlastTypePhiBlastn     = (PATTERN_QUERY_MASK | eBlastTypeBlastn),
    eBlastTypeUndefined     = 0x0
} EBlastProgramType;


typedef struct BlastCompo_SequenceData {
    Uint1 * data;                /**< amino acid or nucleotide data */
    int length;                  /**< the length of data. For amino acid data
                                   &data[-1] is a valid address and
                                   data[-1] == 0. */
    Uint1 * buffer;               /**< if non-nil, points to memory that
                                    must be freed when this instance of
                                    BlastCompo_SequenceData is deleted. */
} BlastCompo_SequenceData;

typedef struct SDustOptions {
    int level;
    int window;
    int linker;  /**< min distance to link segments. */
} SDustOptions;

  
typedef struct SSegOptions {
    int window;     /**< initial window to trigger further work. */
    double locut;
    double hicut;
} SSegOptions;

typedef struct SRepeatFilterOptions {
    char* database;   /**< Nucleotide database for mini BLAST search. */
} SRepeatFilterOptions;

typedef struct SWindowMaskerOptions {
    int          taxid;    /**< Select masking database for this TaxID. */
    const char * database; /**< Use winmasker database at this location. */
} SWindowMaskerOptions;

typedef struct SBlastFilterOptions {
    Boolean mask_at_hash;         /**< mask query only for lookup table creation */
    SDustOptions* dustOptions;    /**< low-complexity filtering for nucleotides. */
    SSegOptions* segOptions;      /**< low-complexity filtering for proteins sequences 
            (includes translated nucleotides). */
    SRepeatFilterOptions* repeatFilterOptions;  /**< for organism specific repeat filtering. */
    SWindowMaskerOptions* windowMaskerOptions;  /**< organism specific filtering with window masker. */
} SBlastFilterOptions;


int
s_DoSegSequenceData(BlastCompo_SequenceData * seqData,
        EBlastProgramType program_name,
        Boolean* is_seq_biased);


Boolean s_preliminaryTestNearIdentical(Int4 queryLength,
        struct ungappedExtension *ungappedExtension,
        double cutoff);


Boolean
s_IsContained(struct ungappedExtension *in_align,
        struct alignment *alignments,
        double lambda);


Boolean
s_WithDistinctEnds(struct gappedExtension *newAlign,
                   struct gappedExtension *align,
                   Boolean is_same_adjustment);
