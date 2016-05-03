#ifndef _global_
#define _global_

#define maximum(a,b) ((a > b) ? a : b)
#define minimum(a,b) ((a < b) ? a : b)

// BLAST datatypes

// This structure contains values for a scoring matrix such as BLOSUM62
struct scoreMatrix
{
	int2** matrix;
	int4 highestValue;
	int4 lowestValue;
    float averageMatchScore;
};

// This structure contains a PSSM (Position Specific Scoring Matrix) which has
// a column for each amino acid in original query sequence (of length "length")
// and 25 rows.
struct PSSMatrix
{
    int4 length;
    int4 strandLength;
	int4 highestValue;
	int4 lowestValue;
	int2** matrix;
	unsigned char* queryCodes;
    unsigned char* bestMatchCodes;
    unsigned char* bytePackedCodes;
	unsigned char* xorCodes;
};

// A query/subject coordinate pair
struct coordinate
{
    int4 queryOffset;
	int4 subjectOffset;
};


// An alignment trace resulting from dynamic programming
struct trace
{
	uint4 length;
	uint4 queryStart;
	uint4 subjectStart;
	unsigned char* traceCodes;
    uint4 traceCodeOff;
};



// Information about an ungapped extension
struct ungappedExtension
{
	struct trace trace;
    uint4 sequenceCount;
    uint4 queryCount;
	struct coordinate start;
	struct coordinate end;
    struct coordinate gap_start;
	struct coordinate seed;
	int4 nominalScore;
    double eValue;
	char status;
	struct ungappedExtension* next;
};

// A region either copied or unpacked by blast
struct unpackRegion
{
	int4 startOffset;
    int4 endOffset;
	unsigned char* unpackedSubject;
    unsigned char* subject;
    int4 subjectLength;
};

// Result of dynamic programming
struct dpResults
{
	struct coordinate best;
	uint4 bestScore;
	unsigned char** traceback;
};

// A gapped alignment
struct gappedExtension
{
	struct trace trace;
	int4 nominalScore;
	int4 queryEnd;
    int4 subjectStart;
	int4 subjectEnd;
	float normalizedScore;
	double eValue;
	struct gappedExtension* next;
	int nextOffset;
};

// Information about the alignment between the query and a subject
struct alignment
{
    int4 descriptionLocation;
    int4 descriptionLength;
    char* description;
    unsigned char* subject;
    int4 subjectLength;
    struct ungappedExtension* ungappedExtensions;
    int4  ungappedExtensionOffset;
    struct gappedExtension* gappedExtensions;
    int4  gappedExtensionOffset;
    unsigned char* edits;
    int4 encodedLength;
    char joinChecked;
    char inMemorySubject;
    struct unpackRegion* unpackRegions;
    uint4 numUnpackRegions;
    uint4 cluster;
    uint4 sequenceCount;
    uint4 queryCount;
    uint4 numExtensions;
    int volumnNumber;
    double best_eValue;
};

/*typedef struct alignment align_t;*/

/*typedef struct*/
/*{*/
/*int4 descriptionLocation;*/
/*int4 descriptionLength;*/
/*unsigned char* subject;*/
/*int4 subjectLength;*/
/*struct ungappedExtension* ungappedExtensions;*/
/*int4  ungappedExtensionOffset;*/
/*struct gappedExtension* gappedExtensions;*/
/**//*unsigned char* edits;*/
/*int4 encodedLength;*/
/**//*char joinChecked;*/
/*char inMemorySubject;*/
/**//*struct unpackRegion* unpackRegions;*/
/**//*uint4 numUnpackRegions;*/
/**//*uint4 cluster;*/
/*uint4 sequenceCount;*/
/*uint4 queryCount;*/
/*uint4 numExtensions;*/
/*double best_eValue;*/
/*} align_t;*/

// A final alignment above cutoff
struct finalAlignment
{
	int4 highestNominalScore;
    char* description;
	struct alignment* alignment;
    int4 thread_id;
};

extern int rank, num_procs;
extern uint8 total_numberOfLetters; 

// Timing variables
extern int4 blast_prepTime, blast_searchTime;
extern int4 blast_gappedScoreTime, blast_gappedExtendTime, blast_finalizeTime;
extern int4 blast_semiGappedScoreTime, blast_copyTime, blast_unpackTime;

// BLAST statistics
extern uint4 blast_numHits;
extern uint4 blast_numUngappedExtensions, blast_numTriggerExtensions, blast_numTriggerSequences;
extern uint4 blast_numGapped;
extern uint4 blast_numSemiGapped;
extern uint4 blast_numExtensionsPruned;
extern uint4 blast_numAttemptedJoin, blast_numSuccessfullyJoined;
extern uint4 blast_numGoodAlignments;
extern uint4 blast_numGoodExtensions;
extern uint4 blast_totalUnpacked;
extern uint4 blast_totalCopied;
extern uint4 blast_numExpandedSequences;

//DB_INDEX
extern uint4 longestQueryLength;

// BLAST global variables
extern int4 blast_ungappedNominalTrigger;
extern int4 blast_gappedNominalCutoff;
extern int4 blast_gappedNominalCutoff, blast_nominalR1cutoff, blast_nominalR2cutoff;
extern int4 blast_dynamicGappedNominalCutoff, blast_dynamicNominalR1cutoff;
extern int4 blast_dloc;
extern char* blast_queryDescription;
extern char* blast_queryDescription_multi[];

extern uint4 blast_numQuery;
extern uint4 blast_numBlocks;

// Initialize global variables
void global_initialize();
// Convert a 32-bit int4eger int4o a string with commas
char* global_int4toString(uint4 number);
// Convert a 64-bit int4eger int4o a string with commas
char* global_int8toString(uint8 number);

extern uint8 global_totalMalloc;

// Malloc new memory, check that malloc was successful
void* global_malloc(size_t size);
// Realloc memory, check that realloc was successful
void* global_realloc(void* ptr, size_t size);

// Free the global convert string
void global_free();

//DB_INDEX
extern int proteinLookup_numBlock;
extern uint4 dbIdx_block_size;

// Initialize global variables
void global_initialize_multi();
// BLAST statistics
extern uint4 blast_numHits_multi[];
extern uint4 blast_numUngappedExtensions_multi[], blast_numTriggerExtensions_multi[], blast_numTriggerSequences_multi[];
extern uint4 blast_numGapped_multi[];
extern uint4 blast_numGappedExtension_multi[];
extern uint4 blast_numSemiGapped_multi[];
extern uint4 blast_numExtensionsPruned_multi[];
extern uint4 blast_numAttemptedJoin_multi[], blast_numSuccessfullyJoined_multi[];
extern uint4 blast_numGoodAlignments_multi[];
extern uint4 blast_numGoodExtensions_multi[];
extern uint4 blast_totalUnpacked_multi[];
extern uint4 blast_totalCopied_multi[];
extern uint4 blast_numExpandedSequences_multi[];
extern uint2 blast_tgSize;
extern int4 blast_dynamicGappedNominalCutoff_multi[], blast_dynamicNominalR1cutoff_multi[];
extern int4 blast_dloc_multi[];

extern int4 blast_ungappedNominalTrigger_multi[];
extern int4 blast_gappedNominalCutoff_multi[], blast_nominalR1cutoff_multi[], blast_nominalR2cutoff_multi[];

/*extern long long blast_hitDetectCycle, blast_pesudoCycle,  blast_ungappedExtCycle, blast_sortCycle, blast_findGoodCycle, blast_findFinalCycle, blast_getFinalCycle, blast_getTrackCycle; */


#endif

