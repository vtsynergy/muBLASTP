#ifndef _blast_
#define _blast_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <errno.h>
#include <ctype.h>

#define NEIGHBOR_INDEX
#define NCBI_BLAST
#define MAX_NUM_THREADS 56 

#define BATCH_SIZE 1000

#define MAX_NUM_UNGAPPED_EXT 40000000 
#define MAX_ALIGNMENTS_PER_QUERY 2000
#define MAX_EXTENSIONS_PER_QUERY 64000
#define MAX_NUM_GAPEXT_PER_QUERY 5000

#define SEMIGAPPED_ROWSIZE 35213
#define INIT_TREE_SIZE 8000
#define DP_MEM_SIZE 8000
#define MAX_NUM_HSP 8000
#define EDIT_SCRIPT_MAX_NUM_ROWS 200


#define int4 int32_t
#define uint4 uint32_t
#define int8 int64_t
#define uint8 uint64_t
#define int2 int16_t
#define uint2 uint16_t
#define hit_t uint32_t
#define byte unsigned char

#define INT2_MAX 32767
#define INT2_MIN (-32768)

#define INT4_MAX 2147483647
#define INT4_MIN (-2147483648)

#include "global.h"
#include "constants.h"
#include "vbyte.h"
#include "memSingleBlock.h"
#include "memBlocks.h"
#include "readFile.h"
#include "encoding.h"
#include "writedb.h"
#include "readNcbidb.h"
#include "readdb.h"
#include "parameters.h"
#include "scoreMatrix.h"
#include "PSSMatrix.h"
#include "karlin.h"
#include "statistics.h"
#include "qPosList.h"
#include "wordLookupDFA.h"
#include "nucleotideLookup.h"
#include "readFasta.h"
#include "hitMatrix.h"
#include "ungappedExtension.h"
#include "gappedExtension.h"
#include "fasterGappedExtension.h"
#include "descriptions.h"
#include "bytepackGappedScoring.h"
#include "fasterBytepackGappedScoring.h"
#include "tableGappedScoring.h"
#include "oldSemiGappedScoring.h"
#include "oldGappedScoring.h"
#include "semiGappedScoring.h"
#include "gappedScoring.h"
#include "nuGappedScoring.h"
#include "smithWatermanScoring.h"
#include "smithWatermanTraceback.h"
#include "alignments.h"
#include "unpack.h"
#include "print.h"
#include "index.h"
#include "postings.h"
#include "hashcounter.h"
#include "search.h"
#include "wildcards.h"
#include "dust.h"
#include "seg.h"


#include "utils.h"
#include "dbLookup.h"
#include "dbLookup_compress.h"
#include "queryLookup.h"
#include "ungappedExtension_multi.h"
#include "wordLookupDFA_multi.h"
#include "hitMatrix_multi.h"
#include "fasterGappedExtension_multi.h"
#include "semiGappedScoring_multi.h"
#include "gappedScoring_multi.h"
#include "unpack_multi.h"
#include "ncbi_gapalign.h"
#include "ncbi_stat.h"
#include "ncbi_tree.h"
#include "ncbi_traceback.h"
#include "search_dbIdx.h"
#include "ncbi_matrix.h"
#include "ncbi_composition_adjustment.h"
#include "ncbi_filter.h"
#include "search_queryIdx.h"
#include "search_dbIdx.h"
#include "alignments_dbIdx.h"
#include "alignments_queryIdx.h"
#include "alignments_multi.h"
#include "gappedExtension_multi.h"

#define ASSERT(x) if(!(x)) {fprintf(stderr, "ERR in line %d at file: %s\n", __LINE__, __FILE__); exit(1);}

#define MAX(a,b) ((a)>=(b)?(a):(b))
#define MIN(a,b) ((a)>(b)?(b):(a))


#define TRUE 1
#define FALSE 0

#endif
