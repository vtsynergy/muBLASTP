#ifndef _constants_
#define _constants_

#define constants_minSignedInt -2147483647
#define constants_maxSignedInt 2147483647
#define constants_max_int2 32767
#define constants_gappedExtensionDummyValue -666

#ifdef NCBI_BLAST
#define constants_sentinalScore -32768
#else
#define constants_sentinalScore -9999
#endif

#define constants_initialAllocUngappedExtensions 1000
#define constants_initialAllocAlignments 10000
#define constants_initialAllocGoodAlignments 10000
#define constants_initialAllocFinalAlignments 10000
#define constants_initialAllocCodewordQueryPositions 5
#define constants_initialTracebackAlloc 10000
#define constants_volumeMaxSize 2000000000
#define constants_maxNumVolumes 100
#define constants_maximumTracebackSize 10000000
#define constants_unpackRegionExtend 1000
#define constants_initialAllocUnpackRegions 10000
#define constants_databaseVersion 3
#define constants_initialSequenceData 100000

#define constants_initialSelectHits MAX_NUM_UNGAPPED_EXT

extern float Robinson_prob[];
extern float Nucleotide_prob[];

#endif
