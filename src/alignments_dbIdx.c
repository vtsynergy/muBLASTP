#include "blast.h"
#include <sys/time.h>

struct alignment *goodAlignQuery[BATCH_SIZE];
struct alignment *finalAlignQuery[BATCH_SIZE];
int4 numGoodAlignQuery[BATCH_SIZE];

int4 goodAlignCount_arr[MAX_NUM_THREADS] = {0};
int4 goodExtensionCount_arr[MAX_NUM_THREADS] = {0};
int4 gappedExtensionCount_arr[MAX_NUM_THREADS] = {0};
int4 traceCodeCount_arr[MAX_NUM_THREADS] = {0};
int4 extHitsCount_arr[BATCH_SIZE] = {0};
int4 numFinalAlignQuery[BATCH_SIZE] = {0};

struct alignment *goodAlignBuf_arr[MAX_NUM_THREADS];
struct ungappedExtension *goodExtensionBuf_arr[MAX_NUM_THREADS];
struct gappedExtension *gappedExtension_arr[MAX_NUM_THREADS];
unsigned char *traceCodeBuf_arr[MAX_NUM_THREADS];
HitPair *extHitsBuff_arr[BATCH_SIZE];

size_t goodAlignBufSize_arr[MAX_NUM_THREADS], 
       goodExtensionBufSize_arr[MAX_NUM_THREADS],
       gapExtensionBufSize_arr[MAX_NUM_THREADS],
       traceCodeBufSize_arr[MAX_NUM_THREADS],
       extHitsBuffSize_arr[BATCH_SIZE];

size_t finalAlignQueryBufSize[BATCH_SIZE];

void prelim_search_dbIdx(
        struct PSSMatrix *PSSMatrix_arr,
        struct scoreMatrix scoreMatrix, 
        int numQuery);

void prelim_search_dbIdx2(
        struct PSSMatrix *PSSMatrix_arr,
        struct scoreMatrix scoreMatrix, 
        int numQuery);

void traceback(
        struct PSSMatrix *PSSMatrix_arr,
        struct scoreMatrix scoreMatrix, 
        char *query_arr[],
        char *queryDescription_arr[],
        int numQuery);

void loadDbIdx_nextVolumn()
{

#ifndef COMPRESSED_INDEX
    free_dbindex();
    read_dbLookup(parameters_subjectDatabaseFile);
#else
    free_indexdb_cp();
    read_dbLookup_cp(parameters_subjectDatabaseFile);
#endif

}

void loadSubject(struct alignment *alignment)
{
    // Make copy of sequence
    unsigned char *subject = 
        (unsigned char *)global_malloc(sizeof(unsigned char) *
                alignment->encodedLength);
    subject++;
    memcpy(subject - 1, alignment->subject - 1, alignment->encodedLength);
    alignment->subject = subject;
}

void merge(int numQuery);

void alignments_dbIdx(
        struct PSSMatrix *PSSMatrix_arr,
        struct scoreMatrix scoreMatrix, 
        char *query_arr[],
        char *queryDescription_arr[],
        int numQuery)
{



    int numQueriesPerThread = ceil((float)numQuery / parameters_num_threads);


    int ii;
    for(ii = 0; ii < parameters_num_threads; ii++)
    {
        goodAlignBufSize_arr[ii] = MAX_ALIGNMENTS_PER_QUERY * numQueriesPerThread;
        goodExtensionBufSize_arr[ii] = MAX_EXTENSIONS_PER_QUERY * numQueriesPerThread;
        gapExtensionBufSize_arr[ii] = MAX_NUM_GAPEXT_PER_QUERY * numQueriesPerThread;
        traceCodeBufSize_arr[ii] = 1000 * parameters_numDisplayAlignments * numQueriesPerThread;

        goodAlignBuf_arr[ii] = 
            (struct alignment *)malloc(
                    sizeof(struct alignment) * goodAlignBufSize_arr[ii]);


        gappedExtension_arr[ii] = (struct gappedExtension *)malloc(
                sizeof(struct gappedExtension) * 
                gapExtensionBufSize_arr[ii]);


        goodExtensionBuf_arr[ii] = 
            (struct ungappedExtension *)malloc(
                    sizeof(struct ungappedExtension) * goodExtensionBufSize_arr[ii]);

        traceCodeBuf_arr[ii] = 
            (unsigned char *)malloc(traceCodeBufSize_arr[ii]);
    }

    for(ii = 0; ii < numQuery; ii++)
    {
        extHitsBuffSize_arr[ii] = MAX_EXTENSIONS_PER_QUERY;
        extHitsBuff_arr[ii] = (HitPair *)global_malloc(sizeof(HitPair) * 
                extHitsBuffSize_arr[ii]); 
    }

    for(ii = 0; ii < numQuery; ii++)
    {
        finalAlignQuery[ii] = (struct alignment *)global_malloc(
                sizeof(struct alignment) * parameters_numDisplayAlignments);
        finalAlignQueryBufSize[ii] = parameters_numDisplayAlignments;
    }

    int goodAlignOffset[MAX_NUM_THREADS] = {0};

    while(1)
    {
        prelim_search_dbIdx(PSSMatrix_arr, scoreMatrix, numQuery);

        merge(numQuery);

        if(readdb_volume + 1 < readdb_numberOfVolumes)
        {
            readdb_nextVolume_mem();
            loadDbIdx_nextVolumn();
        }
        else
        {
            break;
        }

    }


#if 1
    traceback(PSSMatrix_arr, scoreMatrix, query_arr, queryDescription_arr, numQuery);



    struct timeval start, end;
    gettimeofday(&start, NULL);

    for (ii = 0; ii < numQuery; ii++) {

        printf("Query= %s\n\n", queryDescription_arr[ii]);
        printf("Length=%ld\n", strlen(query_arr[ii]));

        struct finalAlignment *currentAlignment = 
            alignments_finalAlignments_multi[ii]->block;

        int jj;
        for(jj = 0; 
                jj < alignments_finalAlignments_multi[ii]->numEntries; 
                jj++, currentAlignment++)
        {

            int thread_id = currentAlignment->thread_id;
            struct alignment *alignment = currentAlignment->alignment;
            currentAlignment->description = alignment->description;
            ASSERT(alignment->gappedExtensionOffset != -1);
            alignment->gappedExtensions = 
                gappedExtension_arr[thread_id] 
                + alignment->gappedExtensionOffset;

            struct gappedExtension *currentExtension = alignment->gappedExtensions; 

            while(currentExtension != NULL)
            {
                if(currentExtension->nextOffset != -1)
                {
                    currentExtension->next = gappedExtension_arr[thread_id] 
                        + currentExtension->nextOffset;
                    currentExtension = currentExtension->next;
                }
                else
                {
                    currentExtension = NULL;
                }

            }

        }


        print_gappedAlignmentsBrief_multi(ii);
        print_gappedAlignmentsFull_multi(query_arr[ii], PSSMatrix_arr[ii], ii);

        alignments_free_multi2(ii);

        for(jj = 0; jj < numFinalAlignQuery[ii]; jj++)
        {
            struct alignment *alignment = &finalAlignQuery[ii][jj];
            free(alignment->subject - 1);
            free(alignment->description);
            free(alignment->ungappedExtensions);
        }

        free(extHitsBuff_arr[ii]);
        //free(goodAlignQuery[ii]);
        free(finalAlignQuery[ii]);
    }

    //int tid;
    //for(tid = 0; tid < parameters_num_threads; tid++)
    //{
    //for(ii = 0; ii < goodAlignCount_arr[tid]; ii++)
    //{
    //struct alignment *alignment = &goodAlignBuf_arr[tid][ii];
    //free(alignment->subject - 1);
    ////free(alignment->description);
    //}
    //}

    gettimeofday(&end, NULL);

    long print_time = ((end.tv_sec * 1000000 + end.tv_usec) -
            (start.tv_sec * 1000000 + start.tv_usec));

    fprintf(stderr, "print time: %f\n", (float)print_time * 1e-6);
#endif

    //int tid;
    //for(tid = 0; tid < parameters_num_threads; tid++)
    //{
    //for(ii = 0; ii < goodAlignCount_arr[tid]; ii++)
    //{
    //struct alignment alignment = goodAlignBuf_arr[tid][ii];
    //free(alignment.subject - 1);
    //}
    //}

    for(ii = 0; ii < parameters_num_threads; ii++)
    {
        free(goodExtensionBuf_arr[ii]);
        free(gappedExtension_arr[ii]);
        free(traceCodeBuf_arr[ii]);
        free(goodAlignBuf_arr[ii]);
        //free(subjectBuf_arr[ii]);
    }


}

void prelim_search_dbIdx2(
        struct PSSMatrix *PSSMatrix_arr,
        struct scoreMatrix scoreMatrix, 
        int numQuery)
{

    fprintf(stderr, "prelim search...");

    struct ungappedExtension **ungappedExtension_new_arr[parameters_num_threads];
    HitPair *selectHits1_arr[parameters_num_threads], 
            *selectHits2_arr[parameters_num_threads];

    BlastGapDP *dp_mem_arr = (BlastGapDP *)malloc(sizeof(BlastGapDP) * DP_MEM_SIZE * parameters_num_threads);
    BlastHSP *BlastHSP_arr = (BlastHSP *)malloc(sizeof(BlastHSP) * MAX_NUM_HSP * parameters_num_threads);

    size_t maxNumSecondBins = 0;

    int maxMaxNumSeq = 0, maxMaxDiag = 0;

    int ii, jj;
    for(jj = 0; jj < proteinLookup_numBlocks; jj++)
    {
#ifndef COMPRESSED_INDEX
        size_t numSeqBlk = proteinLookup_db_b[jj].numSeqBlk;
        size_t maxDiag = longestQueryLength + proteinLookup_db_b[jj].dbIdxblk_longestSeq;
#else
        size_t numSeqBlk = proteinLookup_db_blk_cp[jj].numSeqBlk;
        size_t maxDiag = longestQueryLength + proteinLookup_db_blk_cp[jj].dbIdxblk_longestSeq;
#endif
        size_t numSeqBins = (((numSeqBlk + 1) * maxDiag)) + 1; 
        if(numSeqBins > maxNumSecondBins)
        {
            maxNumSecondBins = numSeqBins;
            maxMaxNumSeq = numSeqBlk;
            maxMaxDiag = maxDiag;
        }
        //maxNumSecondBins = MAX(numSeqBins, maxNumSecondBins);
    }

    //fprintf(stderr, "maxNumSecondBins: %d\n", maxNumSecondBins);
    //

    for(ii = 0; ii < parameters_num_threads; ii++)
    {

        selectHits2_arr[ii] = (HitPair *)global_malloc(sizeof(HitPair) * 
                MAX_EXTENSIONS_PER_QUERY); 

        ungappedExtension_new_arr[ii] = 
            (struct ungappedExtension **)global_malloc(sizeof(struct ungappedExtension *) 
                    * MAX_EXTENSIONS_PER_QUERY);
    }

    uint2 *lastHits_arr = (uint2 *)global_malloc(sizeof(uint2) * maxNumSecondBins * parameters_num_threads);

    //fprintf(stderr, "maxNumSecondBins: %d maxMaxNumSeq: %d maxMaxDiag: %d lasthit_arr: %d selectHits_arr: %d\n",
    //maxNumSecondBins, maxMaxNumSeq, maxMaxDiag,
    //sizeof(uint2) * maxNumSecondBins * parameters_num_threads >> 20,
    //sizeof(HitPair) * maxNumSecondBins * 2 >> 20
    //);

    if(lastHits_arr == NULL)
    {
        fprintf(stderr, "failed to malloc lasthit_arr: %lu\n", 
                sizeof(uint2) * maxNumSecondBins * parameters_num_threads);
    }

    struct timeval start, end;
    gettimeofday(&start, NULL);


#pragma omp parallel num_threads(parameters_num_threads) default(shared) private(ii, jj)
    {
        int4 thread_id = omp_get_thread_num();

        uint2 *lastHits = lastHits_arr + maxNumSecondBins * thread_id;

        //int goodAlignCount = 0, goodExtensionCount = 0, gappedExtensionCount = 0;

        BlastGapDP *dp_mem = dp_mem_arr + DP_MEM_SIZE * thread_id;

        uint4 blast_numHits_t = 0;
        uint4 blast_numUngappedExtensions_t = 0;
        uint4 blast_numTriggerExtensions_t = 0;
        uint4 blast_numTriggerSequences_t = 0;

        BlastHSP *BlastHSP = BlastHSP_arr + MAX_NUM_HSP * thread_id;


        BlastIntervalTree *tree = Blast_IntervalTreeInit(0, longestQueryLength + 1,
                0, SEMIGAPPED_ROWSIZE + 1);
        BlastIntervalTree *private_tree = Blast_IntervalTreeInit(0, 
                longestQueryLength + 1, 0, SEMIGAPPED_ROWSIZE + 1);


        int bid;
        for (bid = 0; bid < proteinLookup_numBlocks; bid++) 
        {

            //#pragma omp single
            //read_dbIdxBlock(parameters_subjectDatabaseFile, bid);

            int kk;
#pragma omp for schedule(dynamic)
            for (kk = 0; kk < numQuery; kk++)
            {
                search_protein2hit_dbIdx_hitDetect(
                        thread_id, PSSMatrix_arr, 
                        kk, bid,
                        lastHits, 
                        &extHitsBuff_arr[kk],
                        &extHitsCount_arr[kk],
                        &extHitsBuffSize_arr[kk]);
            }


#if 1

#pragma omp for schedule(dynamic)
            for (kk = 0; kk < numQuery; kk++)
            {
                search_protein2hit_dbIdx_ungapExt(
                        thread_id, PSSMatrix_arr, scoreMatrix, 
                        kk, readdb_sequenceData, bid,
                        extHitsBuff_arr[kk],
                        extHitsCount_arr[kk],
                        selectHits2_arr[thread_id], 
                        &goodAlignBuf_arr[thread_id], 
                        &goodAlignCount_arr[thread_id], 
                        &goodAlignBufSize_arr[thread_id], 
                        &goodExtensionBuf_arr[thread_id], 
                        &goodExtensionCount_arr[thread_id], 
                        &goodExtensionBufSize_arr[thread_id], 
                        ungappedExtension_new_arr[thread_id], dp_mem, 
                        tree, private_tree, 
                        &blast_numHits_t, 
                        &blast_numUngappedExtensions_t, 
                        &blast_numTriggerExtensions_t, 
                        &blast_numTriggerSequences_t, 
                        BlastHSP);
            }

#endif
        }

        Blast_IntervalTreeFree2(private_tree);
        Blast_IntervalTreeFree2(tree);

    }

    gettimeofday(&end, NULL);

    long prelim_time = ((end.tv_sec * 1000000 + end.tv_usec) -
            (start.tv_sec * 1000000 + start.tv_usec));

    fprintf(stderr, "time: %f\n", (float)prelim_time * 1e-6);



    for(ii = 0; ii < parameters_num_threads; ii++)
    {
        //free(selectHits1_arr[ii]);
        free(selectHits2_arr[ii]);
        free(ungappedExtension_new_arr[ii]);
    }

    free(lastHits_arr);
    free(dp_mem_arr);
    free(BlastHSP_arr);
}


void prelim_search_dbIdx(
        struct PSSMatrix *PSSMatrix_arr,
        struct scoreMatrix scoreMatrix, 
        int numQuery)
{

    fprintf(stderr, "prelim search...");

    struct ungappedExtension **ungappedExtension_new_arr[parameters_num_threads];
    HitPair *selectHits1_arr[parameters_num_threads], 
            *selectHits2_arr[parameters_num_threads];

    BlastGapDP *dp_mem_arr = (BlastGapDP *)malloc(sizeof(BlastGapDP) * DP_MEM_SIZE * parameters_num_threads);
    BlastHSP *BlastHSP_arr = (BlastHSP *)malloc(sizeof(BlastHSP) * MAX_NUM_HSP * parameters_num_threads);

    size_t maxNumSecondBins = 0;

    int maxMaxNumSeq = 0, maxMaxDiag = 0;

    int ii, jj;
    for(jj = 0; jj < proteinLookup_numBlocks; jj++)
    {
#ifndef COMPRESSED_INDEX
        size_t numSeqBlk = proteinLookup_db_b[jj].numSeqBlk;
        size_t maxDiag = longestQueryLength + proteinLookup_db_b[jj].dbIdxblk_longestSeq;
#else
        size_t numSeqBlk = proteinLookup_db_blk_cp[jj].numSeqBlk;
        size_t maxDiag = longestQueryLength + proteinLookup_db_blk_cp[jj].dbIdxblk_longestSeq;
#endif
        size_t numSeqBins = (((numSeqBlk + 1) * maxDiag)) + 1; 
        if(numSeqBins > maxNumSecondBins)
        {
            maxNumSecondBins = numSeqBins;
            maxMaxNumSeq = numSeqBlk;
            maxMaxDiag = maxDiag;
        }
        //maxNumSecondBins = MAX(numSeqBins, maxNumSecondBins);
    }

    //fprintf(stderr, "maxNumSecondBins: %d\n", maxNumSecondBins);

    for(ii = 0; ii < parameters_num_threads; ii++)
    {
        selectHits1_arr[ii] = (HitPair *)global_malloc(sizeof(HitPair) * 
                maxNumSecondBins); 

        selectHits2_arr[ii] = (HitPair *)global_malloc(sizeof(HitPair) * 
                maxNumSecondBins); 

        if(selectHits1_arr == NULL || selectHits2_arr == NULL)
        {
            fprintf(stderr, "failed to malloc selectHits_arr: %lu\n", 
                    maxNumSecondBins);
        }

        ungappedExtension_new_arr[ii] = 
            (struct ungappedExtension **)global_malloc(sizeof(struct ungappedExtension *) 
                    * MAX_EXTENSIONS_PER_QUERY);
    }


    uint2 *lastHits_arr = (uint2 *)global_malloc(sizeof(uint2) * maxNumSecondBins * parameters_num_threads);

    //fprintf(stderr, "maxNumSecondBins: %d maxMaxNumSeq: %d maxMaxDiag: %d lasthit_arr: %d selectHits_arr: %d\n",
    //maxNumSecondBins, maxMaxNumSeq, maxMaxDiag,
    //sizeof(uint2) * maxNumSecondBins * parameters_num_threads >> 20,
    //sizeof(HitPair) * maxNumSecondBins * 2 >> 20
    //);

    if(lastHits_arr == NULL)
    {
        fprintf(stderr, "failed to malloc lasthit_arr: %lu\n", 
                sizeof(uint2) * maxNumSecondBins * parameters_num_threads);
    }

    struct timeval start, end;
    gettimeofday(&start, NULL);


#pragma omp parallel num_threads(parameters_num_threads) default(shared) private(ii, jj)
    {
        int4 thread_id = omp_get_thread_num();

        uint2 *lastHits = lastHits_arr + maxNumSecondBins * thread_id;

        //int goodAlignCount = 0, goodExtensionCount = 0, gappedExtensionCount = 0;

        BlastGapDP *dp_mem = dp_mem_arr + DP_MEM_SIZE * thread_id;

        BlastIntervalTree *tree = Blast_IntervalTreeInit(0, longestQueryLength + 1,
                0, SEMIGAPPED_ROWSIZE + 1);
        BlastIntervalTree *private_tree = Blast_IntervalTreeInit(0, 
                longestQueryLength + 1, 0, SEMIGAPPED_ROWSIZE + 1);

        uint4 blast_numHits_t = 0;
        uint4 blast_numUngappedExtensions_t = 0;
        uint4 blast_numTriggerExtensions_t = 0;
        uint4 blast_numTriggerSequences_t = 0;

        BlastHSP *BlastHSP = BlastHSP_arr + MAX_NUM_HSP * thread_id;

        int bid;
        for (bid = 0; bid < proteinLookup_numBlocks; bid++) 
        {

            //#pragma omp single
            //read_dbIdxBlock(parameters_subjectDatabaseFile, bid);

            int kk;
#pragma omp for schedule(dynamic) nowait
            for (kk = 0; kk < numQuery; kk++)
            {
#ifndef COMPRESSED_INDEX
                search_protein2hit_dbIdx_lasthit_radix(
#else
                        search_protein2hit_dbIdx_lasthit_radix_cp(
#endif
                            thread_id, PSSMatrix_arr, scoreMatrix, 
                            kk, readdb_sequenceData, bid,
                            lastHits, 
                            selectHits1_arr[thread_id], 
                            selectHits2_arr[thread_id], 
                            &goodAlignBuf_arr[thread_id], 
                            &goodAlignCount_arr[thread_id], 
                            &goodAlignBufSize_arr[thread_id], 
                            &goodExtensionBuf_arr[thread_id], 
                            &goodExtensionCount_arr[thread_id], 
                            &goodExtensionBufSize_arr[thread_id], 
                            ungappedExtension_new_arr[thread_id], dp_mem, 
                            tree, private_tree, 
                            &blast_numHits_t, 
                            &blast_numUngappedExtensions_t, 
                            &blast_numTriggerExtensions_t, 
                            &blast_numTriggerSequences_t, 
                            BlastHSP);
            }

            //#pragma omp single
            //free_dbIdxBlock(bid);

        }

        Blast_IntervalTreeFree2(private_tree);
        Blast_IntervalTreeFree2(tree);
    }

    gettimeofday(&end, NULL);

    long prelim_time = ((end.tv_sec * 1000000 + end.tv_usec) -
            (start.tv_sec * 1000000 + start.tv_usec));

    fprintf(stderr, "time: %f\n", (float)prelim_time * 1e-6);



    for(ii = 0; ii < parameters_num_threads; ii++)
    {
        free(selectHits1_arr[ii]);
        free(selectHits2_arr[ii]);
        free(ungappedExtension_new_arr[ii]);
    }

    free(lastHits_arr);
    free(dp_mem_arr);
    free(BlastHSP_arr);
}

void merge(int numQuery)
{
    struct timeval start, end;
    gettimeofday(&start, NULL);

    //memset(numGoodAlignQuery, 0, sizeof(int4) * numQuery);

    int tid, ii, jj;
    for(ii = 0; ii < numQuery; ii++)
    {
        numGoodAlignQuery[ii] = 0;

    }

    for(tid = 0; tid < parameters_num_threads; tid++)
    {
        for(ii = 0; ii < goodAlignCount_arr[tid]; ii++)
        {
            struct alignment alignment = goodAlignBuf_arr[tid][ii];
            numGoodAlignQuery[alignment.queryCount]++;
        }
    }

    for(ii = 0; ii < numQuery; ii++)
    {
        goodAlignQuery[ii] = (struct alignment *)global_malloc(
                sizeof(struct alignment) * numGoodAlignQuery[ii]);
        numGoodAlignQuery[ii] = 0;
    }

    //memset(numGoodAlignQuery, 0, sizeof(int4) * numQuery);

    for(tid = 0; tid < parameters_num_threads; tid++)
    {
        for(ii = 0; ii < goodAlignCount_arr[tid]; ii++)
        {
            struct alignment alignment = goodAlignBuf_arr[tid][ii];
            alignment.ungappedExtensions = 
                goodExtensionBuf_arr[tid] + alignment.ungappedExtensionOffset;
            goodAlignQuery[alignment.queryCount]
                [numGoodAlignQuery[alignment.queryCount]] = alignment;
            numGoodAlignQuery[alignment.queryCount]++;
        }

    }

    for(ii = 0; ii < numQuery; ii++)
    {

        if(numFinalAlignQuery[ii] + numGoodAlignQuery[ii] >= finalAlignQueryBufSize[ii])
        {
            finalAlignQueryBufSize[ii] = numFinalAlignQuery[ii] + numGoodAlignQuery[ii];
            finalAlignQuery[ii] = (struct alignment *)realloc(finalAlignQuery[ii], 
                    sizeof(struct alignment) * finalAlignQueryBufSize[ii]);
        }

        for(jj = 0; jj < numGoodAlignQuery[ii]; jj++)
        {
                finalAlignQuery[ii][numFinalAlignQuery[ii]] = goodAlignQuery[ii][jj];
                numFinalAlignQuery[ii]++;
        }

        for(jj = 0; jj < numFinalAlignQuery[ii]; jj++)
        {
            //if(jj < parameters_numDisplayAlignments)
            {
                struct alignment *alignment = &finalAlignQuery[ii][jj];
                if(alignment->inMemorySubject == 0)
                {
                    loadSubject(alignment);
                    alignment->description = 
                        descriptions_getDescription_mem(alignment->descriptionLocation, 
                                alignment->descriptionLength);
                    struct ungappedExtension *ungappedExtensions = 
                        (struct ungappedExtension*)global_malloc(sizeof(struct ungappedExtension) * 
                                alignment->numExtensions);
                    memcpy(ungappedExtensions, alignment->ungappedExtensions, 
                            sizeof(struct ungappedExtension) * alignment->numExtensions);
                    alignment->ungappedExtensions = ungappedExtensions; 

                    alignment->inMemorySubject = 1;
                }
            }
            //else
            //{
                //struct alignment *alignment = &finalAlignQuery[ii][jj];
                //if(alignment->inMemorySubject != 0)
                //{
                    //free(alignment->subject - 1);
                    //free(alignment->description);
                    //free(alignment->ungappedExtensions);
                //}
            //}
        }

        //if(numFinalAlignQuery[ii] > parameters_numDisplayAlignments)
        //numFinalAlignQuery[ii] = parameters_numDisplayAlignments;

        free(goodAlignQuery[ii]);
    }

    for(tid = 0; tid < parameters_num_threads; tid++)
    {
        goodAlignCount_arr[tid] = 0;
        goodExtensionCount_arr[tid] = 0;
    }

    gettimeofday(&end, NULL);
    long merge_time = ((end.tv_sec * 1000000 + end.tv_usec) -
            (start.tv_sec * 1000000 + start.tv_usec));

    fprintf(stderr, "merge time: %f\n", (float)merge_time * 1e-6);
}

void traceback(
        struct PSSMatrix *PSSMatrix_arr,
        struct scoreMatrix scoreMatrix, 
        char *query_arr[],
        char *queryDescription_arr[],
        int numQuery)
{

    Uint8 **query_words_arr = 
        (Uint8 **)malloc(sizeof(Uint8 *) * parameters_num_threads);

    Blast_AminoAcidComposition *query_composition_arr = 
        (Blast_AminoAcidComposition *)malloc(sizeof(Blast_AminoAcidComposition) * numQuery);

    Blast_CompositionWorkspace **NRrecord_arr = 
        (Blast_CompositionWorkspace **)malloc(sizeof(Blast_CompositionWorkspace *) * parameters_num_threads); ;

    Blast_MatrixInfo *matrixInfo_arr = 
        (Blast_MatrixInfo *)malloc(sizeof(Blast_MatrixInfo) * parameters_num_threads);

    BlastGapAlignStruct **gap_align_arr = 
        (BlastGapAlignStruct **)malloc(sizeof(BlastGapAlignStruct *) * parameters_num_threads);

    void **matrix_arr = (void **)malloc(sizeof(void *) * parameters_num_threads);

    unsigned char **seq_data_arr = 
        (unsigned char **)malloc(sizeof(unsigned char *) * parameters_num_threads);

    ReNewtonSystem **newton_system_arr = 
        (ReNewtonSystem **)malloc(sizeof(ReNewtonSystem *) * parameters_num_threads);

    double ***grads_arr = (double ***)malloc(sizeof(double **) * parameters_num_threads);
    double ***Scores_arr = (double ***)malloc(sizeof(double **) * parameters_num_threads);

    //ReNewtonSystem *newton_system;
    double **z_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    double **resids_x_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    double **resids_z_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    double **old_scores_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    double **workspace_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);

    char **querySeq_arr = (char **)malloc(sizeof(char *) * numQuery);

    int ii, jj;
    for(ii = 0; ii < numQuery; ii++)
    {
        querySeq_arr[ii] = (char *)malloc(sizeof(char) * PSSMatrix_arr[ii].length);
        for(jj = 0; jj < PSSMatrix_arr[ii].length; jj++)
        {
            ASSERT(PSSMatrix_arr[ii].queryCodes[jj] < FSA_AA_SIZE);
            querySeq_arr[ii][jj] = to_ncbi[PSSMatrix_arr[ii].queryCodes[jj]];
        }

        Blast_ReadAaComposition_fsa(
                &query_composition_arr[ii], 
                BLASTAA_SIZE, 
                PSSMatrix_arr[ii].queryCodes, 
                PSSMatrix_arr[ii].length);

    }

    int alphsize = COMPO_NUM_TRUE_AA;
    int n  = alphsize * alphsize;
    int mA = 2 * alphsize - 1;
    int m  = mA + 1;

    for(ii = 0; ii < parameters_num_threads; ii++)
    {
        NRrecord_arr[ii] = Blast_CompositionWorkspaceNew();
        Blast_MatrixInfoNew2(&matrixInfo_arr[ii], BLASTAA_SIZE, BLASTAA_SIZE, FALSE);
        BLAST_GapAlignStructNew(readdb_longestSequenceLength, &gap_align_arr[ii]);
        matrix_arr[ii] = (void *)_PSIAllocateMatrix(BLASTAA_SIZE, BLASTAA_SIZE, sizeof(Int4)); 
        seq_data_arr[ii] = (unsigned char *)malloc((readdb_longestSequenceLength + 10) * sizeof(unsigned char));
        newton_system_arr[ii] = ReNewtonSystemNew(alphsize); 
        grads_arr[ii] = Nlm_DenseMatrixNew(2, n);
        Scores_arr[ii] = Nlm_DenseMatrixNew(BLASTAA_SIZE, BLASTAA_SIZE);
        resids_x_arr[ii] = (double *) malloc(n * sizeof(double));
        resids_z_arr[ii] = (double *) malloc((mA + 1) * sizeof(double));
        z_arr[ii] = (double *) calloc( mA + 1,   sizeof(double));
        old_scores_arr[ii] = (double *) malloc(n * sizeof(double));
        workspace_arr[ii] = (double *) malloc(n * sizeof(double));

        edit_script_num_rows_t[ii] = EDIT_SCRIPT_MAX_NUM_ROWS;
        edit_script_t[ii] = (Uint1**) malloc(sizeof(Uint1*) * edit_script_num_rows_t[ii]);
        edit_start_offset_t[ii] = 
            (Int4*) malloc(sizeof(Int4) * edit_script_num_rows_t[ii]);

        query_words_arr[ii] = (Uint8 *)malloc(sizeof(Uint8) * longestQueryLength);
    }


    fprintf(stderr, "traceback search...");

    struct timeval start, end;
    gettimeofday(&start, NULL);

#pragma omp parallel num_threads(parameters_num_threads) default(shared) private(ii, jj)
    {

        struct timeval start_t, end_t;
        gettimeofday(&start_t, NULL);

        int4 thread_id = omp_get_thread_num();

        int qid;
#pragma omp for schedule(dynamic) nowait
        for (qid = 0; qid < numQuery; qid++) {
            alignments_getTracebacks_ncbi(
                    thread_id,
                    finalAlignQuery[qid],
                    numFinalAlignQuery[qid],     
                    PSSMatrix_arr[qid], 
                    scoreMatrix, &query_composition_arr[qid], 
                    qid, NRrecord_arr[thread_id], 
                    querySeq_arr[qid], 
                    &matrixInfo_arr[thread_id], 
                    gap_align_arr[thread_id], 
                    matrix_arr[thread_id], 
                    seq_data_arr[thread_id], 
                    newton_system_arr[thread_id], 
                    grads_arr[thread_id], 
                    Scores_arr[thread_id], 
                    z_arr[thread_id], 
                    resids_x_arr[thread_id], 
                    resids_z_arr[thread_id], 
                    old_scores_arr[thread_id], 
                    workspace_arr[thread_id], 
                    &gappedExtension_arr[thread_id], 
                    &gappedExtensionCount_arr[thread_id], 
                    &gapExtensionBufSize_arr[thread_id],
                    &traceCodeBuf_arr[thread_id],
                    &traceCodeCount_arr[thread_id],
                    &traceCodeBufSize_arr[thread_id],
                    query_words_arr[thread_id]);
        }

        gettimeofday(&end_t, NULL);
        long traceback_time_t = ((end_t.tv_sec * 1000000 + end_t.tv_usec) -
                (start_t.tv_sec * 1000000 + start_t.tv_usec));

        //fprintf(stderr, "tid: %d traceback time: %f\n", 
        //thread_id,
        //(float)traceback_time_t * 1e-6);

    }

    gettimeofday(&end, NULL);
    long traceback_time = ((end.tv_sec * 1000000 + end.tv_usec) -
            (start.tv_sec * 1000000 + start.tv_usec));

    fprintf(stderr, "time: %f\n", (float)traceback_time * 1e-6);

    for(ii = 0; ii < numQuery; ii++)
    {
        free(querySeq_arr[ii]);
    }

    for(ii = 0; ii < parameters_num_threads; ii++)
    {
        free(workspace_arr[ii]);
        free(old_scores_arr[ii]);
        free(z_arr[ii]);
        free(resids_z_arr[ii]);
        free(resids_x_arr[ii]);
        Nlm_DenseMatrixFree(&grads_arr[ii]);
        Nlm_DenseMatrixFree(&Scores_arr[ii]);
        ReNewtonSystemFree(&newton_system_arr[ii]);
        free(seq_data_arr[ii]);
        _PSIDeallocateMatrix((void **)matrix_arr[ii], BLASTAA_SIZE);
        BLAST_GapAlignStructFree(gap_align_arr[ii]);
        Blast_CompositionWorkspaceFree(&NRrecord_arr[ii]);
        Blast_MatrixInfoFree2(&matrixInfo_arr[ii]);
        free(edit_start_offset_t[ii]);
        free(edit_script_t[ii]);
        free(query_words_arr[ii]);
    }

    free(query_words_arr);
    free(gap_align_arr);
    free(query_composition_arr);
    free(NRrecord_arr);
    free(matrix_arr);
    free(matrixInfo_arr);
    free(seq_data_arr);
    free(newton_system_arr);
    free(grads_arr);
    free(Scores_arr);
    free(z_arr);
    free(resids_z_arr);
    free(resids_x_arr);
    free(old_scores_arr);
    free(workspace_arr);
    free(querySeq_arr);
}
