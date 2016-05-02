#include "blast.h"

void alignments_queryIdx(
        struct PSSMatrix *PSSMatrix_arr, 
        struct scoreMatrix scoreMatrix,
        char **query_arr,
        char **queryDescription_arr,
        int numQuery
        ) 
{

    struct timeval start, end;
    gettimeofday(&start, NULL);

    int numQueriesPerThread = ceil((float)numQuery / parameters_num_threads);
    size_t hitBufSize = MAX_NUM_HIT_PER_QUERY * numQueriesPerThread; 

    hit_t *secondBin_arr = (hit_t *)malloc(sizeof(hit_t) * 
            hitBufSize * parameters_num_threads);
    //fprintf(stderr, "secondBin_arr: %p size: %lu\n", secondBin_arr, 
    //hitBufSize * parameters_num_threads * sizeof(hit_t) >> 20);
    uint4 *numExtHit_arr = (uint4 *)malloc(sizeof(uint4) * numQuery *
            parameters_num_threads);

    int maxDiag = longestQueryLength + readdb_longestSequenceLength;
    int maxNumSecondBins = (((numQuery + 1) * maxDiag)) + 1;

    uint2 *lastHits_arr = (uint2 *)malloc(sizeof(uint2) * maxNumSecondBins * parameters_num_threads);

    size_t goodExtensionBufSize = MAX_EXTENSIONS_PER_QUERY * numQueriesPerThread;
    struct ungappedExtension *goodExtensionBuf_arr = (struct ungappedExtension *)
        malloc(sizeof(struct ungappedExtension) * 
                goodExtensionBufSize * parameters_num_threads);

    size_t goodAlignBufSize = MAX_ALIGNMENTS_PER_QUERY * numQueriesPerThread;
    struct alignment *goodAlignBuf_arr = (struct alignment *)malloc(sizeof(struct alignment) * goodAlignBufSize * parameters_num_threads);

    struct ungappedExtension **ungappedExtension_new_arr = 
        (struct ungappedExtension **)malloc(sizeof(struct ungappedExtension *) 
                * MAX_NUM_TRIGGER_EXT * parameters_num_threads);

    BlastGapDP *dp_mem_arr = (BlastGapDP *)malloc(sizeof(BlastGapDP) * 
            DP_MEM_SIZE * parameters_num_threads);

    BlastHSP *BlastHSPs_arr = (BlastHSP *)malloc(sizeof(BlastHSP) * MAX_NUM_HSP *
            parameters_num_threads);

    uint4 *numHitsPerQPos_arr = (uint4 *)malloc(sizeof(uint4) * 
            longestQueryLength * parameters_num_threads);

    Blast_MatrixInfo *matrixInfo_arr = (Blast_MatrixInfo *)malloc(sizeof(Blast_MatrixInfo) * parameters_num_threads);

    Blast_CompositionWorkspace **NRrecord = (Blast_CompositionWorkspace **)malloc(sizeof(Blast_CompositionWorkspace *) * parameters_num_threads); ;

    BlastGapAlignStruct **gap_align_arr = (BlastGapAlignStruct **)malloc(sizeof(BlastGapAlignStruct *) * parameters_num_threads);

    void **matrix_arr = (void **)malloc(sizeof(void *) * parameters_num_threads);

    unsigned char **seq_data_arr = (unsigned char **)malloc(sizeof(unsigned char *) * parameters_num_threads);

    ReNewtonSystem **newton_system_arr = (ReNewtonSystem **)malloc(sizeof(ReNewtonSystem *) * parameters_num_threads);

    double ***grads_arr = (double ***)malloc(sizeof(double **) * parameters_num_threads);
    double ***Scores_arr = (double ***)malloc(sizeof(double **) * parameters_num_threads);

    //ReNewtonSystem *newton_system;
    double **z_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    double **resids_x_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    double **resids_z_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    double **old_scores_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    double **workspace_arr = (double **)malloc(sizeof(double *) * parameters_num_threads);
    //double **grads;
    //double **Scores;


    int alphsize = COMPO_NUM_TRUE_AA;
    int n  = alphsize * alphsize;
    int mA = 2 * alphsize - 1;
    int m  = mA + 1;

    int ii, jj;
    for(ii = 0; ii < parameters_num_threads; ii++)
    {
        NRrecord[ii] = Blast_CompositionWorkspaceNew();
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
    }

    char **querySeq_arr = (char **)malloc(sizeof(char *) * numQuery);

    Blast_AminoAcidComposition *query_composition_arr = (Blast_AminoAcidComposition *)malloc(sizeof(Blast_AminoAcidComposition) * numQuery);

    for(ii = 0; ii < numQuery; ii++)
    {
        querySeq_arr[ii] = (char *)malloc(sizeof(char) * PSSMatrix_arr[ii].length);
        for(jj = 0; jj < PSSMatrix_arr[ii].length; jj++)
        {
            ASSERT(PSSMatrix_arr[ii].queryCodes[jj] < FSA_AA_SIZE);
            querySeq_arr[ii][jj] = to_ncbi[PSSMatrix_arr[ii].queryCodes[jj]];
        }
        Blast_ReadAaComposition_fsa(&query_composition_arr[ii], BLASTAA_SIZE, PSSMatrix_arr[ii].queryCodes, PSSMatrix_arr[ii].length);
    }

    struct gappedExtension *gappedExtension_arr = (struct gappedExtension *)malloc(sizeof(struct gappedExtension) * 
            MAX_NUM_GAPEXT_PER_QUERY * numQueriesPerThread * parameters_num_threads);

    gettimeofday(&end, NULL);
    long prepare_time = ((end.tv_sec * 1000000 + end.tv_usec) -
            (start.tv_sec * 1000000 + start.tv_usec));
    fprintf(stderr, "Prepare time: %f\n", (float)prepare_time * 1e-6);


#pragma omp parallel num_threads(parameters_num_threads) default(shared) private(start, end)
        //secondBin_arr, numExtHit_arr, lastHits_arr,\
        //goodAlignBuf_arr, goodExtensionBuf_arr, ungappedExtension_new_arr,\
        //dp_mem_arr, BlastHSPs_arr, numQuery, \
        //goodAlignBufSize, goodExtensionBufSize, scoreMatrix, \
        //readdb_numberOfSequences, readdb_sequenceData, PSSMatrix_arr,\
        //numHitsPerQPos_arr, longestQueryLength, \
        //gappedExtension_arr, numQueriesPerThread, workspace_arr,\
        //old_scores_arr, resids_z_arr, resids_x_arr, z_arr,\
        //Scores_arr, grads_arr, newton_system_arr, \
        //seq_data_arr, matrix_arr, gap_align_arr, \
        //matrixInfo_arr, querySeq_arr, NRrecord, \
        //query_composition_arr, \
        //stderr)
    {
        int goodExtensionCount = 0;
        int goodAlignCount = 0;
        int gappedExtension_cnt = 0;

        int4 thread_id = omp_get_thread_num();
        hit_t *secondBin = secondBin_arr + hitBufSize * thread_id;

        uint4 *numExtHit = numExtHit_arr + numQuery * thread_id;
        uint2 *lastHits = lastHits_arr + maxNumSecondBins * thread_id;

        struct alignment *goodAlignBuf = goodAlignBuf_arr + goodAlignBufSize * thread_id;

        struct ungappedExtension *goodExtensionBuf = goodExtensionBuf_arr +
            goodExtensionBufSize * thread_id;

        struct ungappedExtension **ungappedExtension_new = ungappedExtension_new_arr +
            MAX_NUM_TRIGGER_EXT * thread_id;

        BlastGapDP *dp_mem = dp_mem_arr + DP_MEM_SIZE * thread_id;

        BlastIntervalTree *tree = Blast_IntervalTreeInit(0, longestQueryLength + 1,
                0, SEMIGAPPED_ROWSIZE + 1);

        BlastIntervalTree *private_tree = Blast_IntervalTreeInit(0, 
                longestQueryLength + 1, 0, SEMIGAPPED_ROWSIZE + 1);

        BlastHSP *BlastHSPs = BlastHSPs_arr + MAX_NUM_HSP * thread_id; 

        uint4 *numHitsPerQPos = numHitsPerQPos_arr + longestQueryLength * thread_id;

        gettimeofday(&start, NULL);

        int ii, jj;
        for(jj = 0; jj < numQueryBlk + 1; jj++)
        {
            int kk;
#pragma omp for schedule(dynamic)
            for (kk = 0; kk < readdb_numberOfSequences; kk++) {
                search_protein2hit_queryIdx(
                        thread_id,
                        kk,
                        jj,
                        PSSMatrix_arr, 
                        readdb_sequenceData,
                        readdb_numberOfSequences, 
                        numQuery,
                        secondBin, 
                        numExtHit, 
                        lastHits, 
                        scoreMatrix, 
                        goodExtensionBuf, 
                        &goodExtensionCount,
                        goodAlignBuf,
                        &goodAlignCount,
                        ungappedExtension_new,
                        dp_mem,
                        tree,
                        private_tree,
                        BlastHSPs
                        );
                if(goodAlignCount >= MAX_ALIGNMENTS_PER_QUERY * numQueriesPerThread)
                {
                    printf("ERROR! goodAlignCount = %d >= MAX_ALIGNMENTS_PER_QUERY * numQueriesPerThread\n", goodAlignCount);
                    exit(0);
                }
                if(goodExtensionCount >= MAX_EXTENSIONS_PER_QUERY * numQueriesPerThread)
                {
                    printf("ERROR! goodExtensionCount = %d >= MAX_EXTENSIONS_PER_QUERY * numQueriesPerThread\n", goodExtensionCount);
                    exit(0);
                }
            }
        }
        

        gettimeofday(&end, NULL);
        long hit_time = ((end.tv_sec * 1000000 + end.tv_usec) -
                (start.tv_sec * 1000000 + start.tv_usec));


        Blast_IntervalTreeFree2(private_tree);
        Blast_IntervalTreeFree2(tree);

        struct gappedExtension *gappedExtension_t = gappedExtension_arr + MAX_NUM_GAPEXT_PER_QUERY * thread_id * numQueriesPerThread;

        gettimeofday(&start, NULL);
#if 0
#pragma omp for schedule(dynamic)
        for (jj = 0; jj < numQuery; jj++) {
            alignments_getTracebacks_ncbi_multi(PSSMatrix_arr[jj], 
                    scoreMatrix, &query_composition_arr[jj], jj, 
                    NRrecord[thread_id], querySeq_arr[jj], 
                    matrixInfo_arr + thread_id, gap_align_arr[thread_id], 
                    matrix_arr[thread_id], seq_data_arr[thread_id], 
                    newton_system_arr[thread_id], grads_arr[thread_id], 
                    Scores_arr[thread_id], z_arr[thread_id], 
                    resids_x_arr[thread_id], resids_z_arr[thread_id], 
                    old_scores_arr[thread_id], 
                    workspace_arr[thread_id], 
                    gappedExtension_t, 
                    &gappedExtension_cnt);

            if(gappedExtension_cnt >= MAX_NUM_GAPEXT_PER_QUERY * numQueriesPerThread)
            {
                fprintf(stderr, "gappedExtension_arr overflow\n");
                exit(1);
            }
        }
#endif
        gettimeofday(&end, NULL);

        long traceback_time = ((end.tv_sec * 1000000 + end.tv_usec) -
                (start.tv_sec * 1000000 + start.tv_usec));

        fprintf(stderr, "thread_id: %d Prelim time: %f Traceback time: %f numExt: %d numAlign: %d numGapExt: %d\n", 
                thread_id, (float)hit_time * 1e-6, (float)traceback_time * 1e-6, goodExtensionCount, goodAlignCount, gappedExtension_cnt);

    }
    free(numHitsPerQPos_arr);
    free(numExtHit_arr);
    free(lastHits_arr);
    free(secondBin_arr);
    free(BlastHSPs_arr);
    free(dp_mem_arr);
    free(ungappedExtension_new_arr);

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
        Blast_CompositionWorkspaceFree(&NRrecord[ii]);
        Blast_MatrixInfoFree2(&matrixInfo_arr[ii]);
        //Blast_IntervalTreeFree2(private_tree[ii]);
        //Blast_IntervalTreeFree2(tree[ii]);
    }

    free(workspace_arr);
    free(old_scores_arr);
    free(z_arr);
    free(resids_z_arr);
    free(resids_x_arr);
    free(grads_arr);
    free(Scores_arr);
    free(newton_system_arr);
    free(gap_align_arr);
    free(matrix_arr);
    free(seq_data_arr);

    free(NRrecord);

    free(matrixInfo_arr);

    for (jj = 0; jj < numQuery; jj++)
    {
        free(querySeq_arr[jj]);
    }
    free(query_composition_arr);

    free(querySeq_arr);

    for (ii = 0; ii < numQuery; ii++) {

        printf("Query= %s\n\n", queryDescription_arr[ii]);
        printf("Length=%ld\n\n", strlen(query_arr[ii]));

        //print_info(ii, PSSMatrix_arr);
        print_gappedAlignmentsBrief_multi(ii);
        print_gappedAlignmentsFull_multi(query_arr[ii], PSSMatrix_arr[ii], ii);

#if 0
        printf("Query %d alignment info:\n", ii);
        printf("Query length: %d\n", PSSMatrix_arr[ii].length);
        printf("Number of hits: %d\n", blast_numHits_multi[ii]);
        printf("Number of extensions: %u\n", blast_numUngappedExtensions_multi[ii]);
        printf("Number of successful extensions: %u\n",
                blast_numTriggerExtensions_multi[ii]);
        printf("Number of sequences with successful extensions: %u\n",
                blast_numTriggerSequences_multi[ii]);
        printf("Number of sequences with semi-gapped score above cutoff: %u\n",
                blast_numGoodAlignments_multi[ii]);
        printf("Number of sequences better than %g: %u\n\n", parameters_cutoff,
                alignments_finalAlignments_multi[ii]->numEntries);
        printf("%d numHits: %d numExts: %d numTriggerExts: %d numGapExts: %d numGoodExts: %d numSeqPassed: %d\n", ii, blast_numHits_multi[ii], blast_numUngappedExtensions_multi[ii], blast_numTriggerExtensions_multi[ii], blast_numGappedExtension_multi[ii], blast_numGoodExtensions_multi[ii], blast_numGoodAlignments_multi[ii]);
#endif
        //fasterGappedExtension_free_multi(ii);
        //gappedExtension_free_multi(ii);
        //PSSMatrix_free(PSSMatrix_arr[ii]);
        alignments_free_multi2(ii);
    }

    free(gappedExtension_arr);
    free(goodExtensionBuf_arr);
    free(goodAlignBuf_arr);
}


