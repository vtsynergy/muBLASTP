#include "blast.h"

#define SCALING_FACTOR 32
#define EVALUE_STRETCH 5

#define BLAST_CMP(a,b) ((a)>(b) ? 1 : ((a)<(b) ? -1 : 0))

static int score_compare_match(const void *v1, const void *v2)
{
    struct ungappedExtension *h1, *h2;
    int result = 0;

    h1 = *((struct ungappedExtension **)v1);
    h2 = *((struct ungappedExtension **)v2);

    int4 l1 = h1->end.subjectOffset - h1->start.subjectOffset;
    int4 l2 = h2->end.subjectOffset - h2->start.subjectOffset;

    if (0 == (result = BLAST_CMP(h2->nominalScore,
                    h1->nominalScore)) &&
            0 == (result = BLAST_CMP(h1->start.subjectOffset,
                    h2->start.subjectOffset)) &&
            0 == (result = BLAST_CMP(l2,
                    l1)) &&
            0 == (result = BLAST_CMP(h1->start.queryOffset,
                    h2->start.queryOffset))) {
        result = BLAST_CMP(l2,
                l1);
    }

    return result;
}

static int score_compare_match3(const void *v1, const void *v2)
{
    struct ungappedExtension *h1, *h2;
    int result = 0;

    h1 = ((struct ungappedExtension *)v1);
    h2 = ((struct ungappedExtension *)v2);

    int4 l1 = h1->end.subjectOffset - h1->start.subjectOffset;
    int4 l2 = h2->end.subjectOffset - h2->start.subjectOffset;

    if (0 == (result = BLAST_CMP(h2->nominalScore,
                    h1->nominalScore)) &&
            0 == (result = BLAST_CMP(h1->start.subjectOffset,
                    h2->start.subjectOffset)) &&
            0 == (result = BLAST_CMP(l2,
                    l1)) &&
            0 == (result = BLAST_CMP(h1->start.queryOffset,
                    h2->start.queryOffset))) {
        result = BLAST_CMP(l2,
                l1);
    }

    return result;
}

void alignments_sortUngapedExtension(struct ungappedExtension **ungappedExtensions, int numExts)
{
    qsort(ungappedExtensions,
            numExts,
            sizeof(struct ungappedExtension *), score_compare_match);
}

void alignments_sortUngapedExtension3(struct ungappedExtension *ungappedExtensions, int numExts)
{
    qsort(ungappedExtensions,
            numExts,
            sizeof(struct ungappedExtension), score_compare_match3);
}

void alignments_getTracebacks_ncbi_multi(
        struct PSSMatrix PSSMatrix, 
        struct scoreMatrix scoreMatrix, 
        Blast_AminoAcidComposition *query_composition, 
        int queryNum, 
        Blast_CompositionWorkspace *NRrecord, 
        char *querySeq, 
        Blast_MatrixInfo *matrixInfo,
        BlastGapAlignStruct *gap_align,
        Int4 **matrix, unsigned char *seq_data,
        ReNewtonSystem *newton_system, 
        double **grads, double **Scores, 
        double *z,
        double *resids_x,
        double *resids_z,
        double *old_scores,
        double *workspace, struct gappedExtension *gappedExtension_arr, int *gappedExtension_cnt)
{
    struct finalAlignment *goodAlignment;
    struct alignment *alignment;
    BlastScoringParameters score_params;
    struct ungappedExtension *ungappedExtension;

    alignments_sortGoodAlignments_multi(queryNum);

    score_params.gap_open = 352;
    score_params.gap_extend = 32;
    score_params.shift_pen = 32767;
    score_params.scale_factor = 32;
    int ii, jj;
    double cutoff = 149.53287490731404;
    Boolean subject_maybe_biased = FALSE;
    gap_align->gap_x_dropoff = 2076; 

    s_MatrixInfoInit(matrixInfo, 0.31760595763573063, 32);

    memSingleBlock_resetCurrent(alignments_goodAlignments_multi[queryNum]);

    struct ungappedExtension *hsp_array[MAX_NUM_ALIGNMENTS];

    int hspcnt = 0;
    int num_goodAlignment = 0;

    ASSERT(alignments_finalAlignments_multi[queryNum]->numEntries == 0);
    while ((goodAlignment = memSingleBlock_getCurrent(
                    alignments_goodAlignments_multi[queryNum])) != NULL) {
        alignment = goodAlignment->alignment;

        Boolean oldNearIdenticalStatus = FALSE;

        ungappedExtension = alignment->ungappedExtensions;


        int hsp_index = 0;
        int num_adjustments = 0;

        double ecutoff = 0.002;
        if(num_goodAlignment >= parameters_numDisplayAlignments && alignment->best_eValue > EVALUE_STRETCH * ecutoff)
        {
            continue;
        }

        num_goodAlignment++;
        BlastCompo_SequenceData seqData;
        seqData.buffer = seq_data; 
        memset(seqData.buffer, 0, sizeof(unsigned char) * alignment->subjectLength + 2);
        seqData.length = alignment->subjectLength;
        seqData.data = seqData.buffer + 1;

        for(ii = 0; ii < alignment->subjectLength; ii++)
        {
            seqData.data[ii] = to_ncbi[alignment->subject[ii]];
        }

        while (ungappedExtension != NULL) 
        {
            if(ungappedExtension->status == ungappedExtension_SEMIGAPPED){

                double lambda = 0.008344;
                Boolean isContained = s_IsContained(ungappedExtension, alignment, lambda);

                if(isContained)
                {
                    ungappedExtension->status == ungappedExtension_DELETED;
                    ungappedExtension = ungappedExtension->next;
                    continue;
                }

                Blast_AminoAcidComposition subject_composition;

                ungappedExtension->nominalScore *= SCALING_FACTOR;

                Boolean shouldTestIdentical =  
                    s_preliminaryTestNearIdentical(PSSMatrix.length,
                            ungappedExtension,
                            cutoff);

                if(hsp_index == 0 || (shouldTestIdentical != oldNearIdenticalStatus))
                {
                    if(!shouldTestIdentical)
                    {
                        s_DoSegSequenceData(&seqData, eBlastTypeBlastp,
                                &subject_maybe_biased);
                        num_adjustments++;
                    }
                }



                if(hsp_index == 0 || (shouldTestIdentical != oldNearIdenticalStatus))
                {

                    Blast_CompositionWorkspaceInit(NRrecord,
                            "BLOSUM62");
                    s_GetComposition(&subject_composition, BLASTAA_SIZE, seqData.data, alignment->subjectLength);

                    int compositionTestIndex = 0;
                    int RE_pseudocounts = 20;
                    ECompoAdjustModes compo_adjust_mode = eCompositionMatrixAdjust;
                    double pvalueForThisPair = (-1);
                    double LambdaRatio;
                    EMatrixAdjustRule matrix_adjust_rule = eDontAdjustMatrix;

                    Blast_AdjustScores2(matrix, query_composition,
                            PSSMatrix.length,
                            &subject_composition,
                            alignment->subjectLength,
                            matrixInfo, compo_adjust_mode,
                            RE_pseudocounts, NRrecord,
                            &matrix_adjust_rule,
                            &pvalueForThisPair,
                            compositionTestIndex,
                            &LambdaRatio,
                            newton_system, z, resids_x, resids_z, old_scores, workspace, grads, Scores);

                }

                BLAST_GappedAlignmentWithTraceback(querySeq, seqData.data, gap_align, &score_params, ungappedExtension->gap_start.queryOffset, ungappedExtension->gap_start.subjectOffset, PSSMatrix.length, alignment->subjectLength, NULL, matrix);

                int cutoff_s = 1;
                if(ungappedExtension->nominalScore > cutoff_s)
                {
                    ungappedExtension->start.queryOffset = gap_align->query_start;
                    ungappedExtension->start.subjectOffset = gap_align->subject_start;
                    ungappedExtension->end.queryOffset = gap_align->query_stop;
                    ungappedExtension->end.subjectOffset = gap_align->subject_stop;
                    ungappedExtension->nominalScore = gap_align->score;
                    ungappedExtension->status = ungappedExtension_GAPPED;

#if 1
                    struct trace trace;
                    trace.length = gap_align->fwd_prelim_tback->total_num_ops + gap_align->rev_prelim_tback->total_num_ops;
                    trace.queryStart = gap_align->query_start;
                    trace.subjectStart = gap_align->subject_start;
                    unsigned char* traceCodes;
                    traceCodes = (unsigned char*)global_malloc(sizeof(unsigned char) * trace.length);
                    trace.traceCodes = traceCodes;


                    int ii;
                    int traceCode_cnt = 0;
                    for(ii = 0; ii < gap_align->rev_prelim_tback->num_ops; ii++)
                    {
                        for(jj = 0; jj < gap_align->rev_prelim_tback->edit_ops[ii].num; jj++)
                        {
                            switch(gap_align->rev_prelim_tback->edit_ops[ii].op_type) {
                                case eGapAlignDel: 
                                    traceCodes[traceCode_cnt] = 1; 
                                    break;
                                case eGapAlignSub: 
                                    traceCodes[traceCode_cnt] = 0; 
                                    break;
                                case eGapAlignIns:
                                    traceCodes[traceCode_cnt] = 2; 
                                    break;
                                default:
                                    fprintf(stderr, "Err: alien opt code %d\n", traceCode_cnt);
                                    exit(1);
                            }
                            traceCode_cnt++;
                        }
                    }

                    for(ii = gap_align->fwd_prelim_tback->num_ops - 1; ii >= 0; ii--)
                    {
                        for(jj = 0; jj < gap_align->fwd_prelim_tback->edit_ops[ii].num; jj++)
                        {
                            switch(gap_align->fwd_prelim_tback->edit_ops[ii].op_type) {
                                case eGapAlignDel: 
                                    traceCodes[traceCode_cnt] = 1; 
                                    break;
                                case eGapAlignSub: 
                                    traceCodes[traceCode_cnt] = 0; 
                                    break;
                                case eGapAlignIns:
                                    traceCodes[traceCode_cnt] = 2; 
                                    break;
                                default:
                                    fprintf(stderr, "Err: alien opt code %d\n", traceCode_cnt);
                                    exit(1);
                            }
                            traceCode_cnt++;
                        }
                    }

                    ungappedExtension->trace = trace;
#endif
                }
                else
                {
                    ungappedExtension->status = ungappedExtension_DELETED;
                }
                hsp_index++;
                oldNearIdenticalStatus = shouldTestIdentical;

            }


            ungappedExtension = ungappedExtension->next;

        }


        ungappedExtension = alignment->ungappedExtensions;
        hspcnt = 0;
        while (ungappedExtension != NULL) 
        {

            if(ungappedExtension->status == ungappedExtension_GAPPED)
            {
                hsp_array[hspcnt] = ungappedExtension;
                hspcnt++;
                if(hspcnt > MAX_NUM_ALIGNMENTS)
                {
                    fprintf(stderr, "hspcnt exceed MAX_NUM_ALIGNMENTS\n");
                    exit(1);

                }
            }

            ungappedExtension = ungappedExtension->next;
        }

        int bestScore = 0;
        double best_eValue = 0;
        int hspcnt_new = s_HitlistEvaluateAndPurge(&bestScore, &best_eValue, hsp_array, hspcnt, PSSMatrix.length, alignment->subjectLength);

        for(ii = hspcnt_new; ii < hspcnt; ii++)
        {
            free(hsp_array[ii]->trace.traceCodes);
        }

        if(hspcnt_new)
        {
            int ii;
            for(ii = 0; ii < hspcnt_new; ii++)
            {
                ASSERT(hsp_array[ii]->status == ungappedExtension_GAPPED); 

                struct gappedExtension *gappedExtension = gappedExtension_arr + *gappedExtension_cnt;
                (*gappedExtension_cnt)++;

                gappedExtension->trace = hsp_array[ii]->trace;
                gappedExtension->next = NULL;
                gappedExtension->queryEnd = hsp_array[ii]->end.queryOffset;
                gappedExtension->subjectEnd = hsp_array[ii]->end.subjectOffset;
                gappedExtension->subjectStart = hsp_array[ii]->start.subjectOffset;
                gappedExtension->nominalScore = hsp_array[ii]->nominalScore;
                gappedExtension->eValue = hsp_array[ii]->eValue;

#if 1
                gappedExtension_score_ncbi_multi(gappedExtension, queryNum);
                alignments_addGappedExtension(alignment, gappedExtension);
#endif
            }

            //alignments_addGappedExtension(alignment, gappedExtension);
            struct finalAlignment *finalAlignment = alignments_addFinalAlignment_multi(bestScore, alignment, queryNum);
#if 1
            finalAlignment->description = descriptions_getDescription_mem(
                    finalAlignment->alignment->descriptionLocation,
                    finalAlignment->alignment->descriptionLength);
#endif
        }
    }

}


// Perform initial scoring of all ungapped extensions to find "good" alignments
// that may
// score above the cutoff
int alignments_findGoodAlignments_ncbi_multi3(struct alignment *alignment,
        struct PSSMatrix PSSMatrix,
        struct scoreMatrix scoreMatrix,
        int queryNum,
        struct ungappedExtension **ungappedExtension_new,
        BlastGapDP *dp_mem,
        BlastIntervalTree *tree,
        BlastIntervalTree *private_tree, 
        BlastHSP *BlastHSP_arr) {
    struct ungappedExtension *ungappedExtension;
    int4 bestScore, hasChildren;
    Int4 num_ungappedExt = alignment->numExtensions;

    Int4 BlastHSP_cnt = 0;
    BlastGapAlignStruct gap_align;
    gap_align.gap_x_dropoff = statistics_gappedNominalDropoff_multi[queryNum] ;
    gap_align.dp_mem_alloc = DP_MEM_SIZE;
    gap_align.dp_mem = dp_mem; 

    BlastScoringParameters score_params;
    //score_params.gap_open = parameters_semiGappedOpenGap ;
    score_params.gap_open = 11;
    score_params.gap_extend = parameters_semiGappedExtendGap;


    blast_dloc_multi[queryNum] = alignment->descriptionLocation;

    // Record if subject has children
    if (encoding_alphabetType == encoding_protein &&
            alignment->encodedLength > alignment->subjectLength + 2)
    {
        hasChildren = 1;
        fprintf(stderr, "hasChildren: %d\n", hasChildren);
    }
    else
    {
        hasChildren = 0;
    }

    tree = Blast_IntervalTreeInit2(0, PSSMatrix.length + 1, 0, alignment->subjectLength + 1, tree);

    Boolean restricted_align;

    const double kRestrictedMult = 0.68;

    ungappedExtension = alignment->ungappedExtensions;

    alignments_sortUngapedExtension3(ungappedExtension, num_ungappedExt);

    if (ungappedExtension->nominalScore <
            (Int4)(kRestrictedMult *
                blast_gappedNominalCutoff_multi[queryNum] )) {
        restricted_align = TRUE;
    }
    else {
        restricted_align = FALSE;
    }

    Int4 index = 0;
    Int4 redo_query = -1;
    Int4 redo_index = -1;
    Int4 query_index = queryNum;
    bestScore = 0;
    for(index = 0; index < num_ungappedExt; index++) {
        BlastHSP tmp_hsp;
        tmp_hsp.score = ungappedExtension[index].nominalScore; 
        tmp_hsp.context = 0;
        tmp_hsp.query.offset = ungappedExtension[index].start.queryOffset;
        tmp_hsp.query.end = ungappedExtension[index].end.queryOffset;
        tmp_hsp.query.frame = 0;
        tmp_hsp.subject.offset = ungappedExtension[index].start.subjectOffset;
        tmp_hsp.subject.end = ungappedExtension[index].end.subjectOffset;
        tmp_hsp.subject.frame = 0;

        Int4 min_diag_separation = 0;

        if ((index <= redo_index && 
                    !BlastIntervalTreeContainsHSP(private_tree, &tmp_hsp, min_diag_separation)) ||
                !BlastIntervalTreeContainsHSP(tree, &tmp_hsp, min_diag_separation)) {

            // Find the seed
            ungappedExtension_findSeed(ungappedExtension + index, PSSMatrix,
                    alignment->subject);

            BLAST_SequenceBlk query_blk;
            query_blk.sequence = PSSMatrix.queryCodes;
            query_blk.length = PSSMatrix.length;

            BLAST_SequenceBlk subject_blk;

            subject_blk.sequence = alignment->subject;
            subject_blk.length = alignment->subjectLength;

            Int4 q_start, s_start, q_length, s_length, q_off, s_off;
            q_off = q_start = ungappedExtension[index].start.queryOffset;
            s_off = s_start = ungappedExtension[index].start.subjectOffset;
            s_length = q_length = ungappedExtension[index].end.queryOffset - ungappedExtension[index].start.queryOffset;
            Int4 max_offset = BlastGetStartForGappedAlignment(PSSMatrix.queryCodes, alignment->subject, scoreMatrix.matrix, q_start, q_length, s_start, s_length);
            s_off += max_offset - q_off;
            q_off = max_offset;

            blast_numGappedExtension_multi[queryNum]++;
            s_BlastProtGappedAlignment(&query_blk, q_off, &subject_blk, s_off, &gap_align, &score_params, restricted_align, ungappedExtension + index, scoreMatrix);

            Int4 cutoff, restricted_cutoff = 0;
            cutoff = blast_gappedNominalCutoff_multi[queryNum];

            if (restricted_align)
                restricted_cutoff = (Int4)(kRestrictedMult * cutoff);


            if (restricted_align &&
                    gap_align.score < cutoff &&
                    gap_align.score >= restricted_cutoff) {

                if (!private_tree) {
                    private_tree = Blast_IntervalTreeInit(0, 
                            PSSMatrix.length + 1, 0, alignment->subjectLength + 1);
                }
                else {
                    private_tree = Blast_IntervalTreeInit2(0, PSSMatrix.length + 1, 0, alignment->subjectLength + 1, private_tree);
                }

                redo_index = index;
                redo_query = query_index;
                index = -1;
                restricted_align = FALSE;
                continue;
            }

            if(gap_align.score >= blast_gappedNominalCutoff_multi[queryNum])
            {
                if(gap_align.score > bestScore)
                    bestScore = gap_align.score;

                ungappedExtension[index].start.queryOffset = gap_align.query_start; 
                ungappedExtension[index].start.subjectOffset = gap_align.subject_start; 
                ungappedExtension[index].end.queryOffset = gap_align.query_stop; 
                ungappedExtension[index].end.subjectOffset = gap_align.subject_stop; 

                ungappedExtension[index].gap_start.queryOffset = q_off; 
                ungappedExtension[index].gap_start.subjectOffset = s_off; 
                ungappedExtension[index].nominalScore = gap_align.score;
                ungappedExtension[index].status = ungappedExtension_SEMIGAPPED;


                BlastHSP *new_hsp = Blast_HSPNew2(BlastHSP_arr, &BlastHSP_cnt);
                //BlastHSP *new_hsp = Blast_HSPNew();

                new_hsp->score = ungappedExtension[index].nominalScore; 
                new_hsp->context = 0;
                new_hsp->query.offset = ungappedExtension[index].start.queryOffset;
                new_hsp->query.end = ungappedExtension[index].end.queryOffset;
                new_hsp->query.frame = 0;
                new_hsp->subject.offset = ungappedExtension[index].start.subjectOffset;
                new_hsp->subject.end = ungappedExtension[index].end.subjectOffset;
                new_hsp->subject.frame = 0;

                BlastIntervalTreeAddHSP(new_hsp, tree, eQueryAndSubject);

                if (index < redo_index) {
                    BlastIntervalTreeAddHSP(new_hsp, private_tree, eQueryAndSubject);
                }
            }

        }

    }

    double prelim_evalue = 50;

    if(num_ungappedExt >= 10 * MAX_EXTENSIONS_PER_QUERY)
    {
        printf("ERROR! queryNum = %d, num_ungappedExt = %d >= 10 * MAX_EXTENSIONS_PER_QUERY\n", queryNum, num_ungappedExt);
        exit(0);
    }
    int num_GappedExt = 0;
    int ii;
    for(ii = 0; ii < num_ungappedExt; ii++)
    {
        if(ungappedExtension[ii].status == ungappedExtension_SEMIGAPPED)
        {
            ASSERT(num_GappedExt < 10 * MAX_EXTENSIONS_PER_QUERY)
                ungappedExtension_new[num_GappedExt] = ungappedExtension + ii;
            num_GappedExt++;
        }
    }


    int new_numGappedExt = Blast_HSPListPurgeHSPsWithCommonEndpoints(ungappedExtension_new, num_GappedExt, TRUE);
    int4 numExtensions = 0;
    if(new_numGappedExt > 0)
    {
        //fprintf(stderr, "%d : %d > %d\n", alignment->sequenceCount, num_GappedExt,  new_numGappedExt);
        alignments_get_eValue(ungappedExtension_new, new_numGappedExt, PSSMatrix.length, alignment->subjectLength, &(alignment->best_eValue));
        s_Blast_HSPListReapByPrelimEvalue(ungappedExtension_new, new_numGappedExt, prelim_evalue);

        for(index = 0; index < new_numGappedExt; index++){
            if(ungappedExtension_new[index]->status == ungappedExtension_SEMIGAPPED)
            {
                //fprintf(stderr, "findGoodAlignment output: subject id: %d index: %d\n", alignment->sequenceCount, index);
                numExtensions++;
            }
        } 
    }

    //free(ungappedExtension_new);
    for(ii = 0; ii < num_ungappedExt; ii++)
    {
        ungappedExtension = alignment->ungappedExtensions + ii;
        ungappedExtension->next = alignment->ungappedExtensions + ii + 1;
    }

    ungappedExtension->next = NULL;
    // If this alignment contains gapped extensions that could score above
    // cutoff
    if (numExtensions > 0) {
        // If a single sequence add to list of "good" alignments
        alignments_addGoodAlignment_multi(bestScore, alignment, queryNum);

        blast_numGoodExtensions_multi[queryNum] += numExtensions;
        blast_numGoodAlignments_multi[queryNum]++;

        //fprintf(stderr, "%d %d\n", alignment->sequenceCount, numExtensions);
        return 1;
    }
    else {
        return 0;
    }
}

