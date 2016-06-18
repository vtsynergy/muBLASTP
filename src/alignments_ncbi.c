#include "blast.h"
#include <sys/time.h>

#define SCALING_FACTOR 32
#define EVALUE_STRETCH 5

#define BLAST_CMP(a,b) ((a)>(b) ? 1 : ((a)<(b) ? -1 : 0))

void gappedExtension_score_ncbi_multi(struct gappedExtension *gappedExtension,
                                 int queryNum) {
  gappedExtension->normalizedScore =
      statistics_gappedNominal2normalized_ncbi(gappedExtension->nominalScore);

  gappedExtension->eValue = gappedExtension->eValue;
}


// Add a high-scoring gapped extension to this alignment
void alignments_addGappedExtension2(struct alignment *alignment,
        struct gappedExtension *newExtension, 
        int4 newExtensionOffset, struct gappedExtension *gappedExtensionBuf) 
{
    struct gappedExtension *currentExtension, *previousExtension;
    int4 currentExtensionOffset, previousExtensionOffset;

    // If this is the first high-scoring gapped extension for this alignment
    if (alignment->gappedExtensionOffset == -1) {
        // Make this the first gapped extension in the alignment
        //alignment->gappedExtensions = newExtension;
        alignment->gappedExtensionOffset = newExtensionOffset;
    } else {
        // Start at beginning of list of gapped extensions (one with highest score)
        currentExtensionOffset = alignment->gappedExtensionOffset;
        previousExtensionOffset = -1;

        currentExtension = gappedExtensionBuf + currentExtensionOffset;

        // Move through list of existing extensions until we either reach the
        // end or reach one with a score less than newExtension
        while (currentExtensionOffset != -1 &&
                (currentExtension->nominalScore > newExtension->nominalScore)) {

            previousExtensionOffset = currentExtensionOffset;
            previousExtension = gappedExtensionBuf + previousExtensionOffset;

            currentExtensionOffset = currentExtension->nextOffset;
            currentExtension = gappedExtensionBuf + currentExtensionOffset;
        }

        if (previousExtensionOffset == -1) {
            // This is the highest scoring extension, insert at front of the queue
            alignment->gappedExtensionOffset = newExtensionOffset;
            newExtension->nextOffset = currentExtensionOffset;
        } else {
            // Insert between higher and lower scoring extensions
            previousExtension->nextOffset = newExtensionOffset;
            newExtension->nextOffset = currentExtensionOffset;
        }
    }
}

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

// Perform initial scoring of all ungapped extensions to find "good" alignments
// that may
// score above the cutoff
int alignments_findGoodAlignments_ncbi(
        struct alignment *alignment,
        struct ungappedExtension *goodExtensionBuf,
        struct PSSMatrix PSSMatrix,
        struct scoreMatrix scoreMatrix,
        int queryNum,
        struct ungappedExtension **ungappedExtension_new,
        BlastGapDP *dp_mem,
        BlastIntervalTree *tree,
        BlastIntervalTree *private_tree, 
        BlastHSP *BlastHSP_arr) 
{
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
    //if (encoding_alphabetType == encoding_protein &&
            //alignment->encodedLength > alignment->subjectLength + 2)
    //{
        //hasChildren = 1;
    //}
    //else
    //{
        //hasChildren = 0;
    //}

    tree = Blast_IntervalTreeInit2(0, PSSMatrix.length + 1, 0, alignment->subjectLength + 1, tree);

    Boolean restricted_align;

    const double kRestrictedMult = 0.68;

    ungappedExtension = goodExtensionBuf + alignment->ungappedExtensionOffset;

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
            //ungappedExtension_findSeed(ungappedExtension + index, PSSMatrix,
            //alignment->subject);

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

            //blast_numGappedExtension_multi[queryNum]++;
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


                BlastHSP *new_hsp = Blast_HSPInit(BlastHSP_arr, &BlastHSP_cnt);
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

    int num_GappedExt = 0;
    int ii;
    for(ii = 0; ii < num_ungappedExt; ii++)
    {
        if(ungappedExtension[ii].status == ungappedExtension_SEMIGAPPED)
        {
            ASSERT(num_GappedExt < MAX_EXTENSIONS_PER_QUERY);
            ungappedExtension_new[num_GappedExt] = ungappedExtension + ii;
            num_GappedExt++;
        }
    }


    int new_numGappedExt = Blast_HSPListPurgeHSPsWithCommonEndpoints(ungappedExtension_new, num_GappedExt, TRUE);
    int4 numExtensions = 0;
    if(new_numGappedExt > 0)
    {
        //fprintf(stderr, "%d : %d > %d\n", alignment->sequenceCount, num_GappedExt,  new_numGappedExt);
        alignments_get_eValue(ungappedExtension_new, 
                new_numGappedExt, PSSMatrix.length, 
                alignment->subjectLength, 
                &(alignment->best_eValue));
        s_Blast_HSPListReapByPrelimEvalue(ungappedExtension_new, 
                new_numGappedExt, 
                prelim_evalue);

        for(index = 0; index < new_numGappedExt; index++){
            if(ungappedExtension_new[index]->status == ungappedExtension_SEMIGAPPED)
            {
                //fprintf(stderr, "findGoodAlignment output: subject id: %d index: %d\n", alignment->sequenceCount, index);
                numExtensions++;
            }
        } 
    }

    //free(ungappedExtension_new);
    //for(ii = 0; ii < num_ungappedExt; ii++)
    //{
        //ungappedExtension = alignment->ungappedExtensions + ii;
        //ungappedExtension->next = alignment->ungappedExtensions + ii + 1;
    //}

    //ungappedExtension->next = NULL;
    // If this alignment contains gapped extensions that could score above
    // cutoff
    if (numExtensions > 0) {
        //alignments_addGoodAlignment_multi(bestScore, alignment, queryNum);
        return 1;
    }
    else {
        return 0;
    }
}

int4 alignments_compareGoodAlignments(const void *alignment1,
                                       const void *alignment2) {
  const struct alignment *a1, *a2;

  a1 = (struct alignment *)alignment1;
  a2 = (struct alignment *)alignment2;

  if (a1->best_eValue < a2->best_eValue) {
    return -1;
  } else if (a1->best_eValue > a2->best_eValue) {
    return 1;
  } else {
    // Resolve conflicts using subject length
    if (a1->subjectLength > a2->subjectLength)
      return 1;
    else if (a1->subjectLength < a2->subjectLength)
      return -1;
    else
      return 0;
  }
}

alignments_sortGoodAlignments_multi(struct alignment *goodAlignQuery,
        int4 numGoodAlignQuery)
{
    qsort(goodAlignQuery, numGoodAlignQuery, sizeof(struct alignment),
            alignments_compareGoodAlignments);
}

int EarlyTermination(double evalue, int numAlign, double worse_eValue)
{
    double ecutoff = 0.002;

    if(numAlign >= parameters_numDisplayAlignments && 
            worse_eValue <= evalue)
    {
        if(evalue <= EVALUE_STRETCH * ecutoff)
            return FALSE;
    }
    else{
        return FALSE;
    }

    return TRUE;
}

void alignments_getTracebacks_ncbi(
        int tid,
        struct alignment *goodAlignQuery,
        int4 numGoodAlignQuery,
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
        double *workspace, 
        struct gappedExtension **gappedExtensionBuf, 
        int4 *gappedExtensionCount,
        int4 *gapExtensionBufSize,
        unsigned char **traceCodeBuf,
        int4 *traceCodeCount,
        int4 *traceCodeBufSize,
        Uint8 *query_words
        )
{
    struct finalAlignment *goodAlignment;
    struct alignment *alignment;
    BlastScoringParameters score_params;

    struct timeval start, end;

    alignments_sortGoodAlignments_multi(goodAlignQuery, numGoodAlignQuery);

    score_params.gap_open = 352;
    score_params.gap_extend = 32;
    score_params.shift_pen = 32767;
    score_params.scale_factor = 32;
    int ii, jj;
    double cutoff = 149.53287490731404;
    Boolean subject_maybe_biased = FALSE;
    gap_align->gap_x_dropoff = 2076; 

    s_MatrixInfoInit(matrixInfo, 0.31760595763573063, 32);

    //memSingleBlock_resetCurrent(alignments_goodAlignments_multi[queryNum]);

    int gappedExtensionCount_t = (*gappedExtensionCount);

    int hspcnt = 0;
    int num_goodAlignment = 0;
    //double worse_eValue = INT2_MAX;

    int alignId;
    for(alignId = 0; alignId < numGoodAlignQuery; alignId++)
    {
        alignment = &goodAlignQuery[alignId];

        Boolean oldNearIdenticalStatus = FALSE;

        //fprintf(stderr, "trace: seqId: %d queryId: %d extOff: %d numExt: %d\n",
        //alignment->sequenceCount, alignment->queryCount, 
        //alignment->ungappedExtensionOffset,
        //alignment->numExtensions);

        int hsp_index = 0;
        int num_adjustments = 0;

        //if(EarlyTermination(
                    //alignment->best_eValue, 
                    //num_goodAlignment, 
                    //worse_eValue))
        //{
            //fprintf(stderr, "seqId: %d worse_eValue: %f best_eValue: %f skip\n", alignment->sequenceCount, worse_eValue, alignment->best_eValue);
            //break;
        //}

        BlastCompo_SequenceData seqData;
        seqData.buffer = seq_data; 
        for(ii = 0; ii < (alignment->subjectLength + 2); ii++)
        {
            seqData.buffer[ii] = 0;
        }
        //memset(seqData.buffer, 0, sizeof(unsigned char) * alignment->subjectLength + 2);
        seqData.length = alignment->subjectLength;
        seqData.data = seqData.buffer + 1;

        //ASSERT(alignment->volumnNumber == readdb_volume);
        //fprintf(stderr, "subject: %p\n", alignment->subject);

        for(ii = 0; ii < alignment->subjectLength; ii++)
        {
            seqData.data[ii] = to_ncbi[alignment->subject[ii]];
        }

        struct ungappedExtension *ungappedExtension;
        ungappedExtension = alignment->ungappedExtensions;

        int extCount = 0;
        for(extCount = 0; extCount < alignment->numExtensions; extCount++)
        {

            if(ungappedExtension == NULL)
            {
                break;
            }
            //ASSERT(ungappedExtension != NULL);

            //fprintf(stderr, "ext status: %d\n", ungappedExtension->status);
            //fprintf(stderr, "ungappedExtension: %d\n", ungappedExtension->status);

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
                    if(!shouldTestIdentical || (shouldTestIdentical) 
                            && (!s_TestNearIdentical(&seqData, 0, 
                                    querySeq, 0, query_words, 
                                    ungappedExtension)))
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
                    s_GetComposition(&subject_composition, 
                            BLASTAA_SIZE, seqData.data, alignment->subjectLength);

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
                            newton_system, z, 
                            resids_x, resids_z, old_scores, 
                            workspace, grads, Scores);

                }

                BLAST_GappedAlignmentWithTraceback(
                        tid,
                        querySeq, 
                        seqData.data, 
                        gap_align, &score_params, 
                        ungappedExtension->gap_start.queryOffset, 
                        ungappedExtension->gap_start.subjectOffset, 
                        PSSMatrix.length, alignment->subjectLength, NULL, matrix);

                int cutoff_s = 1;
                if(ungappedExtension->nominalScore > cutoff_s)
                {
                    num_goodAlignment++;
                    ungappedExtension->start.queryOffset = gap_align->query_start;
                    ungappedExtension->start.subjectOffset = gap_align->subject_start;
                    ungappedExtension->end.queryOffset = gap_align->query_stop;
                    ungappedExtension->end.subjectOffset = gap_align->subject_stop;
                    ungappedExtension->nominalScore = gap_align->score;
                    ungappedExtension->status = ungappedExtension_GAPPED;

                    struct trace trace;
                    trace.length = 
                        gap_align->fwd_prelim_tback->total_num_ops + 
                        gap_align->rev_prelim_tback->total_num_ops;

                    trace.queryStart = gap_align->query_start;
                    trace.subjectStart = gap_align->subject_start;
                    trace.traceCodeOff = *traceCodeCount;

                    if(*traceCodeCount + trace.length >= *traceCodeBufSize)
                    {

                        gettimeofday(&start, NULL);
                        *traceCodeBufSize *= 2;
                        *traceCodeBuf = (unsigned char *)global_realloc(*traceCodeBuf, *traceCodeBufSize);
                        gettimeofday(&end, NULL);

                        long malloc_time = ((end.tv_sec * 1000000 + end.tv_usec) -
                                (start.tv_sec * 1000000 + start.tv_usec));

                        //fprintf(stderr, "traceCodeBuf resize to %d (%ld ms)\n", 
                        //*traceCodeBufSize,
                        //malloc_time);
                    }

                    unsigned char* traceCodes = *traceCodeBuf + *traceCodeCount;
                    *traceCodeCount += trace.length;

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
                }
                else
                {
                    ungappedExtension->status = ungappedExtension_DELETED;
                }
                hsp_index++;
                oldNearIdenticalStatus = shouldTestIdentical;

            }

            ungappedExtension++;
        }

        ungappedExtension = alignment->ungappedExtensions;
        hspcnt = 0;
        struct ungappedExtension *hsp_array[alignment->numExtensions];
        for(extCount = 0; extCount < alignment->numExtensions; extCount++)
        {
            if(ungappedExtension == NULL)
            {
                break;
            }

            if(ungappedExtension->status == ungappedExtension_GAPPED)
            {
                hsp_array[hspcnt] = ungappedExtension;
                hspcnt++;
            }

            ungappedExtension++;
        }

        int bestScore = 0;
        double best_eValue = 0;
        int hspcnt_new = s_HitlistEvaluateAndPurge(
                &bestScore, &best_eValue, hsp_array, 
                hspcnt, PSSMatrix.length, alignment->subjectLength);


        if(hspcnt_new)
        {
            int ii;
            for(ii = 0; ii < hspcnt_new; ii++)
            {
                ASSERT(hsp_array[ii]->status == ungappedExtension_GAPPED); 

                if(gappedExtensionCount_t + 1 >= *gapExtensionBufSize)
                {

                    gettimeofday(&start, NULL);
                    *gapExtensionBufSize *= 2;
                    //fprintf(stderr, "gappedExtensionBuf resized to %d\n", *gapExtensionBufSize);
                    *gappedExtensionBuf = (struct gappedExtension *)
                        global_realloc(*gappedExtensionBuf, 
                                sizeof(struct gappedExtension) * (*gapExtensionBufSize));

                    gettimeofday(&end, NULL);

                    long malloc_time = ((end.tv_sec * 1000000 + end.tv_usec) -
                            (start.tv_sec * 1000000 + start.tv_usec));

                    //fprintf(stderr, "traceCodeBuf resize to %d (%ld ms)\n", 
                    //*traceCodeBufSize,
                    //malloc_time);
                }

                struct gappedExtension *gappedExtension = *gappedExtensionBuf + gappedExtensionCount_t;

                gappedExtension->trace = hsp_array[ii]->trace;
                gappedExtension->next = NULL;
                gappedExtension->nextOffset = -1;
                gappedExtension->queryEnd = hsp_array[ii]->end.queryOffset;
                gappedExtension->subjectEnd = hsp_array[ii]->end.subjectOffset;
                gappedExtension->subjectStart = hsp_array[ii]->start.subjectOffset;
                gappedExtension->nominalScore = hsp_array[ii]->nominalScore;
                gappedExtension->eValue = hsp_array[ii]->eValue;

                if(gappedExtension->eValue < 1.0e-180)
                   gappedExtension->eValue = 0.0;

                gappedExtension_score_ncbi_multi(gappedExtension, queryNum);
                alignments_addGappedExtension2(alignment, gappedExtension, 
                        gappedExtensionCount_t, *gappedExtensionBuf);

                gappedExtensionCount_t++;
            }

            struct finalAlignment *finalAlignment;
            finalAlignment = 
                alignments_addFinalAlignment_multi(bestScore, alignment, queryNum);
            finalAlignment->thread_id = tid;

            //finalAlignment->description = descriptions_getDescription_mem(
                    //finalAlignment->alignment->descriptionLocation,
                    //finalAlignment->alignment->descriptionLength);
        }
    }

    (*gappedExtensionCount) = gappedExtensionCount_t;
}
