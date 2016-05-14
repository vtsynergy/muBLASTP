/* ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 */


#include "blast.h"
#include <math.h>

#define NCBIMATH_LN2 0.69314718055994530941723212145818

//#ifndef INT4_MIN
///** Smallest (most negative) number represented by signed int */
//#define INT4_MIN    (-2147483647-1)
//#endif

#define MININT INT4_MIN/2

#define DBL_MAX 1.79769e+308

#define TRUE 1
#define FALSE 0

#ifndef MAX
/** returns larger of a and b. */
#define MAX(a,b)    ((a)>=(b)?(a):(b))
#endif

#ifndef MIN
/** returns smaller of a and b. */
#define MIN(a,b)    ((a)>(b)?(b):(a))
#endif


#define FENCE_SENTRY 201


Int4 edit_script_num_rows_t[MAX_NUM_THREADS];
Uint1** edit_script_t[MAX_NUM_THREADS];
Int4 *edit_start_offset_t[MAX_NUM_THREADS];

/** Values for the editing script operations in traceback */
enum {
    SCRIPT_SUB           = eGapAlignSub,     /**< Substitution */
    SCRIPT_GAP_IN_A      = eGapAlignDel,     /**< Deletion */
    SCRIPT_GAP_IN_B      = eGapAlignIns,     /**< Insertion */
    SCRIPT_OP_MASK       = 0x07, /**< Mask for edit script operations */

    SCRIPT_EXTEND_GAP_A  = 0x10, /**< continue a gap in A */
    SCRIPT_EXTEND_GAP_B  = 0x40  /**< continue a gap in B */
};

    static BlastCompo_GappingParams *
s_GappingParamsNew(BlastKappa_GappingParamsContext * context,
        const BlastExtensionParameters* extendParams,
        int num_queries)
{
    int i;
    double min_lambda = DBL_MAX;   /* smallest gapped Lambda */
    const BlastScoringParameters * scoring = context->scoringParams;
    const BlastExtensionOptions * options = extendParams->options;
    /* The new object */
    BlastCompo_GappingParams * gapping_params = NULL;

    gapping_params = malloc(sizeof(BlastCompo_GappingParams));
    if (gapping_params != NULL) {
        gapping_params->gap_open = scoring->gap_open;
        gapping_params->gap_extend = scoring->gap_extend;
        gapping_params->context = context;
    }

    for (i = 0;  i < num_queries;  i++) {
        if (context->sbp->kbp_gap[i] != NULL &&
                context->sbp->kbp_gap[i]->Lambda < min_lambda) {
            min_lambda = context->sbp->kbp_gap[i]->Lambda;
        }
    }
    gapping_params->x_dropoff = (Int4)
        MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
                extendParams->gap_x_dropoff_final);
    context->gap_align->gap_x_dropoff = gapping_params->x_dropoff;

    return gapping_params;
}

    void
GapPrelimEditBlockReset(GapPrelimEditBlock *edit_block)
{
    if (edit_block) {
        edit_block->num_ops = 0;
        edit_block->total_num_ops = 0;
        edit_block->last_op = eGapAlignInvalid;
    }
}

    static Boolean
s_GapPurgeState(GapStateArrayStruct* state_struct)
{
    while (state_struct)
    {
        /*
           memset(state_struct->state_array, 0, state_struct->used);
           */
        state_struct->used = 0;
        state_struct = state_struct->next;
    }

    return TRUE;
}

    GapEditScript*
GapEditScriptDelete(GapEditScript* old)
{

    if (old)
    {
        free(old->op_type);
        free(old->num);
        free(old);
    }
    return old;
}

    GapPrelimEditBlock *
GapPrelimEditBlockFree(GapPrelimEditBlock *edit_block)
{
    if (edit_block == NULL)
        return NULL;

    free(edit_block->edit_ops);
    free(edit_block);
    return NULL;
}

    GapStateArrayStruct* 
GapStateFree(GapStateArrayStruct* state_struct)

{
    GapStateArrayStruct* next;

    while (state_struct) {
        next = state_struct->next;
        free(state_struct->state_array);
        free(state_struct);
        state_struct = next;
    }

    return NULL;
}




/* Documented in blast_gapalign.h */
    BlastGapAlignStruct* 
BLAST_GapAlignStructFree(BlastGapAlignStruct* gap_align)
{
    if (!gap_align)
        return NULL;

    GapEditScriptDelete(gap_align->edit_script);
    GapPrelimEditBlockFree(gap_align->fwd_prelim_tback);
    GapPrelimEditBlockFree(gap_align->rev_prelim_tback);
    //if (gap_align->greedy_align_mem)
    //s_BlastGreedyAlignsFree(gap_align->greedy_align_mem);
    GapStateFree(gap_align->state_struct);
    free(gap_align->dp_mem);

    free(gap_align);
    return NULL;
}

#define MAX_DBSEQ_LEN 5000000 

    static Int2 
s_GapPrelimEditBlockRealloc(GapPrelimEditBlock *edit_block, Int4 total_ops)
{
    if (edit_block->num_ops_allocated <= total_ops) {
        Int4 new_size = total_ops * 2;
        GapPrelimEditScript *new_ops;

        new_ops = realloc(edit_block->edit_ops, new_size * 
                sizeof(GapPrelimEditScript));
        if (new_ops == NULL)
            return -1;

        edit_block->edit_ops = new_ops;
        edit_block->num_ops_allocated = new_size;
    }
    return 0;
}



    GapPrelimEditBlock *
GapPrelimEditBlockNew(void)
{
    GapPrelimEditBlock *edit_block = malloc(sizeof(GapPrelimEditBlock));
    if (edit_block != NULL) {
        edit_block->edit_ops = NULL;
        edit_block->num_ops_allocated = 0;
        edit_block->num_ops = 0;
        edit_block->last_op = eGapAlignInvalid;
        s_GapPrelimEditBlockRealloc(edit_block, 100);
    }
    return edit_block;
}

Int2
BLAST_GapAlignStructNew( 
        //const BlastExtensionParameters* ext_params, 
        Uint4 max_subject_length,
        //BlastScoreBlk* sbp, 
        BlastGapAlignStruct** gap_align_ptr)
{
    Int2 status = 0;
    BlastGapAlignStruct* gap_align;

    gap_align = (BlastGapAlignStruct*) calloc(1, sizeof(BlastGapAlignStruct));

    *gap_align_ptr = gap_align;

    if (1) {
        /* allocate structures for ordinary dynamic programming */
        gap_align->dp_mem_alloc = 1000;
        gap_align->dp_mem = (BlastGapDP *)malloc(gap_align->dp_mem_alloc *
                sizeof(BlastGapDP));
        if (!gap_align->dp_mem)
            gap_align = BLAST_GapAlignStructFree(gap_align);
    }
    else {
        fprintf(stderr,"Greedy is not supported\n");
    }

    if (!gap_align)
        return -1;

    gap_align->fwd_prelim_tback = GapPrelimEditBlockNew();
    gap_align->rev_prelim_tback = GapPrelimEditBlockNew();

    return status;
}

#define CHUNKSIZE   2097152

    static GapStateArrayStruct*
s_GapGetState(GapStateArrayStruct** head, Int4 length)

{
    GapStateArrayStruct* retval,* var,* last;
    Int4 chunksize = MAX(CHUNKSIZE, length + length/3);

    length += length/3;  /* Add on about 30% so the end will get reused. */
    retval = NULL;
    if (*head == NULL) {
        retval = (GapStateArrayStruct*) 
            malloc(sizeof(GapStateArrayStruct));
        retval->state_array = (Uint1*) malloc(chunksize*sizeof(Uint1));
        retval->length = chunksize;
        retval->used = 0;
        retval->next = NULL;
        *head = retval;
    } else {
        var = *head;
        last = *head;
        while (var) {
            if (length < (var->length - var->used)) {
                retval = var;
                break;
            } else if (var->used == 0) { 
                /* If it's empty and too small, replace. */
                free(var->state_array);
                var->state_array = (Uint1*) malloc(chunksize*sizeof(Uint1));
                var->length = chunksize;
                retval = var;
                break;
            }
            last = var;
            var = var->next;
        }

        if (var == NULL)
        {
            retval = (GapStateArrayStruct*) malloc(sizeof(GapStateArrayStruct));
            retval->state_array = (Uint1*) malloc(chunksize*sizeof(Uint1));
            retval->length = chunksize;
            retval->used = 0;
            retval->next = NULL;
            last->next = retval;
        }
    }

#ifdef ERR_POST_EX_DEFINED
    if (retval->state_array == NULL)
        ErrPostEx(SEV_ERROR, 0, 0, "state array is NULL");
#endif

    return retval;

}




    static Int2 
s_GapPrelimEditBlockAddNew(GapPrelimEditBlock *edit_block, 
        EGapAlignOpType op_type, Int4 num_ops)
{
    if (s_GapPrelimEditBlockRealloc(edit_block, edit_block->num_ops + 2) != 0)
        return -1;

    //ASSERT(op_type != eGapAlignInvalid);

    edit_block->last_op = op_type;
    edit_block->edit_ops[edit_block->num_ops].op_type = op_type;
    edit_block->edit_ops[edit_block->num_ops].num = num_ops;
    edit_block->num_ops++;

    return 0;
}

    void
GapPrelimEditBlockAdd(GapPrelimEditBlock *edit_block, 
        EGapAlignOpType op_type, Int4 num_ops)
{
    if (num_ops == 0)
        return;

    if (edit_block->last_op == op_type)
        edit_block->edit_ops[edit_block->num_ops-1].num += num_ops;
    else
        s_GapPrelimEditBlockAddNew(edit_block, op_type, num_ops);
}


    Int4
ALIGN_EX(int tid, const Uint1* A, const Uint1* B, Int4 M, Int4 N, Int4* a_offset, 
        Int4* b_offset, GapPrelimEditBlock *edit_block, 
        BlastGapAlignStruct* gap_align, 
        const BlastScoringParameters* score_params, Int4 query_offset, 
        Boolean reversed, Boolean reverse_sequence,
        Int4 **matrix,
        Boolean * fence_hit)
{
    /* See Blast_SemiGappedAlign for more general comments on 
       what this code is doing; comments in this function
       only apply to the traceback computations */

    Int4 i; 
    Int4 a_index;
    Int4 b_index, b_size, first_b_index, last_b_index, b_increment;
    const Uint1* b_ptr = NULL;

    BlastGapDP* score_array;

    Int4 gap_open;
    Int4 gap_extend;
    Int4 gap_open_extend;
    Int4 x_dropoff;
    Int4 best_score;

    //Int4** matrix = NULL;
    //Int4** pssm = NULL;
    Int4* matrix_row = NULL;

    Int4 score;
    Int4 score_gap_row;
    Int4 score_gap_col;
    Int4 next_score;

    GapStateArrayStruct* state_struct;
    Uint1* edit_script_row;
    //Uint1** edit_script_t[tid];
    //Int4 *edit_start_offset_t[tid];
    //Int4 edit_script_num_rows_t[tid];
    Int4 orig_b_index;
    Uint1 script, next_script, script_row, script_col;
    Int4 num_extra_cells;

    //matrix = gap_align->sbp->matrix->data;
    //if (gap_align->positionBased) {
    //pssm = gap_align->sbp->psi_matrix->pssm->data;
    //}

    *a_offset = 0;
    *b_offset = 0;
    gap_open = score_params->gap_open;
    gap_extend = score_params->gap_extend;
    gap_open_extend = gap_open + gap_extend;
    x_dropoff = gap_align->gap_x_dropoff;

    if (x_dropoff < gap_open_extend)
        x_dropoff = gap_open_extend;

    if(N <= 0 || M <= 0) 
        return 0;

    /* Initialize traceback information. edit_script_t[tid][] is
       a 2-D array which is filled in row by row as the 
       dynamic programming takes place */

    s_GapPurgeState(gap_align->state_struct);

    /* Conceptually, traceback requires a byte array of size
       MxN. To save on memory, each row of the array is allocated
       separately. edit_script_t[tid][i] points to the storage reserved
       for row i, and edit_start_offset_t[tid][i] gives the offset into
       the B sequence corresponding to element 0 of edit_script_t[tid][i].

       Also make the number of edit script rows grow dynamically */

    //edit_script_num_rows_t[tid] = EDIT_SCRIPT_MAX_NUM_ROWS;
    //edit_script_t[tid] = (Uint1**) malloc(sizeof(Uint1*) * edit_script_num_rows_t[tid]);
    //edit_start_offset_t[tid] = (Int4*) malloc(sizeof(Int4) * edit_script_num_rows_t[tid]);

    /* allocate storage for the first row of the traceback
       array. Because row elements correspond to gaps in A,
       the alignment can only go x_dropoff/gap_extend positions
       at most before failing the X dropoff criterion */

    if (gap_extend > 0)
        num_extra_cells = x_dropoff / gap_extend + 3;
    else
        num_extra_cells = N + 3;

    if (num_extra_cells > gap_align->dp_mem_alloc) {
        gap_align->dp_mem_alloc = MAX(num_extra_cells + 100,
                2 * gap_align->dp_mem_alloc);
        free(gap_align->dp_mem);
        gap_align->dp_mem = (BlastGapDP *)malloc(gap_align->dp_mem_alloc *
                sizeof(BlastGapDP));
    }

    state_struct = s_GapGetState(&gap_align->state_struct, num_extra_cells);

    edit_script_t[tid][0] = state_struct->state_array;
    edit_start_offset_t[tid][0] = 0;
    edit_script_row = state_struct->state_array;

    score = -gap_open_extend;
    score_array = gap_align->dp_mem;

    //memset(score_array, 0, sizeof(BlastGapDP) * N);
    score_array[0].best = 0;
    score_array[0].best_gap = -gap_open_extend;

    for (i = 1; i <= N; i++) {
        if (score < -x_dropoff) 
            break;

        score_array[i].best = score;
        score_array[i].best_gap = score - gap_open_extend; 
        score -= gap_extend;
        edit_script_row[i] = SCRIPT_GAP_IN_A;
    }
    state_struct->used = i + 1;

    b_size = i;
    best_score = 0;
    first_b_index = 0;
    if (reverse_sequence)
        b_increment = -1;
    else
        b_increment = 1;

    for (a_index = 1; a_index <= M; a_index++) {

        /* Set up the next row of the edit script; this involves
           allocating memory for the row, then pointing to it.
           It is not necessary to allocate space for offsets less
           than first_b_index (since the inner loop will never 
           look at them); 

           It is unknown at this point how far to the right the 
           current traceback row will extend; all that's known for
           sure is that the previous row fails the X-dropoff test
           after b_size cells, and that the current row can go at
           most num_extra_cells beyond that before failing the test */

        if (gap_extend > 0)
            state_struct = s_GapGetState(&gap_align->state_struct, 
                    b_size - first_b_index + num_extra_cells);
        else
            state_struct = s_GapGetState(&gap_align->state_struct, 
                    N + 3 - first_b_index);

        if (a_index == edit_script_num_rows_t[tid]) {
            edit_script_num_rows_t[tid] = edit_script_num_rows_t[tid] * 2;
            edit_script_t[tid] = (Uint1 **)realloc(edit_script_t[tid], 
                    edit_script_num_rows_t[tid] *
                    sizeof(Uint1 *));
            edit_start_offset_t[tid] = (Int4 *)realloc(edit_start_offset_t[tid], 
                    edit_script_num_rows_t[tid] *
                    sizeof(Int4));
        }

        edit_script_t[tid][a_index] = state_struct->state_array + 
            state_struct->used + 1;
        edit_start_offset_t[tid][a_index] = first_b_index;

        /* the inner loop assumes the current traceback
           row begins at offset zero of B */

        edit_script_row = edit_script_t[tid][a_index] - first_b_index;
        orig_b_index = first_b_index;

        //if (!(gap_align->positionBased))
        {
            if(reverse_sequence)
                matrix_row = matrix[ A[ M - a_index ] ];
            else
                matrix_row = matrix[ A[ a_index ] ];
        }


        if(reverse_sequence)
            b_ptr = &B[N - first_b_index];
        else
            b_ptr = &B[first_b_index];

        score = MININT;
        score_gap_row = MININT;
        last_b_index = first_b_index;


        for (b_index = first_b_index; b_index < b_size; b_index++) {

            b_ptr += b_increment;
            score_gap_col = score_array[b_index].best_gap;

            next_score = score_array[b_index].best + matrix_row[ *b_ptr ];

            /* script, script_row and script_col contain the
               actions specified by the dynamic programming.
               when the inner loop has finished, 'script' con-
               tains all of the actions to perform, and is
               written to edit_script_t[tid][a_index][b_index]. Otherwise,
               this inner loop is exactly the same as the one
               in Blast_SemiGappedAlign() */

            script = SCRIPT_SUB;
            script_col = SCRIPT_EXTEND_GAP_B;
            script_row = SCRIPT_EXTEND_GAP_A;

            if (score < score_gap_col) {
                script = SCRIPT_GAP_IN_B;
                score = score_gap_col;
            }
            if (score < score_gap_row) {
                script = SCRIPT_GAP_IN_A;
                score = score_gap_row;
            }

            if (best_score - score > x_dropoff) {

                if (first_b_index == b_index)
                    first_b_index++;
                else
                    score_array[b_index].best = MININT;
            }
            else {
                last_b_index = b_index;
                if (score > best_score) {
                    best_score = score;
                    *a_offset = a_index;
                    *b_offset = b_index;
                }

                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;
                if (score_gap_col < (score - gap_open_extend)) {
                    score_array[b_index].best_gap = score - gap_open_extend;
                }
                else {
                    score_array[b_index].best_gap = score_gap_col;
                    script += script_col;
                }

                if (score_gap_row < (score - gap_open_extend)) 
                    score_gap_row = score - gap_open_extend;
                else
                    script += script_row;

                score_array[b_index].best = score;
            }

            score = next_score;
            edit_script_row[b_index] = script;
        }

        if (first_b_index == b_size || (fence_hit && *fence_hit))
            break;

        if (last_b_index + num_extra_cells + 3 >= gap_align->dp_mem_alloc) {

            gap_align->dp_mem_alloc = MAX(last_b_index + num_extra_cells + 100,
                    2 * gap_align->dp_mem_alloc);
            score_array = (BlastGapDP *)realloc(score_array,
                    gap_align->dp_mem_alloc *
                    sizeof(BlastGapDP));
            gap_align->dp_mem = score_array;
        }


        if (last_b_index < b_size - 1) {
            b_size = last_b_index + 1;
        }
        else {
            while (score_gap_row >= (best_score - x_dropoff) && b_size <= N) {

                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                edit_script_row[b_size] = SCRIPT_GAP_IN_A;
                b_size++;
            }
        }

        /* update the memory allocator to reflect the exact number
           of traceback cells this row needed */

        state_struct->used += MAX(b_index, b_size) - orig_b_index + 1;

        if (b_size <= N) {
            score_array[b_size].best = MININT;
            score_array[b_size].best_gap = MININT;
            b_size++;
        }
    }

    /* Pick the optimal path through the now complete
       edit_script_t[tid][][]. This is equivalent to flattening 
       the 2-D array into a 1-D list of actions. */

    a_index = *a_offset;
    b_index = *b_offset;


    //fprintf(stderr, "a_index: %d b_index: %d\n", a_index, b_index);
    script = SCRIPT_SUB;

    if (fence_hit && *fence_hit)
        goto done;

    while (a_index > 0 || b_index > 0) {
        /* Retrieve the next action to perform. Rows of
           the traceback array do not necessarily start
           at offset zero of B, so a correction is needed
           to point to the correct position */

        next_script = 
            edit_script_t[tid][a_index][b_index - edit_start_offset_t[tid][a_index]];

        switch(script) {
            case SCRIPT_GAP_IN_A:
                script = next_script & SCRIPT_OP_MASK;
                if (next_script & SCRIPT_EXTEND_GAP_A)
                    script = SCRIPT_GAP_IN_A;
                break;

            case SCRIPT_GAP_IN_B:
                script = next_script & SCRIPT_OP_MASK;
                if (next_script & SCRIPT_EXTEND_GAP_B)
                    script = SCRIPT_GAP_IN_B;
                break;

            default:
                script = next_script & SCRIPT_OP_MASK;
                break;
        }

        if (script == SCRIPT_GAP_IN_A) {
            b_index--;
        }
        else if (script == SCRIPT_GAP_IN_B) {
            a_index--;
        }
        else {
            a_index--;
            b_index--;
        }
        edit_block->total_num_ops++;
        GapPrelimEditBlockAdd(edit_block, (EGapAlignOpType)script, 1);
    }

done:
    //free(edit_start_offset_t[tid]);
    //free(edit_script_t[tid]);
    return best_score;
}

/* see gapinfo.h for description */
    GapEditScript* 
GapEditScriptNew(Int4 size)

{
    GapEditScript* new;

    if (size <= 0) 
        return NULL;

    new = (GapEditScript*) calloc(1, sizeof(GapEditScript));
    if (new)
    {
        new->size = size;
        new->op_type = (EGapAlignOpType*) calloc(size, sizeof(EGapAlignOpType));
        new->num = (Int4*) calloc(size, sizeof(Int4));
    }
    return new;
}


    GapEditScript*
Blast_PrelimEditBlockToGapEditScript (GapPrelimEditBlock* rev_prelim_tback,
        GapPrelimEditBlock* fwd_prelim_tback)
{
    Boolean merge_ops = FALSE;
    GapEditScript* esp;
    GapPrelimEditScript *op;
    Int4 i;
    Int4 index=0;
    Int4 size = 0;

    if (rev_prelim_tback == NULL || fwd_prelim_tback == NULL)
        return NULL;

    /* The fwd_prelim_tback script will get reversed here as the traceback started from the highest scoring point
       and worked backwards. The rev_prelim_tback script does NOT get reversed.  Since it was reversed when the 
       traceback was produced it's already "forward" */

    if (fwd_prelim_tback->num_ops > 0 && rev_prelim_tback->num_ops > 0 &&
            fwd_prelim_tback->edit_ops[(fwd_prelim_tback->num_ops)-1].op_type == 
            rev_prelim_tback->edit_ops[(rev_prelim_tback->num_ops)-1].op_type)
        merge_ops = TRUE;

    size = fwd_prelim_tback->num_ops+rev_prelim_tback->num_ops;
    if (merge_ops)
        size--;

    esp = GapEditScriptNew(size);

    index = 0;
    for (i=0; i < rev_prelim_tback->num_ops; i++) {
        op = rev_prelim_tback->edit_ops + i;
        esp->op_type[index] = op->op_type;
        esp->num[index] = op->num;
        index++;
    }

    if (fwd_prelim_tback->num_ops == 0)
        return esp;

    if (merge_ops)
        esp->num[index-1] += fwd_prelim_tback->edit_ops[(fwd_prelim_tback->num_ops)-1].num;

    /* If we merge, then we skip the first one. */
    if (merge_ops)
        i = fwd_prelim_tback->num_ops - 2;
    else
        i = fwd_prelim_tback->num_ops - 1;

    for (; i >= 0; i--) {
        op = fwd_prelim_tback->edit_ops + i;
        esp->op_type[index] = op->op_type;
        esp->num[index] = op->num;
        index++;
    }

    return esp;
}


Int2 BLAST_GappedAlignmentWithTraceback( 
        int tid,
        const Uint1* query, const Uint1* subject, BlastGapAlignStruct* gap_align, 
        const BlastScoringParameters* score_params,
        Int4 q_start, Int4 s_start, Int4 query_length, Int4 subject_length,
        Boolean * fence_hit, Int4 **matrix )
{
    Boolean found_start, found_end;
    Int4 score_right, score_left, private_q_length, private_s_length;
    Int4 q_length, s_length;
    //Boolean is_ooframe = score_params->options->is_ooframe;
    Int2 status = 0;
    Boolean switch_seq = FALSE;
    GapPrelimEditBlock *fwd_prelim_tback;
    GapPrelimEditBlock *rev_prelim_tback;

    if (gap_align == NULL)
        return -1;

    fwd_prelim_tback = gap_align->fwd_prelim_tback;
    rev_prelim_tback = gap_align->rev_prelim_tback;
    GapPrelimEditBlockReset(fwd_prelim_tback);
    GapPrelimEditBlockReset(rev_prelim_tback);

    found_start = FALSE;
    found_end = FALSE;

    q_length = query_length;
    s_length = subject_length;

    score_left = 0;
    found_start = TRUE;

    //fprintf(stderr, "q_start: %d s_start: %d gap_open: %d gap_extend: %d dropoff: %d\n", q_start, s_start, score_params->gap_open, score_params->gap_extend, gap_align->gap_x_dropoff);
    score_left = ALIGN_EX(tid, query, subject, q_start+1, s_start+1, 
            &private_q_length, &private_s_length, rev_prelim_tback,
            gap_align, 
            score_params, q_start, FALSE /*reversed*/, TRUE /*reverse_sequence*/,
            matrix,
            fence_hit);


    gap_align->query_start = q_start - private_q_length + 1;
    gap_align->subject_start = s_start - private_s_length + 1;


    //fprintf(stderr, "score_left: %d query_start: %d subject_start: %d\n", score_left, gap_align->query_start, gap_align->subject_start);

    score_right = 0;


    found_end = TRUE;
    score_right = 
        ALIGN_EX(tid, query+q_start, subject+s_start, 
                q_length-q_start-1, s_length-s_start-1, &private_q_length, 
                &private_s_length, fwd_prelim_tback, gap_align, 
                score_params, q_start, FALSE, FALSE,
                matrix,
                fence_hit);

    gap_align->query_stop = q_start + private_q_length + 1;
    gap_align->subject_stop = s_start + private_s_length + 1;


    //fprintf(stderr, "score_right: %d query_stop: %d subject_stop: %d\n", score_right, gap_align->query_stop, gap_align->subject_stop);

    if (found_start == FALSE) {	/* Start never found */
        gap_align->query_start = q_start;
        gap_align->subject_start = s_start;
    }

    if (found_end == FALSE) {
        gap_align->query_stop = q_start - 1;
        gap_align->subject_stop = s_start - 1;
    }

    gap_align->score = score_right + score_left;
    return status;
}


