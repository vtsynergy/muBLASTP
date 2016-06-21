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

//#ifndef INT4_MIN
///** Smallest (most negative) number represented by signed int */
//#define INT4_MIN    (-2147483647-1)
//#endif

#define MININT INT4_MIN/2


//#define INT4_MIN    (-2147483647-1)
/* * Lower bound for scores. Divide by two to prevent underflows. */

#define HSP_MAX_WINDOW 11

    Int4 
BlastGetStartForGappedAlignment (const Uint1* query, const Uint1* subject,
        int2 **matrix, Uint4 q_start, Uint4 q_length, 
        Uint4 s_start, Uint4 s_length)
{
    Int4 index1, max_offset, score, max_score, hsp_end;
    const Uint1* query_var,* subject_var;
    //Boolean positionBased = (sbp->psi_matrix != NULL);

    if (q_length <= HSP_MAX_WINDOW) {
        max_offset = q_start + q_length/2;
        return max_offset;
    }

    hsp_end = q_start + HSP_MAX_WINDOW;
    query_var = query + q_start;
    subject_var = subject + s_start;
    score=0;
    for (index1=q_start; index1<hsp_end; index1++) {
        //if (!(positionBased))
        score += matrix[*query_var][*subject_var];
        //else
        //score += sbp->psi_matrix->pssm->data[index1][*subject_var];
        query_var++; subject_var++;
    }
    max_score = score;
    max_offset = hsp_end - 1;
    hsp_end = q_start + MIN(q_length, s_length);
    for (index1=q_start + HSP_MAX_WINDOW; index1<hsp_end; index1++) {
        //if (!(positionBased)) {
        score -= matrix[*(query_var-HSP_MAX_WINDOW)][*(subject_var-HSP_MAX_WINDOW)];
        score += matrix[*query_var][*subject_var];
        //} else {
        //score -= sbp->psi_matrix->pssm->data[index1-HSP_MAX_WINDOW][*(subject_var-HSP_MAX_WINDOW)];
        //score += sbp->psi_matrix->pssm->data[index1][*subject_var];
        //}
        if (score > max_score) {
            max_score = score;
            max_offset = index1;
        }
        query_var++; subject_var++;
    }
    if (max_score > 0)
        max_offset -= HSP_MAX_WINDOW/2;
    else 
        max_offset = q_start;

    return max_offset;
}

    Int4 
Blast_SemiGappedAlign(const Uint1* A, const Uint1* B, Int4 M, Int4 N,
        Int4* a_offset, Int4* b_offset, 
        BlastGapAlignStruct* gap_align, 
        const BlastScoringParameters* score_params, 
        Int4 query_offset, Boolean reversed, Boolean reverse_sequence,
        Int2** matrix)
{
    Int4 i;                     /* sequence pointers and indices */
    Int4 a_index;
    Int4 b_index, b_size, first_b_index, last_b_index, b_increment;
    const Uint1* b_ptr;
    const Uint1* a_ptr;

    BlastGapDP* score_array;

    Int4 gap_open;              /* alignment penalty variables */
    Int4 gap_extend;
    Int4 gap_open_extend;
    Int4 x_dropoff;

    Int2* matrix_row = NULL;

    Int4 score;                 /* score tracking variables */
    Int4 score_gap_row;
    Int4 score_gap_col;
    Int4 next_score;
    Int4 best_score;
    Int4 num_extra_cells;

    /* do initialization and sanity-checking */
    *a_offset = 0;
    *b_offset = 0;
    gap_open = score_params->gap_open;
    gap_extend = score_params->gap_extend;
    gap_open_extend = gap_open + gap_extend;
    x_dropoff = gap_align->gap_x_dropoff;

    //fprintf(stderr, "gap_open: %d gap_extend: %d x_dropoff: %d M: %d N: %d\n", gap_open, gap_extend, x_dropoff, M, N);

    if (x_dropoff < gap_open_extend)
        x_dropoff = gap_open_extend;

    if(N <= 0 || M <= 0) 
        return 0;

    /* Allocate and fill in the auxiliary bookeeping structures.
       Since A and B could be very large, maintain a window
       of auxiliary structures only large enough to contain to current
       set of DP computations. The initial window size is determined
       by the number of cells needed to fail the x-dropoff test */

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

    score_array = gap_align->dp_mem;
    score = -gap_open_extend;
    score_array[0].best = 0;
    score_array[0].best_gap = -gap_open_extend;

    for (i = 1; i <= N; i++) {
        if (score < -x_dropoff) 
            break;

        score_array[i].best = score;
        score_array[i].best_gap = score - gap_open_extend; 
        score -= gap_extend;
    }

    /* The inner loop below examines letters of B from 
       index 'first_b_index' to 'b_size' */

    b_size = i;
    best_score = 0;
    first_b_index = 0;
    if (reverse_sequence)
        b_increment = -1;
    else
        b_increment = 1;

    for (a_index = 1; a_index <= M; a_index++) {
        /* pick out the row of the score matrix 
           appropriate for A[a_index] */

        if(reverse_sequence)
        {
            matrix_row = matrix[ A[ M - a_index ] ];
            a_ptr = &A[ M - a_index ];
        }
        else
        {
            matrix_row = matrix[ A[ a_index ] ];
            a_ptr = &A[ a_index ];
        }

        if(reverse_sequence)
            b_ptr = &B[N - first_b_index];
        else
            b_ptr = &B[first_b_index];

        /* initialize running-score variables */
        score = MININT;
        score_gap_row = MININT;
        last_b_index = first_b_index;

        for (b_index = first_b_index; b_index < b_size; b_index++) {

            b_ptr += b_increment;
            score_gap_col = score_array[b_index].best_gap;
            next_score = score_array[b_index].best + matrix_row[ *b_ptr ];

            //fprintf(stderr, "next_score: %d b_index: %d score_array[b_index].best: %d *b_ptr: %c matrix_row[ *b_ptr  ]: %d a_index: %d *a_ptr: %c\n", next_score, b_index, score_array[b_index].best, encoding_getLetter(*b_ptr), matrix_row[ *b_ptr  ], a_index, encoding_getLetter(*a_ptr));

            if (score < score_gap_col)
                score = score_gap_col;

            if (score < score_gap_row)
                score = score_gap_row;

            if (best_score - score > x_dropoff) {

                /* the current best score failed the X-dropoff
                   criterion. Note that this does not stop the
                   inner loop, only forces future iterations to
                   skip this column of B. 

                   Also, if the very first letter of B that was
                   tested failed the X dropoff criterion, make
                   sure future inner loops start one letter to 
                   the right */

                if (b_index == first_b_index)
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
                    //fprintf(stderr, "score: %d best_score: %d a_offset: %d b_offset: %d\n", score, best_score, a_index, b_index);
                }

                /* If starting a gap at this position will improve
                   the best row, or column, score, update them to 
                   reflect that. */

                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;
                score_array[b_index].best_gap = MAX(score - gap_open_extend,
                        score_gap_col);
                score_gap_row = MAX(score - gap_open_extend, score_gap_row);
                score_array[b_index].best = score;
            }

            score = next_score;
        }

        /* Finish aligning if the best scores for all positions
           of B will fail the X-dropoff test, i.e. the inner loop 
           bounds have converged to each other */

        if (first_b_index == b_size)
            break;

        /* enlarge the window for score data if necessary */

        if (last_b_index + num_extra_cells + 3 >= gap_align->dp_mem_alloc) {

#if 0
            gap_align->dp_mem_alloc = MAX(last_b_index + num_extra_cells + 100,
                    2 * gap_align->dp_mem_alloc);
            score_array = (BlastGapDP *)realloc(score_array,
                    gap_align->dp_mem_alloc *
                    sizeof(BlastGapDP));
            gap_align->dp_mem = score_array;
#endif
            fprintf(stderr, "DP_MEM overflow! required = %d, allocated = %d\n", last_b_index + num_extra_cells + 3, gap_align->dp_mem_alloc);
            exit(0);
        }

        if (last_b_index < b_size - 1) {
            /* This row failed the X-dropoff test earlier than
               the last row did; just shorten the loop bounds
               before doing the next row */

            b_size = last_b_index + 1;
        }
        else {
            /* The inner loop finished without failing the X-dropoff
               test; initialize extra bookkeeping structures until
               the X dropoff test fails or we run out of letters in B. 
               The next inner loop will have larger bounds */

            while (score_gap_row >= (best_score - x_dropoff) && b_size <= N) {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                b_size++;
            }
        }

        if (b_size <= N) {
            score_array[b_size].best = MININT;
            score_array[b_size].best_gap = MININT;
            b_size++;
        }
    }

    return best_score;
}


#define RESTRICT_SIZE 10

    Int4 
s_RestrictedGappedAlign(const Uint1* A, const Uint1* B, Int4 M, Int4 N,
        Int4* a_offset, Int4* b_offset,
        BlastGapAlignStruct* gap_align, 
        const BlastScoringParameters* score_params, 
        Int4 query_offset, Boolean reverse_sequence, Int2 **matrix)
{
    /* see Blast_SemiGappedAlign for general details */

    Int4 i;
    Int4 a_index;
    Int4 b_index, b_size, first_b_index, last_b_index, b_increment;
    const Uint1* b_ptr;
    Int4 b_gap;

    BlastGapDP* score_array;

    Int4 gap_open;
    Int4 gap_extend;
    Int4 gap_open_extend;
    Int4 x_dropoff;

    Int2* matrix_row = NULL;

    Int4 score;
    Int4 score_gap_row;
    Int4 score_gap_col;
    Int4 next_score;
    Int4 best_score;
    Int4 num_extra_cells;

    //matrix = gap_align->sbp->matrix->data;

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

    score_array = gap_align->dp_mem;
    score = -gap_open_extend;
    score_array[0].best = 0;
    score_array[0].best_gap = -gap_open_extend;

    for (i = 1; i <= N; i++) {
        if (score < -x_dropoff) 
            break;

        score_array[i].best = score;
        score_array[i].best_gap = score - gap_open_extend; 
        score -= gap_extend;
    }

    b_size = i;
    best_score = 0;
    first_b_index = 0;
    b_gap = 0;
    if (reverse_sequence)
        b_increment = -1;
    else
        b_increment = 1;

    for (a_index = 1; a_index <= M; a_index++) {

        if(reverse_sequence)
            matrix_row = matrix[ A[ M - a_index ] ];
        else
            matrix_row = matrix[ A[ a_index ] ];

        if(reverse_sequence)
            b_ptr = &B[N - first_b_index];
        else
            b_ptr = &B[first_b_index];

        score = MININT;
        score_gap_row = MININT;
        last_b_index = first_b_index;

        /* The double loop that computes the alignment is essentially
           unchanged from that of Blast_SemiGappedAlign. The only
           real difference is that a gap in A or B is not allowed to 
           start unless the offset of the gap is divisible by 
           RESTRICT_SIZE. This allows the ordinary dynamic programming
           recurrence relations to be simplified */

        if (a_index % RESTRICT_SIZE != 0) {

            /* a gap will never start in A; do not bother checking
               or updating score_gap_row */

            for (b_index = first_b_index; b_index < b_size; b_index++) {

                b_ptr += b_increment;
                next_score = score_array[b_index].best + matrix_row[ *b_ptr ];

                if (b_index != b_gap) {

                    /* the majority of cases fall here; a gap
                       may not start in either A or B */

                    if (best_score - score > x_dropoff) {
                        score_array[b_index].best = MININT;
                        if (b_index == first_b_index)
                            first_b_index++;
                    }
                    else {
                        last_b_index = b_index;
                        if (score > best_score) {
                            best_score = score;
                            *a_offset = a_index;
                            *b_offset = b_index;
                        }
                        score_array[b_index].best = score;
                    }
                }
                else {

                    /* a gap may start in B. Update b_gap, the
                       offset when this will next happen, and
                       compute the two-term recurrence */

                    b_gap += RESTRICT_SIZE;
                    score_gap_col = score_array[b_index].best_gap;

                    if (score < score_gap_col)
                        score = score_gap_col;

                    if (best_score - score > x_dropoff) {
                        score_array[b_index].best = MININT;
                        if (b_index == first_b_index)
                            first_b_index++;
                    }
                    else {
                        last_b_index = b_index;
                        if (score > best_score) {
                            best_score = score;
                            *a_offset = a_index;
                            *b_offset = b_index;
                        }

                        score_gap_col -= gap_extend;
                        score_array[b_index].best_gap = 
                            MAX(score - gap_open_extend, score_gap_col);
                        score_array[b_index].best = score;
                    }
                }
                score = next_score;
            }

            score_gap_row = score;
        }
        else {

            /* gap may start in A */

            for (b_index = first_b_index; b_index < b_size; b_index++) {

                b_ptr += b_increment;
                next_score = score_array[b_index].best + matrix_row[ *b_ptr ];

                if (b_index != b_gap) {

                    /* gap may not start in B. Compute
                       the resulting two-term recurrence */

                    if (score < score_gap_row)
                        score = score_gap_row;

                    if (best_score - score > x_dropoff) {
                        score_array[b_index].best = MININT;
                        if (b_index == first_b_index)
                            first_b_index++;
                    }
                    else {
                        last_b_index = b_index;
                        if (score > best_score) {
                            best_score = score;
                            *a_offset = a_index;
                            *b_offset = b_index;
                        }

                        score_gap_row -= gap_extend;
                        score_gap_row = MAX(score - gap_open_extend, 
                                score_gap_row);
                        score_array[b_index].best = score;
                    }
                }
                else {

                    /* the ordinary case: a gap may start in A
                       or B and the full three-term recurrence
                       must be computed. This happens once every
                       RESTRICT_SIZE*RESTRICT_SIZE cells */

                    b_gap += RESTRICT_SIZE;
                    score_gap_col = score_array[b_index].best_gap;

                    if (score < score_gap_col)
                        score = score_gap_col;
                    if (score < score_gap_row)
                        score = score_gap_row;

                    if (best_score - score > x_dropoff) {
                        score_array[b_index].best = MININT;
                        if (b_index == first_b_index)
                            first_b_index++;
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
                        score_array[b_index].best_gap = 
                            MAX(score - gap_open_extend, score_gap_col);
                        score_gap_row = MAX(score - gap_open_extend, 
                                score_gap_row);
                        score_array[b_index].best = score;
                    }
                }
                score = next_score;
            }
        }

        if (first_b_index == b_size)
            break;

        /* compute the first offset of the next row where a
           gap in B may start */

        b_index = first_b_index % RESTRICT_SIZE;
        b_gap = first_b_index;
        if (b_index > 0)
            b_gap += RESTRICT_SIZE - b_index;

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
            /* The inner loop finished without failing the X-dropoff
               test; initialize extra bookkeeping structures until
               the X dropoff test fails or we run out of letters in B. 
               The next inner loop will have larger bounds.

               Note that if gaps in A are not allowed for this row, then
               score_gap_row is -infinity and the loop below will not 
               extend the region to explore */

            while (score_gap_row >= (best_score - x_dropoff) && b_size <= N) {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                b_size++;
            }
        }

        if (b_size <= N) {
            score_array[b_size].best = MININT;
            score_array[b_size].best_gap = MININT;
            b_size++;
        }
    }

    return best_score;
}


#define MAX_SUBJECT_OFFSET 90000
#define MAX_TOTAL_GAPS 3000

    void 
AdjustSubjectRange(Int4* subject_offset_ptr, Int4* subject_length_ptr, 
        Int4 query_offset, Int4 query_length, Int4* start_shift)
{
    Int4 s_off;
    Int4 subject_length = *subject_length_ptr;
    Int4 max_extension_left, max_extension_right;

    /* If subject sequence is not too long, leave everything as is */
    if (subject_length < MAX_SUBJECT_OFFSET) {
        *start_shift = 0;
        return;
    }

    s_off = *subject_offset_ptr;
    /* Maximal extension length is the remaining length in the query, plus 
       an estimate of a maximal total number of gaps. */
    max_extension_left = query_offset + MAX_TOTAL_GAPS;
    max_extension_right = query_length - query_offset + MAX_TOTAL_GAPS;

    if (s_off <= max_extension_left) {
        *start_shift = 0;
    } else {
        *start_shift = s_off - max_extension_left;
        *subject_offset_ptr = max_extension_left;
    }

    *subject_length_ptr = 
        MIN(subject_length, s_off + max_extension_right) - *start_shift;
}

Int2 
s_BlastProtGappedAlignment(
        BLAST_SequenceBlk* query_blk, 
        Int4 q_off,
        BLAST_SequenceBlk* subject_blk, 
        Int4 s_off,
        BlastGapAlignStruct* gap_align,
        const BlastScoringParameters* score_params,
        Boolean restricted_alignment,
        struct ungappedExtension *ungappedExtension,
        struct scoreMatrix scoreMatrix)
{
    Boolean found_start, found_end;
    Int4 q_length=0, s_length=0, score_right, score_left;
    Int4 private_q_start, private_s_start;
    Uint1* query=NULL,* subject=NULL;
    Boolean switch_seq = FALSE;
    Int4 query_length = query_blk->length;
    Int4 subject_length = subject_blk->length;
    Int4 subject_shift = 0;

    Int2 **matrix = scoreMatrix.matrix;

    q_length = q_off + 1;
    s_length = s_off + 1;
    query = query_blk->sequence;
    subject = subject_blk->sequence;

    AdjustSubjectRange(&s_length, &subject_length, q_length, query_length, 
            &subject_shift);

    //fprintf(stderr, "q_off: %d s_off: %d query_length: %d q_length: %d subject_length: %d s_length: %d gap_x_dropoff: %d gap_open: %d gap_extend: %d subject_shift: %d\n ", 
    //q_off, s_off, query_length, q_length, subject_length, s_length, gap_align->gap_x_dropoff, score_params->gap_open, score_params->gap_extend, subject_shift);
    //fprintf(stderr, "q_length: %d s_length: %d query_length: %d subject_length: %d subject_shift: %d\n", q_length, s_length, query_length, subject_length, subject_shift);

    found_start = FALSE;
    found_end = FALSE;

    /* Looking for "left" score */
    score_left = 0;
    if (q_length != 0 && s_length != 0) {
        found_start = TRUE;
        if (restricted_alignment) {
            score_left = s_RestrictedGappedAlign(query, subject+subject_shift, 
                    q_length, s_length,
                    &private_q_start, &private_s_start,
                    gap_align, score_params, 
                    q_off, 
                    TRUE, 
                    matrix);
        }
        else {
            score_left = Blast_SemiGappedAlign(query, subject+subject_shift, 
                    q_length, s_length,
                    &private_q_start, &private_s_start,  
                    gap_align, score_params, 
                    q_off, 
                    FALSE, 
                    TRUE,
                    matrix);
        }

        gap_align->query_start = q_length - private_q_start;
        gap_align->subject_start = s_length - private_s_start + subject_shift;
    }

    score_right = 0;
    if (q_length < query_length && s_length < subject_length) {
        found_end = TRUE;
        if (restricted_alignment) {
            score_right = s_RestrictedGappedAlign(
                    query + q_off, 
                    subject + s_off, 
                    query_length-q_length, 
                    subject_length-s_length, 
                    &(gap_align->query_stop), 
                    &(gap_align->subject_stop), 
                    gap_align, score_params, 
                    q_off, 
                    FALSE, 
                    matrix);
        }
        else {
            score_right = Blast_SemiGappedAlign(
                    query + q_off, 
                    subject + s_off, 
                    query_length-q_length, 
                    subject_length-s_length, 
                    &(gap_align->query_stop), 
                    &(gap_align->subject_stop), 
                    gap_align, score_params, 
                    q_off, 
                    FALSE, 
                    FALSE,
                    matrix);
        }
        /* Make end offsets point to the byte after the end of the 
           alignment */
        gap_align->query_stop += q_off + 1;
        gap_align->subject_stop += s_off + 1;
    }


    if (found_start == FALSE) {  /* impossible for in-frame */
        gap_align->query_start = q_off;
        gap_align->subject_start = s_off;
    }
    if (found_end == FALSE) {    /* impossible for out-of-frame */
        gap_align->query_stop = q_off;
        gap_align->subject_stop = s_off;
    }

    gap_align->score = score_right+score_left;
    //fprintf(stderr, "%d %d %d\n", gap_align->score, score_right, score_left);

    return 0;
}

BlastHSP* Blast_HSPNew(void)
{
     BlastHSP* new_hsp = (BlastHSP*) calloc(1, sizeof(BlastHSP));
     return new_hsp;
}

BlastHSP* Blast_HSPInit(BlastHSP *Blast_HSP_arr, int *Blast_HSP_cnt)
{
    BlastHSP* new_hsp = Blast_HSP_arr + *Blast_HSP_cnt;
    (*Blast_HSP_cnt)++;
    if(*Blast_HSP_cnt > MAX_NUM_HSP)
    {
        fprintf(stderr, "HSP_array overflow\n");
    }
    return new_hsp;
}
