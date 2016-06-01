#include "blast.h"
#include "stdint.h"


/** 64-bit integers */
#ifndef NCBI_CONST_INT8 /* C Toolkit */
#  ifdef INT64_C /* stdint.h should have this */
#    define NCBI_CONST_INT8(v)   INT64_C(v)
#    define NCBI_CONST_UINT8(v)  UINT64_C(v)
#  elif defined(_MSC_VER)
#    define NCBI_CONST_INT8(v)   v##i64
#    define NCBI_CONST_UINT8(v)  v##ui64
#  else /* Try treating as (unsigned) long long */
#    define NCBI_CONST_INT8(v)   v##LL
#    define NCBI_CONST_UINT8(v)  v##ULL
#  endif
#endif
/**
 * Do a simple gapped extension to the right from the beginning of query and
 * subject ranges examining only matches and mismatches. The extension stops
 * when there are more than max_shift mismatches or mismatches or gaps are not
 * followed by two identical matches. This is a simplified version of the
 * Danielle and Jean Thierry-Miegs' jumper
 * alignment implemented in NCBI Magic
 * http://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/Download/Downloads.html
 *
 * @param query_seq Query sequence [in]
 * @param query_len Query length [in]
 * @param subject_seq Subject sequence [in]
 * @param subject_len Subject length [in]
 * @param max_shift Maximum number of mismatches or gaps, extension stops if
 *        this number is reached [in]
 * @param query_ext_len Extension length on the query [out]
 * @param subject_ext_len Extension length on the subject [out]
 * @param align_len Alignment length [out]
 * @return Number of identical residues
 */
static int s_ExtendRight(Uint1* query_seq, int query_len,
                         Uint1* subject_seq, int subject_len,
                         int max_shift,
                         int* query_ext_len, int* subject_ext_len,
                         int* align_len)
{
    int num_identical = 0;
    int q_pos, s_pos;
    int gaps_in_query = 0;
    int gaps_in_subject = 0;
    q_pos = 0;
    s_pos = 0;
    while (q_pos < query_len && s_pos < subject_len) {
        int n;
        int match = 0;

        while (q_pos < query_len && s_pos < subject_len
               && query_seq[q_pos] == subject_seq[s_pos]) {

            num_identical++;
            q_pos++;
            s_pos++;
        }

        /* try to skip mismatches or gaps */
        for (n=1; n < max_shift && q_pos + n + 1 < query_len
                 && s_pos + n + 1 < subject_len && !match; n++) {

            /* mismatches */
            if (query_seq[q_pos + n] == subject_seq[s_pos + n]
                && query_seq[q_pos + n + 1] == subject_seq[s_pos + n + 1]) {
                
                /* we have already checked that two positions behind mismatches
                   match so we can advance further */
                q_pos += n + 2;
                s_pos += n + 2;
                num_identical += 2;
                match = 1;
            }

            /* gap in subject */
            if (!match && query_seq[q_pos + n] == subject_seq[s_pos]
                && query_seq[q_pos + n + 1] == subject_seq[s_pos + 1]) {
                
                q_pos += n + 2;
                s_pos += 2;
                num_identical += 2;
                gaps_in_subject += n;
                match = 1;
            }

            /* gap in query */
            if (!match && query_seq[q_pos] == subject_seq[s_pos + n]
                && query_seq[q_pos + 1] == subject_seq[s_pos + n + 1]) {

                q_pos += 2;
                s_pos += n + 2;
                num_identical += 2;
                gaps_in_query += n;
                match = 1;
            }
        }

        if (match) {
            continue;
        }

        /* exit the loop */
        break;
    }
    *query_ext_len = q_pos;
    *subject_ext_len = s_pos;
    *align_len = q_pos > s_pos ? q_pos + gaps_in_query : s_pos + gaps_in_subject;

    return num_identical;
}


/** 
 * Extend left from the end of the sequence and subject ranges and count
 * identities. The extension stops when there are more than max_shift
 * mismatches or mismatches or gaps are not followed by two identical matches.
 * See description for s_ExtendRight for more details.
 *
 * @param query_seq Query sequence [in]
 * @param query_len Query length [in]
 * @param subject_seq Subject sequence [in]
 * @param subject_len Subject length [in]
 * @param max_shift Maximum number of mismatches or gaps, extension stops if
 *        this number is reached [in]
 * @param query_ext_len Extension length on the query [out]
 * @param subject_ext_len Extension length on the subject [out]
 * @param align_len Alignment length [out]
 * @return Number of identical residues
 */
static int s_ExtendLeft(Uint1* query_seq, int query_len,
                        Uint1* subject_seq, int subject_len,
                        int max_shift,
                        int* query_ext_len, int* subject_ext_len,
                        int* align_len)
{
    int q_pos = query_len - 1;
    int s_pos = subject_len - 1;
    int num_identical = 0;
    int gaps_in_query = 0;
    int gaps_in_subject = 0;
    while (q_pos >= 0 && s_pos >= 0) {
        int n;
        int match = 0;

        /* process identies */
        while (q_pos > 0 && s_pos > 0 && query_seq[q_pos] == subject_seq[s_pos]) {
            num_identical++;
            q_pos--;
            s_pos--;
        }

        /* try to skip mismatches or gaps */
        for (n=1;n < max_shift && q_pos - n - 1 > 0 && s_pos - n - 1 > 0
                 && !match; n++) {

            /* mismatch */
            if (query_seq[q_pos - n] == subject_seq[s_pos - n]
                && query_seq[q_pos - n - 1] == subject_seq[s_pos - n - 1]) {
                q_pos -= n + 2;
                s_pos -= n + 2;
                num_identical += 2;
                match = 1;
            }

            /* gap in subject */
            if (!match && query_seq[q_pos - n] == subject_seq[s_pos]
                && query_seq[q_pos - n - 1] == subject_seq[s_pos - 1]) {
                q_pos -= n + 2;
                s_pos -= 2;
                num_identical += 2;
                gaps_in_subject += n;
                match = 1;
            }

            /* gap in query */
            if (!match && query_seq[q_pos] == subject_seq[s_pos - n]
                && query_seq[q_pos - 1] == subject_seq[s_pos - n - 1]) {
                q_pos -= 2;
                s_pos -= n + 2;
                num_identical += 2;
                gaps_in_query += n;
                match = 1;
            }
        }

        if (match) {
            continue;
        }

        break;
    }
    *query_ext_len = query_len - q_pos - 1;
    *subject_ext_len = subject_len - s_pos - 1;
    *align_len += *query_ext_len > *subject_ext_len ?
        *query_ext_len + gaps_in_query : *subject_ext_len + gaps_in_subject;

    return num_identical;
}


/** 
 * Get hash for a word of word_size residues assuming 28-letter alphabet
 *
 * @param data Sequence [in]
 * @param word_size Word size [in]
 * @return Hash value
 */
static Uint8 s_GetHash(const Uint1* data, int word_size)
{
    Uint8 hash = 0;
    int k;
    for (k=0;k < word_size;k++) {
        hash <<= 5;
        hash += (Int8)data[k];
    }
    return hash;
}

/** 
 * Find a local number of identical residues in two aligned sequences by
 * finding word matches and doing a simple gapped extensions from the word hits
 *
 * @param query_seq Query sequence [in]
 * @param query_hashes Array of query words with index of each word
 *        corresponding to word position in the query [in]
 * @param query_len Query length [in]
 * @param subject_seq Subject sequence [in]
 * @param subject_len Subject length [in]
 * @param max_shift Maximum number of local mismatches or gaps for extensions
 *        [in]
 * @return Number of identical residues
 */
static int s_FindNumIdentical(Uint1* query_seq,
                              const Uint8* query_hashes,
                              int query_len,
                              Uint1* subject_seq,
                              int subject_len,
                              int max_shift)
{
    int word_size = 8;         /* word size for k-mer matching */
    Uint8 hash = 0;
    Uint8 mask = NCBI_CONST_UINT8(0xFFFFFFFFFF); /* mask for computing hash
                                                    values */
    int query_from = 0;
    int subject_from = 0;

    int s_pos;                 /* position in the subject sequence */
    int num_identical = 0;     /* number of identical residues found */
    Boolean match = FALSE;

    /* if query or subject length is smaller than word size, exit */
    if (!query_seq || !query_hashes || !subject_seq
        || query_len < word_size || subject_len < word_size) {

        return 0;
    }

    /* for each subject position */
    for (s_pos = 0; s_pos < subject_len - word_size; s_pos++) {
        int q_pos;

        /* find word hash */
        if (s_pos == 0 || match) {
            hash = s_GetHash(&subject_seq[s_pos], word_size);
        }
        else {
            hash <<= 5;
            hash &= mask;
            hash += subject_seq[s_pos + word_size - 1];
        }

        /* find matching query word; index of hash is position of the word 
           the query */
        for (q_pos = query_from;q_pos < query_len - word_size; q_pos++) {
            if (query_hashes[q_pos] == hash) {
                break;
            }
        }

        /* if match */
        if (q_pos < query_len - word_size) {
            int query_start = q_pos;
            int subject_start = s_pos;

            int query_left_len, query_right_len;
            int subject_left_len, subject_right_len;
            int align_len_left=0, align_len_right=0;
            
            match = TRUE;
            num_identical += word_size;

            /* extend left from word match */
            num_identical += s_ExtendLeft(query_seq + query_start - 1,
                                          query_start - query_from,
                                          subject_seq + subject_start - 1,
                                          subject_start - subject_from,
                                          max_shift,
                                          &query_left_len, &subject_left_len,
                                          &align_len_left);

            /* extend right from word match */
            num_identical += s_ExtendRight(query_seq + query_start + word_size,
                                       query_len - query_start - word_size,
                                       subject_seq + subject_start + word_size,
                                       subject_len - subject_start - word_size,
                                       max_shift,
                                       &query_right_len, &subject_right_len,
                                       &align_len_right);


            /* disregard already matched and extended words when matching
               further positions */

            query_from = query_start + word_size + query_right_len;
            subject_from = subject_start + word_size + subject_right_len;
            /* s_pos will be incremented in the loop */
            s_pos = subject_from - 1;
        }
        else {
            match = FALSE;
        }
    }

    return num_identical;
}

/**
 * Test whether the aligned parts of two sequences that
 * have a high-scoring gapless alignment are nearly identical.
 *
 * First extend from the left end of the query and subject ranges and stop if
 * there are too manu mismatches. Then extend from the right end. Then for the
 * remaining protion of ths sequences find matching words and extend left and
 * right from the word hit. Repeat the last steo until the whole alignment
 * ranges are processed.
 *
 * @params seqData Subject sequence [in]
 * @params seqOffse Starting offset of the subject sequence in alignment data
 *        [in]
 * @params queryData Query sequence [in]
 * @params queryOffset Starting offset of the query sequence in alignment data
 *         [in]
 * @param query_words Array of query words with word index corresponding to
 *        word's position in the query [in]
 * @param align Alignment data [in]
 * @return True if sequence parts are nearly identical, false otherwise
 */
Boolean
s_TestNearIdentical(const BlastCompo_SequenceData* seqData,
                    const int seqOffset,
                    const char* queryData,
                    const int queryOffset,
                    const Uint8* query_words,
                    const struct ungappedExtension* align)
{
    int qStart = align->start.queryOffset;
    /* align->queryEnd points to one position past alignment end */
    int qEnd = align->end.queryOffset;
    int sStart = align->start.subjectOffset;
    int sEnd = align->end.subjectOffset;
    const double kMinFractionNearIdentical = 0.96;
    int max_shift = 8;

    int query_len = qEnd - qStart + 1;
    int subject_len = sEnd - sStart + 1;
    int align_len = MIN(query_len, subject_len);

    int query_left_len = 0;
    int subject_left_len = 0;
    int query_right_len = 0;
    int subject_right_len = 0;
    int align_left_len = 0;
    int align_right_len = 0;

    double fraction_identical;

    /* first find number of identies going from the beginning of the query
       and subject ranges */
    int num_identical = s_ExtendRight(queryData + qStart, query_len,
                                      seqData->data + sStart, subject_len,
                                      max_shift,
                                      &query_right_len, &subject_right_len,
                                      &align_right_len);

    /* if the whole query range was processed return near identical status */
    if (query_right_len >= query_len || subject_right_len >= subject_len) {
        fraction_identical = (double)num_identical / (double)align_len;
        ASSERT(fraction_identical - 1.0 < 1e-10);
        return fraction_identical > kMinFractionNearIdentical;
    }

    /* find the number of identies going from the end of the query and subject
       ranges */
    num_identical += s_ExtendLeft(queryData + qStart + query_right_len,
                                  query_len - query_right_len,
                                  seqData->data + sStart + subject_right_len,
                                  subject_len - subject_right_len,
                                  max_shift,
                                  &query_left_len, &subject_left_len,
                                  &align_left_len);

    /* if the whole alignment ranges where covered, return the near identical
       status */
    if (query_left_len + query_right_len >= query_len
        || subject_left_len + subject_right_len >= subject_len) {

        fraction_identical = (double)num_identical / (double)(align_len);
        ASSERT(fraction_identical - 1.0 < 1e-10);
        return fraction_identical > kMinFractionNearIdentical;
    }

    /* find the number of identical matches in the middle portion of the
       alignment ranges */
    num_identical += s_FindNumIdentical(queryData + qStart + query_right_len,
                            query_words + qStart + query_right_len,
                            query_len - query_left_len - query_right_len,
                            seqData->data + sStart + subject_right_len,
                            subject_len - subject_left_len - subject_right_len,
                            max_shift);

    fraction_identical = (double)num_identical / (double)align_len;
    ASSERT(fraction_identical - 1.0 < 1e-10);
    if (fraction_identical > kMinFractionNearIdentical) {
        return TRUE;
    }
    else {
        return FALSE;
    }
}
