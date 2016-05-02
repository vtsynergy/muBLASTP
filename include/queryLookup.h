#ifndef _queryLookup_
#define _queryLookup_

/*#define QUERY_BLOCK_SIZE 100*/
#define NUM_QUERY_BLOCK 50

#define PV_ARRAY_TYPE uint4     /**< The pv_array 'native' type. */
#define PV_ARRAY_BYTES 4        /**< number of BYTES in 'native' type. */
#define PV_ARRAY_BTS 5          /**< bits-to-shift from lookup_index to pv_array index. */
#define PV_ARRAY_MASK 31        /**< amount to mask off. */

/** Set the bit at position 'index' in the PV 
 *  array bitfield within 'lookup'
 */
#define PV_SET(lookup, index, shift)    \
    lookup[(index) >> (shift)] |= (PV_ARRAY_TYPE)1 << ((index) & PV_ARRAY_MASK)

/** Test the bit at position 'index' in the PV 
 *  array bitfield within 'lookup'
 */
#define PV_TEST(lookup, index, shift)                   \
      ( lookup[(index) >> (shift)] &                    \
        ((PV_ARRAY_TYPE)1 << ((index) & PV_ARRAY_MASK)) )

extern PV_ARRAY_TYPE *pv;
extern PV_ARRAY_TYPE *pv_arr[];

struct initialWord_protein_query
{
    uint4  numQueryPositions;
    uint4  allocQueryPositions;
    /*int2 wordCode;*/
    /*uint4 queryPositionOffset;*/
    uint32_t* querySequencePositions;
    uint32_t embeddedPositions[3];
};

extern struct initialWord_protein_query *proteinLookup_query;

extern struct initialWord_protein_query *proteinLookup_query_blk[];

extern int numQueryBlk;

void proteinLookup_query_initial(int4 numCodes, int wordLength);

void proteinLookup_query_initial_blk(int4 numCodes, int wordLength, int block_id);

void queryIndex_free();

void queryIndex_free_blk();

void proteinLookup_query_add(uint4 sequenceNum, unsigned char *sequence,
        int4 subjectLength, int4 wordLength,
        struct scoreMatrix scoreMatrix);

void proteinLookup_query_add_blk(uint4 sequenceNum, unsigned char *sequence,
        int4 subjectLength, int4 wordLength,
        struct scoreMatrix scoreMatrix, int block_id);

void proteinLookup_sort();

void proteinLookup_sort_blk(int block_id);
#endif
