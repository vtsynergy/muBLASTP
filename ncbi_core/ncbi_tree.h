
/** How HSPs added to an interval tree are indexed */
typedef enum EITreeIndexMethod {
    eQueryOnly,              /**< Index by query offset only */
    eQueryAndSubject,        /**< Index by query and then by subject offset */
    eQueryOnlyStrandIndifferent /**< Index by query offset only.  Also do not
                                     distinguish between query strands for 
                                     region definition. -RMH- */
} EITreeIndexMethod;



typedef struct SIntervalNode {
    Int4 leftend;   /**< The left endpoint of the region this node describes */
    Int4 rightend;  /**< The right endpoint of the region this node describes */
    Int4 leftptr;   /**< Offset to the subtree describing the left half
                         of the region, OR the query start offset (leaf 
                         nodes only) */
    Int4 midptr;    /**< Used for linked list of segments that cross the
                         center of the region */
    Int4 rightptr;  /**< Offset to the subtree describing the right half
                         of the region */
    BlastHSP *hsp;  /**< The HSP contained in this region (only non-NULL
                         for leaf nodes) */
} SIntervalNode;

/** Main structure describing an interval tree. */
typedef struct BlastIntervalTree {
    SIntervalNode *nodes;    /**< Pool of tree nodes to allocate from */
    Int4 num_alloc;          /**< Number of nodes allocated */
    Int4 num_used;           /**< Number of nodes actually in use */
    Int4 s_min;              /**< minimum subject offset possible */
    Int4 s_max;              /**< maximum subject offset possible */
} BlastIntervalTree;

BlastIntervalTree* 
Blast_IntervalTreeInit(Int4 q_start, Int4 q_end,
                       Int4 s_start, Int4 s_end);

Boolean
BlastIntervalTreeContainsHSP(const BlastIntervalTree *tree, 
        const BlastHSP *hsp,
        //const BlastQueryInfo *query_info,
        Int4 min_diag_separation);

Int2 
BlastIntervalTreeAddHSP(BlastHSP *hsp, BlastIntervalTree *tree,
        //const BlastQueryInfo *query_info,
                        EITreeIndexMethod index_method);

BlastIntervalTree* 
Blast_IntervalTreeFree(BlastIntervalTree *tree);

BlastIntervalTree* 
Blast_IntervalTreeFree2(BlastIntervalTree *tree);

BlastIntervalTree* 
Blast_IntervalTreeInit2(Int4 q_start, Int4 q_end,
                       Int4 s_start, Int4 s_end, BlastIntervalTree *tree);


    void
Blast_IntervalTreeReset(BlastIntervalTree *tree);
