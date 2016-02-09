
#define BLASTAA_SIZE 28 

typedef struct Blast_MatrixInfo {
    char * matrixName;         /**< name of the matrix */
    int    **startMatrix;     /**< Rescaled values of the original matrix */
    double **startFreqRatios;  /**< frequency ratios used to calculate matrix
                                 scores */
    int      rows;             /**< the number of rows in the scoring
                                 matrix. */
    int      cols;             /**< the number of columns in the scoring
                                 matrix, i.e. the alphabet size. */
    int      positionBased;    /**< is the matrix position-based */
    double   ungappedLambda;   /**< ungapped Lambda value for this matrix
                                 in standard context */
} Blast_MatrixInfo;

int
s_MatrixInfoInit(
        Blast_MatrixInfo * self,
        /*BLAST_SequenceBlk* queryBlk,*/
        /*BlastScoreBlk* sbp,*/
        double Lambda,
        double scale_factor
        //const char * matrixName
        );


Blast_MatrixInfo *
Blast_MatrixInfoNew(int rows, int cols, int positionBased);

Blast_MatrixInfo *gap_matrix;

void Blast_MatrixInfoFree(Blast_MatrixInfo ** ss);


double **
Nlm_DenseMatrixNew(int nrows,
        int ncols);

static void
s_RoundScoreMatrix(int **matrix, int rows, int cols,    
        double **floatScoreMatrix);


Blast_MatrixInfo *
Blast_MatrixInfoNew2(Blast_MatrixInfo *ss, int rows, int cols, int positionBased);


void Blast_MatrixInfoFree2(Blast_MatrixInfo * ss);
