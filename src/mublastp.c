// blast.c

#include "blast.h"


char *query_arr[BATCH_SIZE];
char *queryDescription_arr[BATCH_SIZE];

struct PSSMatrix PSSMatrix_arr[BATCH_SIZE];

void blast_search_db(char *searchDbFile, struct PSSMatrix *PSSMatrix_arr,
        int numQuery, struct scoreMatrix scoreMatrix);

int4 main(int4 argc, char *argv[]) {
    //VTPauseSampling();
    char *query, *queryDescription;
    unsigned char queryAlphabetType, previousQueryAlphabetType = 10;
    struct scoreMatrix scoreMatrix;
    struct PSSMatrix PSSMatrix;
    int queryNum;

#ifdef SSEARCH
    parameters_ssearch = 1;
#endif

    // Process command line arguments
    parameters_processArguments(argc, argv);

    // Read the first sequence from FASTA file (the query)
    readFasta_open(parameters_queryFile);
    if (!(readFasta_readSequence())) {
        fprintf(stderr, "Error reading query from FASTA file %s\n", argv[1]);
        exit(-1);
    }

    // Open sequence data file and read information
    readdb_open_mem(parameters_subjectDatabaseFile);

    read_dbLookup(parameters_subjectDatabaseFile);

#ifdef NEIGHBOR_INDEX
    neighbourLookup_init();
#endif

    encoding_initialize(encoding_protein);
    parameters_loadDefaults(encoding_protein);
    scoreMatrix = scoreMatrix_load(parameters_scoringMatrixPath);
    encoding_free();


    int nextQuery = 0;
    do {

        longestQueryLength = 0;
        queryNum = 0;

        do {

            // Get the longest query in queryBatch
            longestQueryLength = MAX(longestQueryLength, readFasta_sequenceLength);

            // Make copy of the sequence
            query_arr[queryNum] = query =
                (char *)global_malloc(sizeof(char) * readFasta_sequenceLength + 1);
            strcpy(query, readFasta_sequenceBuffer);

            // Make copy of the description
            queryDescription_arr[queryNum] = queryDescription = (char *)global_malloc(
                    sizeof(char) * (readFasta_descriptionLength + 1));
            blast_queryDescription_multi[queryNum] = blast_queryDescription = (char *)global_malloc(
                    sizeof(char) * (readFasta_descriptionLength + 1));
            strcpy(queryDescription, readFasta_descriptionBuffer);
            strcpy(blast_queryDescription, readFasta_descriptionBuffer);

            // Determine the alphabet type of the query
            queryAlphabetType = encoding_determineAlphabetType(query, strlen(query));

            if (queryAlphabetType == encoding_nucleotide) {
                fprintf(stderr, "Error: Processing sequence %s\n", query);
                fprintf(
                        stderr,
                        "Error: DB_INDEX search unable to process nucleotide queries\n");
                fflush(stderr);
                exit(-1);
            }

            previousQueryAlphabetType = queryAlphabetType;

            // Initialize encoding
            encoding_initialize(queryAlphabetType);

            // Filter the query using DUST or SEG
            if (parameters_filterEnabled) {
                if (queryAlphabetType == encoding_protein)
                    seg_segSequence(query);
                else
                    dust_dustSequence(query);
            }

            // Load parameter defaults based on query alphabet type
            //parameters_loadDefaults(queryAlphabetType);

            // If a nucleotide alphabet
            if (queryAlphabetType == encoding_nucleotide) {
                fprintf(stderr, "Do not support nucleotide sequence\n");
                exit(0);
            }
            // If a protein alphabet
            else {

                PSSMatrix_arr[queryNum] = PSSMatrix_create(scoreMatrix, query);

#ifdef NEIGHBOR_INDEX
                neighbourLookup_build(PSSMatrix_arr[queryNum], scoreMatrix,
                        parameters_wordSize);
#endif
            }

            encoding_free();
            queryNum++;
            nextQuery = readFasta_readSequence();
        } while (nextQuery && (queryNum < BATCH_SIZE));

        blast_search_db(parameters_subjectDatabaseFile, PSSMatrix_arr, queryNum,
                scoreMatrix);

        int ii;
        for(ii = 0; ii < queryNum; ii++)
        {
            PSSMatrix_free(PSSMatrix_arr[ii]);
            free(query_arr[ii]);
            free(queryDescription_arr[ii]);
            free(blast_queryDescription_multi[ii]);
        }

    } while (nextQuery);



    // close FASTA reader
    readFasta_close();

    // Free all global data
    global_free();

#ifdef NEIGHBOR_INDEX
    neighbourLookup_free();
#endif

    scoreMatrix_free(scoreMatrix); 

    free_dbindex();
    readdb_close_mem();

    parameters_free();
    return 0;
}

void blast_search_db(char *searchDbFile, struct PSSMatrix *PSSMatrix_arr,
        int numQuery, struct scoreMatrix scoreMatrix) {
    fprintf(stderr, "number of queries: %d longestQueryLength: %d\n", numQuery,
            longestQueryLength);

    global_initialize_multi();
    blast_numQuery = numQuery;
    encoding_initialize(encoding_protein);
    //parameters_loadDefaults(encoding_protein);

    int ii;
    for (ii = 0; ii < numQuery; ii++) {
        statistics_initialize_multi(PSSMatrix_arr[ii], readdb_numberOfLetters,
                readdb_numberOfSequences, ii);

        // Determine the minimum/maximum semi-gapped scores to achieve cutoff
        blast_ungappedNominalTrigger_multi[ii] = get_ncbi_ungappedNominalTrigger(PSSMatrix_arr[ii], scoreMatrix);
        blast_gappedNominalCutoff_multi[ii] = get_ncbi_dropoff_score(PSSMatrix_arr[ii]);


        double Lambda = 0.31760595763573063; 
        double scale_factor = 32;
        gap_matrix = Blast_MatrixInfoNew(BLASTAA_SIZE, BLASTAA_SIZE, FALSE);
        s_MatrixInfoInit(gap_matrix, Lambda, scale_factor);
        Blast_MatrixInfoFree(&gap_matrix);

        blast_nominalR1cutoff_multi[ii] = ceil(
                (float)blast_gappedNominalCutoff_multi[ii] * parameters_semiGappedR1);
        blast_nominalR2cutoff_multi[ii] = ceil(
                (float)blast_gappedNominalCutoff_multi[ii] * parameters_semiGappedR2);

    }

    ungappedExtension_initialize_multi2();

    // Initialize collections of alignments
    alignments_initialize_multi2();

    total_numberOfLetters = readdb_numberOfLetters;

    alignments_dbIdx(PSSMatrix_arr, scoreMatrix, 
            query_arr, queryDescription_arr, numQuery);

    encoding_free();
}




