// blast.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license
// agreement,
// provided that this statement is retained.
//
// Main code for blast

#include "blast.h"
//#define QUERY_BLK
char *query_arr[BATCH_SIZE];
char *queryDescription_arr[BATCH_SIZE];
struct PSSMatrix PSSMatrix_arr[BATCH_SIZE];

void blast_search(char *searchDbFile, struct PSSMatrix PSSMatrix, char *query);
void blast_search_multi(char *searchDbFile, struct PSSMatrix *PSSMatrix_arr,
        int numQuery, struct scoreMatrix scoreMatrix);

int4 main(int4 argc, char *argv[]) {
    char *query, *queryDescription;
    unsigned char queryAlphabetType, previousQueryAlphabetType = 10;
    struct scoreMatrix scoreMatrix;
    struct PSSMatrix PSSMatrix;
    int queryNum;

    // Process command line arguments
    parameters_processArguments(argc, argv);

    if (parameters_outputType != parameters_xml &&
            parameters_outputType != parameters_tabular) {
            printf("FSA-BLAST 1.05\n\n");
    }

    //struct PSSMatrix *PSSMatrix_arr = (struct PSSMatrix *)malloc(
    //parameters_batch_size * sizeof(struct PSSMatrix));

    // Read the first sequence from FASTA file (the query)
    readFasta_open(parameters_queryFile);
    if (!(readFasta_readSequence())) {
        fprintf(stderr, "Error reading query from FASTA file %s\n", argv[1]);
        exit(-1);
    }

    // Open sequence data file and read information
    readdb_open(parameters_subjectDatabaseFile);
    encoding_initialize(encoding_protein);
    parameters_loadDefaults(encoding_protein);
    scoreMatrix = scoreMatrix_load(parameters_scoringMatrixPath);

    int nextQuery = 0;
    do {
        int queryLettersBlk = 0;
        longestQueryLength = 0;
        queryNum = 0;
        do {
            // Initialize global variables
            //global_initialize();

            // Get the longest query in queryBatch
            longestQueryLength = max(longestQueryLength, readFasta_sequenceLength);

            // Make copy of the sequence
            query_arr[queryNum] = query =
                (char *)global_malloc(sizeof(char) * readFasta_sequenceLength + 1);
            strcpy(query, readFasta_sequenceBuffer);

            // Make copy of the description
            queryDescription = (char *)global_malloc(
                    sizeof(char) * (readFasta_descriptionLength + 1));
            blast_queryDescription_multi[queryNum] =  blast_queryDescription = (char *)global_malloc(
                    sizeof(char) * (readFasta_descriptionLength + 1));
            strcpy(queryDescription, readFasta_descriptionBuffer);
            strcpy(blast_queryDescription, readFasta_descriptionBuffer);

            // Determine the alphabet type of the query
            queryAlphabetType = encoding_determineAlphabetType(query, strlen(query));

            queryDescription_arr[queryNum] = queryDescription = print_formatDescription(queryDescription, 7, 0, 70);
            //if (parameters_outputType != parameters_xml &&
                    //parameters_outputType != parameters_tabular) {
                //printf("Query= %s\n", queryDescription);
                //printf("         (%lu letters)\n\n", strlen(query));
            //}

            // If not the same alphabet type as previous query, abort
            if (previousQueryAlphabetType < 10 &&
                    previousQueryAlphabetType != queryAlphabetType) {
                fprintf(stderr, "Error: Processing sequence %s\n", query);
                fprintf(stderr, "Error: Unable to process a mix of both protein and "
                        "nucleotide queries\n");
                fflush(stderr);
                exit(-1);
            }
            previousQueryAlphabetType = queryAlphabetType;

            // Initialize encoding
            //encoding_initialize(queryAlphabetType);

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

                fprintf(stderr, "Nucleotide sequence is not support\n");
                exit(0);
#if 0
                // Create a nucleotide scoring matrix use match and mismatch penalties
                scoreMatrix =
                    scoreMatrix_create(parameters_matchScore, parameters_mismatchScore);
                //		scoreMatrix_print(scoreMatrix);

                // Create the PSSMatrix
                PSSMatrix = PSSMatrix_create(scoreMatrix, query);
                //		PSSMatrix_print(PSSMatrix);

                nucleotideLookup_build(PSSMatrix, parameters_wordTableBytes);
                //		nucleotideLookup_print();
#endif		
            }
            // If a protein alphabet
            else {

                // Load the scoring matrix (eg. BLOSUM)
                //scoreMatrix = scoreMatrix_load(parameters_scoringMatrixPath);
                //		scoreMatrix_print(scoreMatrix);

                // Create the PSSMatrix
                PSSMatrix = PSSMatrix_create(scoreMatrix, query);
                //		PSSMatrix_print(PSSMatrix);

                // Use query sequence to build the word lookup FSA structure
                if (readdb_numberOfSequences != readdb_numberOfClusters) {
                    wordLookupDFA_build(PSSMatrix, encoding_sentinalCode,
                            parameters_wordSize);
                } else {
                    wordLookupDFA_build(PSSMatrix, encoding_numRegularLetters,
                            parameters_wordSize);
                }

                wordLookupDFA_assign(queryNum);
                PSSMatrix_arr[queryNum] = PSSMatrix;
            }

            //free(queryDescription);
            //free(blast_queryDescription);
            queryNum++;
            nextQuery = readFasta_readSequence();
            queryLettersBlk += PSSMatrix.length;

        } while (nextQuery && (queryNum < BATCH_SIZE) && queryLettersBlk < 60000 && queryNum < 1024);
            //} while (nextQuery && (queryNum < BATCH_SIZE) && queryLettersBlk);

        int ii;

        struct timeval start, end;
        gettimeofday(&start, NULL);

#if 0
        proteinLookup_query_initial(encoding_numRegularLetters, parameters_wordSize);
        for(ii = 0; ii < queryNum; ii++)
        {
            //proteinLookup_query_add(ii, PSSMatrix_arr[ii].queryCodes, PSSMatrix_arr[ii].length, parameters_wordSize, scoreMatrix);
            proteinLookup_query_add(ii, PSSMatrix_arr[ii].bestMatchCodes, PSSMatrix_arr[ii].length, parameters_wordSize, scoreMatrix);
        }
        proteinLookup_sort();
#else
        proteinLookup_query_initial_blk(encoding_numRegularLetters, 
                parameters_wordSize, numQueryBlk);

        for(ii = 0; ii < queryNum; ii++)
        {
            proteinLookup_query_add_blk(ii, PSSMatrix_arr[ii].queryCodes, 
                    PSSMatrix_arr[ii].length, parameters_wordSize, 
                    scoreMatrix, numQueryBlk);

            //queryLettersBlk += PSSMatrix_arr[ii].length;

            //if(queryLettersBlk >= 32768)
            //{

            //proteinLookup_sort_blk(numQueryBlk);
                //numQueryBlk++;
                //proteinLookup_query_initial_blk(encoding_numRegularLetters, 
                        //parameters_wordSize, numQueryBlk);

                //queryLettersBlk = 0;
            //}
        }

        proteinLookup_sort_blk(numQueryBlk);

#endif

        gettimeofday(&end, NULL);

        fprintf(stderr, "index building time: %f numQueryBlk: %d numQuery: %d\n", (float)((end.tv_sec * 1000000 + end.tv_usec)
                    - (start.tv_sec * 1000000 + start.tv_usec)) * 1e-6, numQueryBlk + 1, queryNum);


        blast_search_multi(parameters_subjectDatabaseFile, PSSMatrix_arr,
                queryNum, scoreMatrix);

        for(ii = 0; ii < queryNum; ii++)
        {
            PSSMatrix_free(PSSMatrix_arr[ii]);
            free(query_arr[ii]);
            free(queryDescription_arr[ii]);
            free(blast_queryDescription_multi[ii]);
        }

        queryIndex_free_blk();

    } while (nextQuery);


    encoding_free();
    // close FASTA reader
    readFasta_close();

    // Free all global data
    global_free();
    semiGappedScoring_free();
    oldSemiGappedScoring_free();
    oldGappedScoring_free();
    gappedScoring_free();
    nuGappedScoring_free();
    bytepackGappedScoring_free();
    fasterBytepackGappedScoring_free();
    gappedExtension_free();
    fasterGappedExtension_free();
    scoreMatrix_free(scoreMatrix);

    // close database
    readdb_close();
    //free(PSSMatrix_arr);
    parameters_free();
    return 0;
}

void blast_search_multi(char *searchDbFile, struct PSSMatrix *PSSMatrix_arr,
        int numQuery, struct scoreMatrix scoreMatrix) {
    printf("number of queries: %d longestQueryLength: %d\n", numQuery,
            longestQueryLength);

    global_initialize_multi();
    blast_numQuery = numQuery;
    //encoding_initialize(encoding_protein);
    //parameters_loadDefaults(encoding_protein);

    int ii;
    for (ii = 0; ii < numQuery; ii++) {
        statistics_initialize_multi(PSSMatrix_arr[ii], readdb_numberOfLetters,
                readdb_numberOfSequences, ii);

        blast_ungappedNominalTrigger_multi[ii] = get_ncbi_ungappedNominalTrigger(PSSMatrix_arr[ii], scoreMatrix);
        blast_gappedNominalCutoff_multi[ii] = get_ncbi_dropoff_score(PSSMatrix_arr[ii]);


        //double Lambda = 0.31760595763573063; 
        //double scale_factor = 32;
        //gap_matrix = Blast_MatrixInfoNew(BLASTAA_SIZE, BLASTAA_SIZE, FALSE);
        //s_MatrixInfoInit(gap_matrix, Lambda, scale_factor);
        //Blast_MatrixInfoFree(&gap_matrix);

        // Determine the minimum gapped nominal score required for reporting the
        // alignment
        //blast_gappedNominalCutoff_multi[ii] =
        //statistics_gappedEvalue2nominal_multi(parameters_cutoff, ii);

        // Determine the minimum/maximum semi-gapped scores to achieve cutoff
        blast_nominalR1cutoff_multi[ii] = ceil(
                (float)blast_gappedNominalCutoff_multi[ii] * parameters_semiGappedR1);
        blast_nominalR2cutoff_multi[ii] = ceil(
                (float)blast_gappedNominalCutoff_multi[ii] * parameters_semiGappedR2);

        //blast_ungappedNominalTrigger_multi[ii] =
        //statistics_ungappedNormalized2nominal(
        //parameters_ungappedNormalizedTrigger);

        // Gapping trigger cannot be greater than final cutoff
        if (blast_ungappedNominalTrigger_multi[ii] >
                blast_gappedNominalCutoff_multi[ii])
            blast_ungappedNominalTrigger_multi[ii] =
                blast_gappedNominalCutoff_multi[ii];
    }

    ungappedExtension_initialize_multi();

    // Initialize collections of alignments
    alignments_initialize_multi();

    //for (ii = 0; ii < numQuery; ii++) {
        //hitMatrix_initialize_multi(longestQueryLength, readdb_longestSequenceLength,
                //readdb_sequences, ii);
    //}

    total_numberOfLetters = readdb_numberOfLetters;
    alignments_queryIdx(PSSMatrix_arr, scoreMatrix, query_arr, queryDescription_arr, numQuery);
    
    //for (ii = 0; ii < numQuery; ii++) {

#if 0
        printf("Query %d alignment info:\n", ii);
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
#endif

        //print_gappedAlignmentsBrief_multi(ii);
        // print_gappedAlignmentsFull_multi(query_arr[ii], PSSMatrix_arr[ii], ii);

        //hitMatrix_free_multi(ii);
        //wordLookupDFA_free_multi(ii);
        //semiGappedScoring_free_multi(ii);
        //fasterGappedExtension_free_multi(ii);
        //gappedExtension_free_multi(ii);
        //alignments_free_multi(ii);
        //
        //}

    //encoding_free();
}

void blast_search(char *searchDbFile, struct PSSMatrix PSSMatrix, char *query) {
    char *indexFilename;
    int4 tickFrequency;

    // Construct sequence filename
    indexFilename = (char *)global_malloc(strlen(searchDbFile) + 9);
    sprintf(indexFilename, "%s.index", searchDbFile);

    // Check if index file exists. If not, disable use of index
    /*	if ((indexFile = fopen(indexFilename, "r")) != NULL)
        fclose(indexFile);
        else*/
    parameters_useIndex = 0;

    // Check that alphabet type of query and database match
    if (encoding_alphabetType == encoding_protein &&
            readdb_dbAlphabetType == encoding_nucleotide) {
        fprintf(stderr, "Error: database %s contains nucleotide sequences\n",
                searchDbFile);
        fprintf(stderr, "Error: searching a nucleotide database with a protein "
                "query is not supported\n\n");
        exit(-1);
    }
    if (encoding_alphabetType == encoding_nucleotide &&
            readdb_dbAlphabetType == encoding_protein) {
        fprintf(stderr, "Error: database %s contains protein sequences\n",
                searchDbFile);
        fprintf(stderr, "Error: searching a protein database with a nucleotide "
                "query is not supported\n\n");
        exit(-1);
    }

    // Determine tick frequence
    tickFrequency = ceil((float)readdb_numberOfSequences / 50.0);

    // Initialize BLAST statistics (calculate log(2), log(K), nominal drop-offs,
    // etc.)
    statistics_initialize(PSSMatrix, readdb_numberOfLetters,
            readdb_numberOfSequences);

    // Determine the minimum gapped nominal score required for reporting the
    // alignment
    blast_gappedNominalCutoff =
        statistics_gappedEvalue2nominal(parameters_cutoff);

    // Determine the minimum/maximum semi-gapped scores to achieve cutoff
    blast_nominalR1cutoff =
        ceil((float)blast_gappedNominalCutoff * parameters_semiGappedR1);
    blast_nominalR2cutoff =
        ceil((float)blast_gappedNominalCutoff * parameters_semiGappedR2);

    // Determine the minimum ungapped nominal score required to trigger gapping
    if (encoding_alphabetType == encoding_protein) {
        blast_ungappedNominalTrigger = statistics_ungappedNormalized2nominal(
                parameters_ungappedNormalizedTrigger);
    } else {
        blast_ungappedNominalTrigger =
            statistics_ungappedNucleotideTrigger(PSSMatrix);
    }

    // Gapping trigger cannot be greater than final cutoff
    if (blast_ungappedNominalTrigger > blast_gappedNominalCutoff)
        blast_ungappedNominalTrigger = blast_gappedNominalCutoff;

    // Initialize collections of alignments
    alignments_initialize();

    // Initialize collections of ungapped extensions
    ungappedExtension_initialize();

    if (parameters_outputType != parameters_xml &&
            parameters_outputType != parameters_tabular) {
        printf("Database: %s\n", searchDbFile);
        printf("           %s sequences;",
                global_int4toString(readdb_numberOfSequences));
        printf(" %s total letters\n\n",
                global_int8toString(readdb_numberOfLetters));

#ifndef VERBOSE
        printf("Searching...");
        fflush(stdout);
#endif
    }

    // Initialize the hitMatrix
    hitMatrix_initialize(PSSMatrix.length, readdb_longestSequenceLength,
            readdb_sequences);

    blast_prepTime = clock();
    blast_searchTime = -clock();

    while (1) {
        // Only one hit required to trigger ungapped extension
        if (parameters_oneHitTrigger) {
            search_protein1hit(PSSMatrix, readdb_sequenceData,
                    readdb_numVolumeSequences, tickFrequency);
        }
        // Two hits to trigger an ungapped extensions
        else {

            search_protein2hit(PSSMatrix, readdb_sequenceData,
                    readdb_numVolumeSequences, tickFrequency);
        }

        if (readdb_volume + 1 < readdb_numberOfVolumes) {
#ifndef NO_STAGE3
            // Before loading next volume, perform initial semi-gapped or bytepacked
            // alignment
            // on high-scoring ungapped extensions in this volume
            blast_searchTime += clock();
            blast_semiGappedScoreTime -= clock();
            alignments_findGoodAlignments(PSSMatrix);
            blast_semiGappedScoreTime += clock();
            blast_searchTime -= clock();
#endif

            // Copy subject sequences from good alignments into memory
            blast_searchTime += clock();
            blast_copyTime -= clock();
            alignments_loadSubjectsIntoMemory(PSSMatrix);
            blast_copyTime += clock();
            blast_searchTime -= clock();

            // Load the next volume
            readdb_nextVolume();

            // Re-initialize the hitMatrix
            hitMatrix_reinitialize(PSSMatrix.length, readdb_longestSequenceLength,
                    readdb_sequences);
        } else
            break;
    }

    if (parameters_outputType != parameters_xml &&
            parameters_outputType != parameters_tabular) {
#ifndef VERBOSE
        printf("done.\n\n\n\n");
        fflush(stdout);
#endif
    }

    blast_searchTime += clock();

    //	blast_compareScorings(PSSMatrix);
    //	exit(0);

#ifndef NO_STAGE3
        // Perform semi-gapped / bytepacked alignment to find good alignments
        blast_semiGappedScoreTime -= clock();
        alignments_findGoodAlignments(PSSMatrix);
        blast_semiGappedScoreTime += clock();

        // Perform gapped alignment to find final alignments
        blast_gappedScoreTime -= clock();
        alignments_findFinalAlignments(PSSMatrix);
        blast_gappedScoreTime += clock();
#endif

#ifndef NO_STAGE4
    blast_gappedExtendTime -= clock();

    // Read the final alignment subject descriptions
    alignments_getFinalAlignmentDescriptions();

    // Find traceback information
    alignments_getTracebacks(PSSMatrix);
    blast_gappedExtendTime += clock();
#endif

    blast_finalizeTime -= clock();

    // Print alignments
    if (alignments_finalAlignments->numEntries == 0 &&
            parameters_outputType != parameters_xml &&
            parameters_outputType != parameters_tabular) {
        printf("\n ***** No hits found ******\n");
    } else {
#ifndef NO_STAGE4
        if (parameters_outputType == parameters_xml) {
            print_XMLheader(query, PSSMatrix);
            print_gappedAlignmentsFull(query, PSSMatrix);
            print_XMLfooter();
        } else if (parameters_outputType == parameters_tabular) {
            print_gappedAlignmentsFull(query, PSSMatrix);
        } else {
            print_gappedAlignmentsBrief();
            print_gappedAlignmentsFull(query, PSSMatrix);
        }
#endif
    }

    if (parameters_outputType != parameters_xml &&
            parameters_outputType != parameters_tabular) {
        if (readdb_numberOfVolumes > 0)
            printf("  Database: %s  (%d volumes)\n", searchDbFile,
                    readdb_numberOfVolumes);
        else
            printf("  Database: %s\n", searchDbFile);
        //    printf("    Posted date:  Apr 5, 2004  5:12 PM\n");
        printf("  Number of letters in database: %s\n",
                global_int8toString(statistics_databaseSize));
        printf("  Number of sequences in database:  %u\n",
                readdb_numberOfSequences);

        printf("\nLambda     K      H     (ungapped)");
        printf("\n %.3f     %.3f  %.3f", statistics_ungappedLambda,
                statistics_ungappedK, statistics_ungappedH);
        printf("\n\nLambda     K      H     (gapped)");
        printf("\n %.3f     %.3f  %.3f", statistics_gappedParams.lambda,
                statistics_gappedParams.K, statistics_gappedParams.H);
        printf("\n\n\nMatrix: %s", parameters_scoringMatrix);
        printf("\nGap Penalties: Existence: %d, Extension: %d", parameters_startGap,
                parameters_extendGap);
        if ((parameters_semiGappedScoring || parameters_bytepackedScoring))
            printf("\nSemi-Gapped Gap Penalties: Existence: %d, Extension: %d",
                    parameters_semiGappedStartGap, parameters_semiGappedExtendGap);
            printf("\nNumber of Hits to DB: %s", global_int4toString(blast_numHits));
        printf("\nNumber of Sequences: %u\n", readdb_numberOfSequences);
            printf("Number of extensions: %u\n", blast_numUngappedExtensions);
            printf("Number of successful extensions: %u\n",
                    blast_numTriggerExtensions);
            printf("Number of sequences with successful extensions: %u\n",
                    blast_numTriggerSequences);
        if ((parameters_semiGappedScoring || parameters_bytepackedScoring ||
                    parameters_tableScoring))
            printf("Number of sequences with semi-gapped score above cutoff: %u\n",
                    blast_numGoodAlignments);
        printf("Number of sequences better than %g: %u\n", parameters_cutoff,
                alignments_finalAlignments->numEntries);
            if (parameters_semiGappedScoring || parameters_bytepackedScoring ||
                    parameters_tableScoring)
                printf("Number of HSP's that attempted semi-gapping: %u\n",
                        blast_numSemiGapped);
            printf("Number of HSP's that attempted gapping: %u\n", blast_numGapped);
            printf("Number of HSP's contained and not gapped: %u\n",
                    blast_numExtensionsPruned);
            printf("Number of HSP's succeeded/attempted join: %u/%u\n",
                    blast_numSuccessfullyJoined, blast_numAttemptedJoin);
        if (blast_numExpandedSequences)
            printf("Number of cluster members recreated = %d\n",
                    blast_numExpandedSequences);
        printf("Total subject bytes copied/unpacked = %d/%d\n", blast_totalCopied,
                blast_totalUnpacked);
        printf("length of query: %u\n", statistics_querySize);
        printf("length of database: %s\n",
                global_int8toString(statistics_databaseSize));
        printf("effective HSP length: %u\n", statistics_lengthAdjust);
        printf("effective length of query: %u\n", statistics_effectiveQuerySize);
        printf("effective length of database: %s\n",
                global_int8toString(statistics_effectiveDatabaseSize));
        printf("effective search space: %llu\n", statistics_searchSpaceSize);
        printf("effective search space used: %llu\n", statistics_searchSpaceSize);

        if (encoding_alphabetType == encoding_protein) {
            printf("T: %d\n", parameters_T);
            printf("A: %d\n", parameters_A);
        }
        printf("X1: %d\n", statistics_ungappedNominalDropoff);
        printf("X2: %d\n", statistics_gappedNominalDropoff);
        printf("X3: %d\n", statistics_gappedFinalNominalDropoff);
        printf("S1: %d\n", blast_ungappedNominalTrigger);
        printf("S2: %d\n", blast_gappedNominalCutoff);
        if (blast_dynamicGappedNominalCutoff > 0)
            printf("S3: %d\n", blast_dynamicGappedNominalCutoff);
        printf("F2: %d\n", blast_nominalR1cutoff);
        if (blast_dynamicNominalR1cutoff > 0)
            printf("F3: %d\n", blast_dynamicNominalR1cutoff);

        //    	printf("Total malloced=%s\n",
        // global_int4toString(global_totalMalloc));
    }

    // Free memory used by hitMatrix, PSSMatrix, alignments and sequence filename
    hitMatrix_free();
    alignments_free();

    blast_finalizeTime += clock();
}
