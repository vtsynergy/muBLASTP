// blast.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license
// agreement,
// provided that this statement is retained.
//
// Main code for blast

#include "blast.h"


char *query_arr[BATCH_SIZE];
char *queryDescription_arr[BATCH_SIZE];
struct PSSMatrix PSSMatrix_arr[BATCH_SIZE];

void print_info(int queryNum, struct PSSMatrix *PSSMatrix_arr);

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

  //if (parameters_outputType != parameters_xml &&
      //parameters_outputType != parameters_tabular) {
    //if (parameters_ssearch)
      //printf("FSA-SSEARCH 1.05\n\n");
    //else
      //printf("FSA-BLAST 1.05\n\n");
  //}

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

  //read_dbLookupAux(parameters_subjectDatabaseFile);
#ifndef COMPRESS_INDEX
  read_dbLookup(parameters_subjectDatabaseFile);
#else
  read_dbLookup_cp(parameters_subjectDatabaseFile);
#endif

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
      longestQueryLength = max(longestQueryLength, readFasta_sequenceLength);

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

        // Create a nucleotide scoring matrix use match and mismatch penalties
        scoreMatrix =
            scoreMatrix_create(parameters_matchScore, parameters_mismatchScore);
        //		scoreMatrix_print(scoreMatrix);

        // Create the PSSMatrix
        PSSMatrix = PSSMatrix_create(scoreMatrix, query);
        //		PSSMatrix_print(PSSMatrix);

        nucleotideLookup_build(PSSMatrix, parameters_wordTableBytes);
        //		nucleotideLookup_print();
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
      //free(queryDescription);
      nucleotideLookup_free();
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
  semiGappedScoring_free();
  oldSemiGappedScoring_free();
  oldGappedScoring_free();
  gappedScoring_free();
  nuGappedScoring_free();
  bytepackGappedScoring_free();
  fasterBytepackGappedScoring_free();
  gappedExtension_free();
  fasterGappedExtension_free();

#ifdef NEIGHBOR_INDEX
  neighbourLookup_free();
#endif

#ifndef COMPRESS_INDEX
  free_dbindex();
#else
  free_indexdb_cp();
#endif
  scoreMatrix_free(scoreMatrix);

  // close database
  readdb_close();
  //free(PSSMatrix_arr);
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
    // Determine the minimum gapped nominal score required for reporting the
    // alignment

    // Determine the minimum/maximum semi-gapped scores to achieve cutoff
#ifndef NCBI_BLAST
    blast_ungappedNominalTrigger_multi[ii] =
        statistics_ungappedNormalized2nominal(
            parameters_ungappedNormalizedTrigger);

    blast_gappedNominalCutoff_multi[ii] =
        statistics_gappedEvalue2nominal_multi(parameters_cutoff, ii);


        if (blast_ungappedNominalTrigger_multi[ii] >
                blast_gappedNominalCutoff_multi[ii])
            blast_ungappedNominalTrigger_multi[ii] =
                blast_gappedNominalCutoff_multi[ii];


#else
    blast_ungappedNominalTrigger_multi[ii] = get_ncbi_ungappedNominalTrigger(PSSMatrix_arr[ii], scoreMatrix);
    blast_gappedNominalCutoff_multi[ii] = get_ncbi_dropoff_score(PSSMatrix_arr[ii]);


    double Lambda = 0.31760595763573063; 
    double scale_factor = 32;
    gap_matrix = Blast_MatrixInfoNew(BLASTAA_SIZE, BLASTAA_SIZE, FALSE);
    s_MatrixInfoInit(gap_matrix, Lambda, scale_factor);
    Blast_MatrixInfoFree(&gap_matrix);
    //printf("blast_gappedNominalCutoff_multi: %d new_cutoffscore: %d\n", blast_gappedNominalCutoff_multi[ii], cutoffscore);
    //printf("cutoff: %d\n", cutoff);
#endif
    //fprintf(stderr, "blast_ungappedNominalTrigger: %d blast_gappedNominalCutoff: %d\n", blast_ungappedNominalTrigger_multi[ii], blast_gappedNominalCutoff_multi[ii]);

    blast_nominalR1cutoff_multi[ii] = ceil(
        (float)blast_gappedNominalCutoff_multi[ii] * parameters_semiGappedR1);
    blast_nominalR2cutoff_multi[ii] = ceil(
        (float)blast_gappedNominalCutoff_multi[ii] * parameters_semiGappedR2);


    //printf("blast_ungappedNominalTrigger_multi: %d blast_gappedNominalCutoff_multi: %d\n", blast_ungappedNominalTrigger_multi[ii], blast_gappedNominalCutoff_multi[ii]);
  }

  //printf("min_lambda: %f\n", min_lambda);

  ungappedExtension_initialize_multi2();

  // Initialize collections of alignments
  alignments_initialize_multi2();

  if(parameters_num_threads > 1)
  {
      alignments_db_omp(PSSMatrix_arr, scoreMatrix, query_arr, queryDescription_arr, numQuery);
  }
  else
  {
      alignments_db_serial(PSSMatrix_arr, scoreMatrix, query_arr, queryDescription_arr, numQuery);
  }

  
  //alignments_global_free_multi2();
  encoding_free();
}

void print_info(int queryNum, struct PSSMatrix *PSSMatrix_arr) {
  int4 count = 0;
  struct alignment *currentAlignment;
  struct finalAlignment *finalAlignment;
  struct gappedExtension *currentExtension;

  while (count < alignments_finalAlignments_multi[queryNum]->numEntries &&
         count < parameters_numDisplayAlignments) {
    finalAlignment = memSingleBlock_getEntry(
        alignments_finalAlignments_multi[queryNum], count);
    currentAlignment = finalAlignment->alignment;
    uint4 sequenceCount = currentAlignment->sequenceCount;

    printf("%d %d %f\n", readdb_sequenceData[sequenceCount].sequenceLength,
           readdb_sequenceData[sequenceCount].sequenceLength -
               PSSMatrix_arr[queryNum].length,
           (float)((float)readdb_sequenceData[sequenceCount].sequenceLength -
                   PSSMatrix_arr[queryNum].length) /
               (PSSMatrix_arr[queryNum].length));
    count++;
  }
}



