// blast.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license
// agreement,
// provided that this statement is retained.
//
// Main code for blast
#include "blast.h"
int4 main(int4 argc, char *argv[]) {
  // User must provide FASTA format file at command line
  if (argc < 3) {
    fprintf(
        stderr,
        "Useage: indexdb <DB filename> <DB index block size (K letters)>\n");
    exit(-1);
  }

  char *filename = argv[1];
  dbIdx_block_size = atoi(argv[2]) * 1024;

  // Open sequence data file and read information
  readdb_open(filename);

  encoding_initialize(encoding_protein);
  struct scoreMatrix scoreMatrix;

  parameters_wordSize = 3;

  parameters_findScoringMatrix();
  scoreMatrix = scoreMatrix_load(parameters_scoringMatrixPath);
  // scoreMatrix_print(scoreMatrix);

#ifndef COMPRESS_INDEX
  if (readdb_numberOfSequences != readdb_numberOfClusters) {
    proteinLookup_db_build(encoding_sentinalCode, parameters_wordSize,
                           scoreMatrix, filename);
  } else {
    proteinLookup_db_build(encoding_numRegularLetters, parameters_wordSize,
                           scoreMatrix, filename);
  }
#else
  if (readdb_numberOfSequences != readdb_numberOfClusters) {
      proteinLookup_db_cp_build(encoding_sentinalCode, parameters_wordSize,
              scoreMatrix, filename);
  } else {
      proteinLookup_db_cp_build(encoding_numRegularLetters, parameters_wordSize,
              scoreMatrix, filename);
  }

  //write_dbLookup_cp(filename);

  //free_indexdb_cp_temp();
  //free_indexdb_cp();
#endif

  // write_dbLookup(filename);
  scoreMatrix_free(scoreMatrix);
  encoding_free();

  // close database
  readdb_close();
  return 0;
}
