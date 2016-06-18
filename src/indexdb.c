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
    if (argc == 3) {
        dbIdx_block_size = atoi(argv[2]) * 1024;
    }
    else if(argc == 2)
    {
        dbIdx_block_size = 128 * 1024;
    }
    else
    {
        fprintf(stderr,
                "Useage: indexdb <DB filename> <DB index block size (K letters)>\n");
        exit(-1);

    }

  char *filename = argv[1];

  // Open sequence data file and read information
  readdb_open_mem(filename);

  encoding_initialize(encoding_protein);
  struct scoreMatrix scoreMatrix;

  parameters_wordSize = 3;

  parameters_findScoringMatrix();
  scoreMatrix = scoreMatrix_load(parameters_scoringMatrixPath);

  while(1)
  {

#ifndef COMPRESSED_INDEX
      proteinLookup_db_build(encoding_numRegularLetters, parameters_wordSize,
              scoreMatrix, filename);
#else
      proteinLookup_db_cp_build(encoding_numRegularLetters, parameters_wordSize,
              scoreMatrix, filename);
#endif


      if(readdb_volume + 1 < readdb_numberOfVolumes)
      {
          readdb_nextVolume_mem();
      }
      else
      {
          break;
      }
  }


  // write_dbLookup(filename);
  scoreMatrix_free(scoreMatrix);
  encoding_free();

  // close database
  readdb_close_mem();
  return 0;
}
