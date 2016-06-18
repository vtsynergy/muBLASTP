// blast.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license
// agreement,
// provided that this statement is retained.
//
// Main code for blast

#include "blast.h"
#include <unistd.h>

int4 main(int4 argc, char *argv[]) {

    char *filename = NULL;
    dbIdx_block_size = 0;

    int c;
    while((c = getopt(argc, argv, "i:s:")) != -1)
    {
        switch (c)
        {
            case 'i':
                filename = optarg;
                break;
            case 's':
                dbIdx_block_size = atoi(optarg);
                break;
            default:
                fprintf(stderr, "Useage: indexdb -i <Database>"
                        " -s [Block size, default 128(K)]\n");
                exit(-1);
        }
    }

    if(filename == NULL)
    {
        fprintf(stderr, "Useage: indexdb -i <Database>"
                " -s [Block size, default 128(K)]\n");
        exit(-1);
    }

    if(dbIdx_block_size == 0)
    {
        dbIdx_block_size = 128 * 1024;
    }
    else
    {
        dbIdx_block_size *= 1024;
    }

    
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
