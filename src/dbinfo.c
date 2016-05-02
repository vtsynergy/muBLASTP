// blast.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Main code for blast
#include "blast.h"
#include <stdio.h>
#include <errno.h>

char *sequenceBuffer;

int compareSeq(const void *a, const void *b)
{
    int seqId_a = *(int *)a;
    int seqId_b = *(int *)b;
    ASSERT(readdb_sequenceData[seqId_a].sequenceLength > 0 && readdb_sequenceData[seqId_b].sequenceLength > 0);
    return readdb_sequenceData[seqId_a].sequenceLength > readdb_sequenceData[seqId_b].sequenceLength; 
}

int4 main(int4 argc, char* argv[])
{
    // User must provide FASTA format file at command line
	if (argc < 2)
	{
		fprintf(stderr, "Useage: indexdb <DB filename>\n");
		exit(-1);
	}

	char *ifilename = argv[1];

    // Open sequence data file and read information
	encoding_initialize(encoding_protein);
    readdb_open(ifilename);

    fprintf(stderr, "numSeq: %d numLetters: %d\n", readdb_numberOfSequences, readdb_numberOfLetters);

    int *seqId = (int *)malloc(sizeof(int) * readdb_numberOfSequences);


    sequenceBuffer = (char *)global_malloc(
            sizeof(char) * (readdb_longestSequenceLength + 1));
    int ii;
    for(ii = 0; ii < readdb_numberOfSequences; ii++)
    {
        seqId[ii] = ii;
    }

    int seqLenCnt[20] = {0};

    while(1)
    {

        fprintf(stderr, "volumn: %d start: %d end: %d\n", readdb_volume, readdb_volumeOffset, readdb_volumeOffset + readdb_numVolumeSequences);

        if (readdb_volume + 1 < readdb_numberOfVolumes) {
            readdb_nextVolume();
        }
        else
        {
            break;
        }
    }

    qsort(seqId, readdb_numberOfSequences, sizeof(int), compareSeq);

    int longSeqCnt = 0;

    for(ii = 0; ii < readdb_numberOfSequences; ii++)
    {
        if(readdb_sequenceData[seqId[ii]].sequenceLength < 2000)
        {
            seqLenCnt[readdb_sequenceData[seqId[ii]].sequenceLength/100]++;
        }
        else
        {
            longSeqCnt++;
        }
    }

    fprintf(stderr, "min: %d max: %d median: %d aver: %d\n", readdb_sequenceData[0].sequenceLength, readdb_longestSequenceLength, readdb_sequenceData[seqId[readdb_numberOfSequences/2]].sequenceLength, readdb_numberOfLetters/readdb_numberOfSequences);

    for(ii = 0; ii < 20; ii++)
    {
        fprintf(stderr, "%d : %f\n", ii * 100, (float)seqLenCnt[ii]/readdb_numberOfSequences);
    }

    fprintf(stderr, "> %d : %f\n", 2000, (float)longSeqCnt/readdb_numberOfSequences);

    free(sequenceBuffer);
    free(seqId);


    //close database
    readdb_close();
    return 0;
}

