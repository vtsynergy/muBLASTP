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
    return readdb_sequenceData[seqId_a].sequenceLength 
        > readdb_sequenceData[seqId_b].sequenceLength; 
}


char* getSequence(uint4 seqId)
{
    // Declare memory for the sequence
    int ii;
    for(ii = 0; ii < readdb_sequenceData[seqId].sequenceLength; ii++)
    {
        sequenceBuffer[ii] = 
            encoding_lettersArray[readdb_sequenceData[seqId].sequence[ii]];
    }
    sequenceBuffer[ii] = '\0';

    return sequenceBuffer;
}

void print_sequence(int seqId, FILE *output_file)
{

#ifdef DESCIPT_IN_MEM
    char *seqDes = 
        descriptions_getDescription_mem(
                readdb_sequenceData[seqId].descriptionStart, 
                readdb_sequenceData[seqId].descriptionLength);
#else
    char *seqDes = descriptions_getDescription(
            readdb_sequenceData[seqId].descriptionStart, 
            readdb_sequenceData[seqId].descriptionLength);
#endif

    char *seq = getSequence(seqId); 

    fprintf(output_file, ">%s\n%s\n", seqDes, seq);
    //fprintf(output_file, ">gi|%d\n%s\n", seqId, seq);

    free(seqDes); 
}

int4 main(int4 argc, char* argv[])
{
    // User must provide FASTA format file at command line
	if (argc < 3)
	{
		fprintf(stderr, "Useage: indexdb <DB filename> <Output filename>\n");
		exit(-1);
	}

	char *ifilename = argv[1];
	char *ofilename = argv[2];

    // Open sequence data file and read information
	encoding_initialize(encoding_protein);
    readdb_open_mem(ifilename);


    FILE *output_file;

    //output_file = stdout;
    output_file = fopen(ofilename, "w");

    int *seqId = (int *)malloc(sizeof(int) * readdb_numberOfSequences);

    int startSeq = 0, endSeq = readdb_numVolumeSequences;

    sequenceBuffer = (char *)global_malloc(
            sizeof(char) * (readdb_longestSequenceLength + 1));

    while(1)
    {

        fprintf(stderr, "volumn: %d\n", readdb_volume);

        int ii;
        for(ii = startSeq; ii < endSeq; ii++)
        {
            seqId[ii] = ii;
        }

        qsort(seqId + startSeq, readdb_numVolumeSequences, sizeof(int), compareSeq);

        int numSeqPercent = readdb_numVolumeSequences / 100;

        for(ii = startSeq; ii < endSeq; ii++)
        {
            if(!((ii - startSeq) % numSeqPercent))
            {
                fprintf(stderr, "%d.", (ii - startSeq) / numSeqPercent);
            }

            print_sequence(seqId[ii], output_file);
        }

        fprintf(stderr, "\n");

        if (readdb_volume + 1 < readdb_numberOfVolumes) {
            readdb_nextVolume();
            startSeq = endSeq;
            endSeq += readdb_numVolumeSequences;
        }
        else
        {
            break;
        }
    }

    free(sequenceBuffer);
    fclose(output_file);
    free(seqId);


    //close database
    readdb_close_mem();
    return 0;
}
