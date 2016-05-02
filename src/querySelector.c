// blast.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Main code for blast
#include "blast.h"
#include <stdio.h>
#include <errno.h>

char* getSequence(uint4 seqId)
{
	char* sequence;

    // Declare memory for the sequence
	sequence = (char*)global_malloc(sizeof(char) * (readdb_sequenceData[seqId].sequenceLength + 1));
    int ii;
    for(ii = 0; ii < readdb_sequenceData[seqId].sequenceLength; ii++)
    {
        sequence[ii] = encoding_getLetter(readdb_sequenceData[seqId].sequence[ii]);
        if(sequence[ii] == 'U')
        {
            fprintf(stderr, "Selenocysteine (U) at position %d replaced by X\n", ii);
            sequence[ii] = 'X';
        }
    }
    sequence[ii] = '\0';
    return sequence;
}

void print_sequence(int seqId)
{
    char *seqDes = descriptions_getDescription_mem(readdb_sequenceData[seqId].descriptionStart, readdb_sequenceData[seqId].descriptionLength);
    char *seq = getSequence(seqId); 
    printf(">%s\n%s\n", seqDes, seq);
    free(seqDes); free(seq);
}

void fprint_sequence(int seqId, int count)
{
    FILE *queryfile;
    char queryname[40];
    sprintf(queryname, "query_%d", count);
    queryfile = fopen(queryname, "w");
    char *seqDes = descriptions_getDescription_mem(readdb_sequenceData[seqId].descriptionStart, readdb_sequenceData[seqId].descriptionLength);
    char *seq = getSequence(seqId); 
    fprintf(queryfile, ">%s\n%s\n", seqDes, seq);
    free(seqDes); free(seq);
    fclose(queryfile);
}

int4 main(int4 argc, char* argv[])
{
    // User must provide FASTA format file at command line
	if (argc < 4)
	{
		fprintf(stderr, "Useage: indexdb <DB filename> <Query Length> <Number of Sequences>\n");
		exit(-1);
	}

	char *filename = argv[1];
    int queryLen = atoi(argv[2]);
	int numSeqs = atoi(argv[3]);

    // Open sequence data file and read information
	encoding_initialize(encoding_protein);
    readdb_open(filename);

    int ii;
    //int interval = readdb_numberOfSequences / numSeqs; 

    int queryCnt = 0;
    for(ii = 0; ii < readdb_numVolumeSequences && queryCnt < numSeqs; ii++)
    {
        if(readdb_sequenceData[ii].sequenceLength >= queryLen && 
                readdb_sequenceData[ii].sequenceLength < queryLen + (readdb_sequenceData[ii].sequenceLength * 0.05))
        {
            print_sequence(ii);
            //fprint_sequence(ii, queryCnt + 1);
            queryCnt++;
        }
    }

    fprintf(stderr, "maxQueryLength: %d\n", readdb_sequenceData[ii].sequenceLength);

    //close database
    readdb_close();
    return 0;
}

