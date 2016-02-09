// blast.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Main code for blast
#include "blast.h"
#include <stdio.h>
#include <errno.h>

int compareSeq(const void *a, const void *b)
{
    int seqId_a = *(int *)a;
    int seqId_b = *(int *)b;
    return readdb_sequenceData[seqId_a].sequenceLength > readdb_sequenceData[seqId_b].sequenceLength; 
}

void sortdb(int *seqId)
{
    qsort(seqId, readdb_numberOfSequences, sizeof(int), compareSeq);
}

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
    char *seqDes = descriptions_getDescription(readdb_sequenceData[seqId].descriptionStart, readdb_sequenceData[seqId].descriptionLength);
    char *seq = getSequence(seqId); 
    printf(">%s\n%s\n", seqDes, seq);
    free(seqDes); free(seq);
}

int4 main(int4 argc, char* argv[])
{
    // User must provide FASTA format file at command line
	if (argc < 2)
	{
		fprintf(stderr, "Useage: indexdb <DB filename>\n");
		exit(-1);
	}

	char *filename = argv[1];

    // Open sequence data file and read information
	encoding_initialize(encoding_protein);
    readdb_open(filename);

    int *seqId = (int *)malloc(sizeof(int) * readdb_numberOfSequences);
    int ii;
    for(ii = 0; ii < readdb_numberOfSequences; ii++)
    {
        seqId[ii] = ii;
    }

    sortdb(seqId);

    for(ii = 0; ii < readdb_numberOfSequences; ii++)
    {
        print_sequence(seqId[ii]);
    }

    free(seqId);

    //close database
    readdb_close();
    return 0;
}

