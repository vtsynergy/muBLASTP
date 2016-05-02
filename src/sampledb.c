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
    //printf(">gi|%d\n%s\n", seqId, seq);
    printf(">%s\n%s\n", seqDes, seq);
    free(seqDes); free(seq);
}

int4 main(int4 argc, char* argv[])
{
    // User must provide FASTA format file at command line
	if (argc < 3)
	{
		fprintf(stderr, "Useage: indexdb <DB filename> <Number of Sequences>\n");
		exit(-1);
	}

	char *filename = argv[1];
	int numSeqs = atoi(argv[2]);

    // Open sequence data file and read information
	encoding_initialize(encoding_protein);
    readdb_open(filename);

    int ii;
    int interval = readdb_numVolumeSequences / numSeqs; 

    int target_length = 931, target_id = 0;

    for(ii = 0; ii < readdb_numVolumeSequences; ii++)
    {
        if(!(ii % interval))
            print_sequence(ii);

        if(readdb_sequenceData[ii].sequenceLength == target_length)
        {
            target_id = ii;
        }
    }

    print_sequence(target_id);

    //close database
    readdb_close();
    return 0;
}

