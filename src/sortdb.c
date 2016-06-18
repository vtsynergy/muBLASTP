// blast.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Main code for blast
#include "blast.h"
#include <stdio.h>
#include <errno.h>
#include <unistd.h>

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

    ASSERT(readdb_sequenceData[seqId].sequenceLength > 0);
    char *sequenceBuffer = (char *)global_malloc(
            sizeof(char) * (readdb_sequenceData[seqId].sequenceLength + 1));
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
    char *seqDes = 
        descriptions_getDescription_mem(
                readdb_sequenceData[seqId].descriptionStart, 
                readdb_sequenceData[seqId].descriptionLength);

    char *seq = getSequence(seqId); 

    fprintf(output_file, ">%s\n%s\n", seqDes, seq);
    //fprintf(output_file, ">gi|%d\n%s\n", seqId, seq);

    free(seq);
    free(seqDes); 
}

int4 main(int4 argc, char* argv[])
{
    char *ifilename = NULL, *ofilename = NULL;
    int c;
    while((c = getopt(argc, argv, "i:o:")) != -1)
    {
        switch (c)
        {
            case 'i':
                ifilename = optarg;
                break;
            case 'o':
                ofilename = optarg;
                break;
            default:
                fprintf(stderr, "Useage: sortdb -i <Database> -o <Sorted database>\n");
                exit(-1);
        }
    }


    if(ifilename == NULL || ofilename == NULL)
    {
        fprintf(stderr, "Useage: sortdb -i <Database> -o <Sorted database>\n");
        exit(-1);
    }
    // User must provide FASTA format file at command line
	//if (argc < 3)
	//{
		//fprintf(stderr, "Useage: sortdb <DB filename> <Output filename>\n");
		//exit(-1);
	//}

	//char *ifilename = argv[1];
	//char *ofilename = argv[2];

    // Open sequence data file and read information
	encoding_initialize(encoding_protein);
    readdb_open_mem(ifilename);


    FILE *output_file;

    //output_file = stdout;
    output_file = fopen(ofilename, "w");

    int *seqId = (int *)malloc(sizeof(int) * readdb_numberOfSequences);

    uint4 startSeq = 0, endSeq = readdb_numVolumeSequences;


    while(1)
    {

        fprintf(stderr, "sorting volumn %d ", readdb_volume);

        int ii;
        for(ii = startSeq; ii < endSeq; ii++)
        {
            seqId[ii] = ii;
        }

        qsort(seqId + startSeq, readdb_numVolumeSequences, sizeof(int), compareSeq);

        for(ii = startSeq; ii < endSeq; ii++)
        {
            if(!((ii - startSeq) % 10000))
            {
                fprintf(stderr, ".");
            }

            print_sequence(seqId[ii], output_file);
        }

        fprintf(stderr, "\n");

        if (readdb_volume + 1 < readdb_numberOfVolumes) {
            readdb_nextVolume_mem();
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

