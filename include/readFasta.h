#ifndef _readFasta_
#define _readFasta_

char* readFasta_descriptionBuffer;
char* readFasta_sequenceBuffer;
uint4 readFasta_sequenceLength;
uint4 readFasta_descriptionLength;

// Open FASTA file for reading
void readFasta_open(char* filename);

// Read a description and sequence from file, return 0 if end of file, otherwise 1
int4 readFasta_readSequence();

// Close fasta file and free memory used to store sequence
void readFasta_close();

#endif
