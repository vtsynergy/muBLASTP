#ifndef _readNcbidb_
#define _readNcbidb_

extern int readNcbidb_version, readNcbidb_alphabetType;
extern uint4 readNcbidb_numberOfSequences, readNcbidb_longestSequenceLength;
extern uint8 readNcbidb_numberOfLetters;

// Open NCBI collection files for reading
void readNcbidb_open(char* filename);
// Read a sequence from the collection
unsigned char* readNcbidb_getSequence(uint4 sequenceNumber, uint4* length);
// Read a description from the collection
unsigned char* readNcbidb_getDescription(uint4 sequenceNumber, uint4* length);

#endif
