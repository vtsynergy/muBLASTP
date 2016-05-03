#ifndef _readdb_
#define _readdb_

extern uint4 readdb_numberOfSequences, readdb_longestSequenceLength, readdb_dbAlphabetType;
extern uint8 readdb_numberOfLetters;
extern unsigned char *readdb_filename, *readdb_sequences;
extern uint4 readdb_fileSize, readdb_sequenceCount;
extern uint4 readdb_volumeNumber, readdb_numberOfClusters, readdb_numberOfVolumes;
extern uint4 readdb_numVolumeSequences, readdb_volume;
extern uint4 readdb_volumeOffset;
extern uint8 readdb_numVolumeLetters;
extern struct sequenceData* readdb_sequenceData;

// Open formatted database for reading
void readdb_open(char* filename);
void readdb_open_mem(char* filename);

// Read a sequence and description information. Return 0 if end-of-collection.
int readdb_readSequence(unsigned char** sequence, uint4* sequenceLength, uint4* descriptionStart,
                        uint4* descriptionLength, uint4* encodedLength);

// Load the next volume
int readdb_nextVolume();
int readdb_nextVolume_mem();

// Get the children
struct child* readdb_getChildren(unsigned char* sequence, uint4 sequenceLength, uint4 encodedLength,
                                 uint4 descriptionLocation, uint4* numChildren);

// Close the database for reading
void readdb_close();
void readdb_close_mem();

#endif
