// readdb.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license
// agreement,
// provided that this statement is retained.
//
// Functions for reading a formatted database

#include "blast.h"
#include <sys/time.h>

uint4 readdb_numberOfSequences, readdb_longestSequenceLength,
      readdb_dbAlphabetType;
uint8 readdb_numberOfLetters;
unsigned char *readdb_filename, *readdb_data, *readdb_sequences;
uint4 readdb_fileSize, readdb_sequenceCount;
char *readdb_sequenceFilename, *readdb_descriptionsFilename,
     *readdb_dataFilename;
struct readFile readdb_readSequences, readdb_readData;
struct readFile_mem readdb_readSequences_mem, readdb_readData_mem;
uint4 readdb_volumeNumber, readdb_numberOfClusters, readdb_numberOfVolumes;
struct child *readdb_childBuffer = NULL;
uint4 readdb_sizeChildBuffer = 0, readdb_volume, readdb_numVolumeSequences, readdb_volumeOffset = 0;
uint8 readdb_numVolumeLetters, readdb_descriptionStart;
struct sequenceData *readdb_sequenceData;

// Open formatted collection for reading
void readdb_open(char *filename) {
    uint4 databaseVersion, encodedLength, sequenceLength, descriptionLength,
          sequenceCount, offset, oid;
    char *wildcardsFile;

    readdb_filename = filename;
    readdb_sequenceCount = 0;
    readdb_descriptionStart = 0;
    readdb_volumeNumber = 0;
    readdb_volume = 0;

    readdb_numVolumeLetters = 0;

    // Open sequence file for reading, mapping contents to readdb_address
    readdb_sequenceFilename = (char *)global_malloc(strlen(filename) + 15);
    sprintf(readdb_sequenceFilename, "%s.sequences", filename);

    // Report error if .sequences file doesn't exist
    if (!readFile_checkOpen(readdb_sequenceFilename)) {
        fprintf(stderr, "Error: unable to open file %s for reading\n",
                readdb_sequenceFilename);
        fprintf(stderr, "Before searching a collection you must to format it ");
        fprintf(stderr, "using FSA-BLAST's formatdb tool which creates .sequences "
                ".descriptions and ");
        fprintf(stderr, ".data files for the collection.\n");
        exit(-1);
    }

    readdb_readSequences = readFile_open(readdb_sequenceFilename);
    readdb_sequences = (char *)readdb_readSequences.address;

    // Open descriptions file for reading
    readdb_descriptionsFilename = (char *)global_malloc(strlen(filename) + 15);
    sprintf(readdb_descriptionsFilename, "%s.descriptions", filename);


    descriptions_open(readdb_descriptionsFilename);

    // Open data file for reading
    readdb_dataFilename = (char *)global_malloc(strlen(filename) + 13);
    sprintf(readdb_dataFilename, "%s.data", filename);
    readdb_readData = readFile_open(readdb_dataFilename);
    readdb_data = (char *)readdb_readData.address;

    // Read in a set of wildcards
    wildcardsFile = (char *)global_malloc(strlen(filename) + 13);
    sprintf(wildcardsFile, "%s.wildcards", filename);
    wildcards_readWildcards(wildcardsFile);
    free(wildcardsFile);

    // Get the sequences file size
    readdb_fileSize = readdb_readSequences.fileSize;

    // Read database statistics
    vbyte_getVbyte(readdb_data, &databaseVersion);
    if (databaseVersion != constants_databaseVersion) {
        fprintf(stderr, "Error: Invalid formatted database version %d. ",
                databaseVersion);
        fprintf(stderr, "Current supported version is %d.\n",
                constants_databaseVersion);
        fprintf(stderr, "Use formatdb tool to reformat database to version %d.\n",
                constants_databaseVersion);
        exit(-1);
    }
    vbyte_getVbyte(readdb_data, &readdb_numberOfSequences);
    readdb_data = vbyte_get64vbyte(readdb_data, &readdb_numberOfLetters);
    vbyte_getVbyte(readdb_data, &readdb_longestSequenceLength);
    vbyte_getVbyte(readdb_data, &readdb_dbAlphabetType);
    vbyte_getVbyte(readdb_data, &readdb_numberOfClusters);
    vbyte_getVbyte(readdb_data, &readdb_numberOfVolumes);

    readdb_sequenceData = (struct sequenceData *)global_malloc(
            sizeof(struct sequenceData) * readdb_numberOfClusters);

    fprintf(stderr, "readdb_numberOfSequences: %d readdb_numberOfVolumes: %d\n", 
            readdb_numberOfSequences, readdb_numberOfVolumes);

    // For each sequence in first volume
    offset = 1;
    sequenceCount = 0;
    while (sequenceCount < readdb_numberOfClusters && offset < readdb_fileSize) {
        // Read sequence data
        vbyte_getVbyte(readdb_data, &descriptionLength);
        vbyte_getVbyte(readdb_data, &sequenceLength);
        vbyte_getVbyte(readdb_data, &encodedLength);
        vbyte_getVbyte(readdb_data, &oid);

        readdb_sequenceData[sequenceCount].descriptionLength = descriptionLength;
        readdb_sequenceData[sequenceCount].descriptionStart =
            readdb_descriptionStart;
        readdb_sequenceData[sequenceCount].sequenceLength = sequenceLength;
        readdb_sequenceData[sequenceCount].encodedLength = encodedLength;
        readdb_sequenceData[sequenceCount].oid = oid;

        readdb_numVolumeLetters += sequenceLength;

        // Record pointer to sequence
        readdb_sequenceData[sequenceCount].sequence = readdb_sequences + offset;

        // If protein data skip past sentinal byte
        if (readdb_dbAlphabetType == encoding_protein)
            readdb_sequenceData[sequenceCount].sequence++;

        offset += encodedLength;

        readdb_descriptionStart += descriptionLength;

        sequenceCount++;
    }

    readdb_numVolumeSequences = sequenceCount;
}


// Open formatted collection for reading
void readdb_open_mem(char *filename) {

    fprintf(stderr, "opening database...");
    struct timeval start, end;
    gettimeofday(&start, NULL);

    uint4 databaseVersion, encodedLength, sequenceLength, descriptionLength,
          sequenceCount, offset;
    char *wildcardsFile;

    readdb_filename = filename;
    readdb_sequenceCount = 0;
    readdb_descriptionStart = 0;
    readdb_volumeNumber = 0;
    readdb_volume = 0;

    readdb_numVolumeLetters = 0;

    // Open sequence file for reading, mapping contents to readdb_address
    readdb_sequenceFilename = (char *)global_malloc(strlen(filename) + 15);
    sprintf(readdb_sequenceFilename, "%s.sequences", filename);

    // Report error if .sequences file doesn't exist
    if (!readFile_checkOpen(readdb_sequenceFilename)) {
        fprintf(stderr, "Error: unable to open file %s for reading\n",
                readdb_sequenceFilename);
        fprintf(stderr, "Before searching a collection you must to format it ");
        fprintf(stderr, "using FSA-BLAST's formatdb tool which creates .sequences "
                ".descriptions and ");
        fprintf(stderr, ".data files for the collection.\n");
        exit(-1);
    }

    readdb_readSequences_mem = readFile_open_mem(readdb_sequenceFilename);
    readdb_sequences = (char *)readdb_readSequences_mem.address;

    // Open descriptions file for reading
    readdb_descriptionsFilename = (char *)global_malloc(strlen(filename) + 15);
    sprintf(readdb_descriptionsFilename, "%s.descriptions", filename);



    // Open data file for reading
    readdb_dataFilename = (char *)global_malloc(strlen(filename) + 13);
    sprintf(readdb_dataFilename, "%s.data", filename);
    readdb_readData_mem = readFile_open_mem(readdb_dataFilename);
    readdb_data = (char *)readdb_readData_mem.address;

    // Read in a set of wildcards
    wildcardsFile = (char *)global_malloc(strlen(filename) + 13);
    sprintf(wildcardsFile, "%s.wildcards", filename);
    wildcards_readWildcards(wildcardsFile);
    free(wildcardsFile);

    // Get the sequences file size
    readdb_fileSize = readdb_readSequences_mem.fileSize;

    // Read database statistics
    vbyte_getVbyte(readdb_data, &databaseVersion);
    if (databaseVersion != constants_databaseVersion) {
        fprintf(stderr, "Error: Invalid formatted database version %d. ",
                databaseVersion);
        fprintf(stderr, "Current supported version is %d.\n",
                constants_databaseVersion);
        fprintf(stderr, "Use formatdb tool to reformat database to version %d.\n",
                constants_databaseVersion);
        exit(-1);
    }
    vbyte_getVbyte(readdb_data, &readdb_numberOfSequences);
    readdb_data = vbyte_get64vbyte(readdb_data, &readdb_numberOfLetters);
    vbyte_getVbyte(readdb_data, &readdb_longestSequenceLength);
    vbyte_getVbyte(readdb_data, &readdb_dbAlphabetType);
    vbyte_getVbyte(readdb_data, &readdb_numberOfClusters);
    vbyte_getVbyte(readdb_data, &readdb_numberOfVolumes);

    fprintf(stderr, "numOfVolumes: %d "
            "totalNumberSequences: %d "
            "readdb_numberOfLetters: %lu\n",
            readdb_numberOfVolumes,
            readdb_numberOfClusters,
            readdb_numberOfLetters
            );

    fprintf(stderr, "loading database(%d/%d)...", readdb_volume + 1, readdb_numberOfVolumes);

    readdb_sequenceData = (struct sequenceData *)global_malloc(
            sizeof(struct sequenceData) * readdb_numberOfClusters);

    // For each sequence in first volume
    offset = 1;
    sequenceCount = 0;
    uint8 description_offset = 0;
    uint4 oid;
    while (sequenceCount < readdb_numberOfClusters && offset < readdb_fileSize) {
        // Read sequence data
        vbyte_getVbyte(readdb_data, &descriptionLength);
        vbyte_getVbyte(readdb_data, &sequenceLength);
        vbyte_getVbyte(readdb_data, &encodedLength);
        vbyte_getVbyte(readdb_data, &oid);

        readdb_sequenceData[sequenceCount].descriptionLength = descriptionLength;
        readdb_sequenceData[sequenceCount].descriptionStart =
            readdb_descriptionStart;
        readdb_sequenceData[sequenceCount].sequenceLength = sequenceLength;
        readdb_sequenceData[sequenceCount].encodedLength = encodedLength;
        readdb_sequenceData[sequenceCount].oid = oid;

        readdb_numVolumeLetters += sequenceLength;

        // Record pointer to sequence
        readdb_sequenceData[sequenceCount].sequence = readdb_sequences + offset;

        // If protein data skip past sentinal byte
        if (readdb_dbAlphabetType == encoding_protein)
            readdb_sequenceData[sequenceCount].sequence++;

        offset += encodedLength;

        readdb_descriptionStart += descriptionLength;

        sequenceCount++;
    }


    uint4 description_size = readdb_descriptionStart; 
    descriptions_open_mem(readdb_descriptionsFilename, 
            description_offset, 
            description_size);

    readdb_numVolumeSequences = sequenceCount;

    gettimeofday(&end, NULL);
    long readdb_time = ((end.tv_sec * 1000000 + end.tv_usec) -
            (start.tv_sec * 1000000 + start.tv_usec));

    fprintf(stderr, "time: %f numOfSequencesVolume: %d\n", 
            (float)readdb_time * 1e-6,
            readdb_numVolumeSequences
            );

    //fprintf(stderr, "readdb_numberOfSequences: %d readdb_numberOfVolumes: %d\n", 
    //readdb_numberOfSequences, readdb_numberOfVolumes);
}

// Load the next volume
int readdb_nextVolume() {
    uint4 encodedLength, sequenceLength, descriptionLength, sequenceCount = 0,
          offset = 0;


    readdb_volumeOffset += readdb_numVolumeSequences;

    readdb_numVolumeLetters = 0;

    // Return 0 if no more volumes to read
    readdb_volume++;
    if (readdb_volume >= readdb_numberOfVolumes)
        return 0;

    // Close current volume
    readFile_close(readdb_readSequences);

    // Open next volume
    sprintf(readdb_sequenceFilename, "%s.sequences%d", readdb_filename,
            readdb_volume);
    readdb_readSequences = readFile_open(readdb_sequenceFilename);
    readdb_sequences = (unsigned char *)readdb_readSequences.address;

    // Get the sequences file size
    readdb_fileSize = readdb_readSequences.fileSize;

    // For each sequence in next volume volume
    offset = 1;
    sequenceCount = 0;
    uint4 oid;
    while (sequenceCount < readdb_numberOfClusters && offset < readdb_fileSize) {
        // Read sequence data
        vbyte_getVbyte(readdb_data, &descriptionLength);
        vbyte_getVbyte(readdb_data, &sequenceLength);
        vbyte_getVbyte(readdb_data, &encodedLength);
        vbyte_getVbyte(readdb_data, &oid);

        readdb_sequenceData[readdb_volumeOffset + sequenceCount].descriptionLength = descriptionLength;
        readdb_sequenceData[readdb_volumeOffset + sequenceCount].descriptionStart =
            readdb_descriptionStart;
        readdb_sequenceData[readdb_volumeOffset + sequenceCount].sequenceLength = sequenceLength;
        readdb_sequenceData[readdb_volumeOffset + sequenceCount].encodedLength = encodedLength;
        readdb_sequenceData[readdb_volumeOffset + sequenceCount].oid = oid;

        // Record pointer to sequence
        readdb_sequenceData[readdb_volumeOffset + sequenceCount].sequence = readdb_sequences + offset;

        // If protein data skip past sentinal byte
        if (encoding_alphabetType == encoding_protein)
            readdb_sequenceData[readdb_volumeOffset + sequenceCount].sequence++;

        readdb_numVolumeLetters += sequenceLength;

        offset += encodedLength;

        readdb_descriptionStart += descriptionLength;

        sequenceCount++;
    }

    readdb_sequenceCount = 0;
    readdb_numVolumeSequences = sequenceCount;

    return 1;
}

int readdb_nextVolume_mem() {


    struct timeval start, end;
    gettimeofday(&start, NULL);

    uint4 encodedLength, 
          sequenceLength, 
          descriptionLength, 
          sequenceCount = 0,
          offset = 0;

    readdb_volumeOffset += readdb_numVolumeSequences;

    readdb_numVolumeLetters = 0;

    // Return 0 if no more volumes to read
    readdb_volume++;
    if (readdb_volume >= readdb_numberOfVolumes)
        return 0;

    fprintf(stderr, "loading database(%d/%d)...", readdb_volume + 1, readdb_numberOfVolumes);

    // Close current volume
    readFile_close_mem(readdb_readSequences_mem);
    descriptions_close_mem();

    // Open next volume
    sprintf(readdb_sequenceFilename, "%s.sequences%d", readdb_filename,
            readdb_volume);
    readdb_readSequences_mem = readFile_open_mem(readdb_sequenceFilename);
    readdb_sequences = (unsigned char *)readdb_readSequences_mem.address;

    // Get the sequences file size
    readdb_fileSize = readdb_readSequences_mem.fileSize;

    // For each sequence in next volume volume
    offset = 1;
    sequenceCount = 0;
    uint8 description_offset = readdb_descriptionStart;
    readdb_descriptionStart = 0;
    uint4 oid;
    while (sequenceCount < readdb_numberOfClusters && offset < readdb_fileSize) {
        // Read sequence data
        vbyte_getVbyte(readdb_data, &descriptionLength);
        vbyte_getVbyte(readdb_data, &sequenceLength);
        vbyte_getVbyte(readdb_data, &encodedLength);
        vbyte_getVbyte(readdb_data, &oid);

        readdb_sequenceData[readdb_volumeOffset + sequenceCount].descriptionLength = descriptionLength;
        readdb_sequenceData[readdb_volumeOffset + sequenceCount].descriptionStart =
            readdb_descriptionStart;
        readdb_sequenceData[readdb_volumeOffset + sequenceCount].sequenceLength = sequenceLength;
        readdb_sequenceData[readdb_volumeOffset + sequenceCount].encodedLength = encodedLength;
        readdb_sequenceData[readdb_volumeOffset + sequenceCount].oid = oid;

        // Record pointer to sequence
        readdb_sequenceData[readdb_volumeOffset + sequenceCount].sequence = readdb_sequences + offset - 1;

        // If protein data skip past sentinal byte
        if (encoding_alphabetType == encoding_protein)
            readdb_sequenceData[readdb_volumeOffset + sequenceCount].sequence++;

        readdb_numVolumeLetters += sequenceLength;

        offset += encodedLength;

        readdb_descriptionStart += descriptionLength;

        sequenceCount++;
    }

    uint8 description_size = readdb_descriptionStart; 

    //fprintf(stderr, "seqStart: %d seqSize: %d descripStart: %lu descripSize: %lu\n", 
    //readdb_volumeOffset,  readdb_fileSize, description_offset,  description_size);


    descriptions_open_mem(readdb_descriptionsFilename, 
            description_offset, 
            description_size);

    readdb_descriptionStart = description_offset + description_size;

    readdb_sequenceCount = 0;
    readdb_numVolumeSequences = sequenceCount;

    gettimeofday(&end, NULL);

    long readdb_time = ((end.tv_sec * 1000000 + end.tv_usec) -
            (start.tv_sec * 1000000 + start.tv_usec));

    fprintf(stderr, "time: %f numOfSequencesVolume: %d\n", 
            (float)readdb_time * 1e-6,
            readdb_numVolumeSequences
            );

    return 1;
}

// Read a sequence and description information. Return 0 if end-of-collection.
int readdb_readSequence(unsigned char **sequence, uint4 *sequenceLength,
        uint4 *descriptionStart, uint4 *descriptionLength,
        uint4 *encodedLength) {
    if (readdb_sequenceCount < readdb_numVolumeSequences) {
        *descriptionStart =
            readdb_sequenceData[readdb_sequenceCount].descriptionStart;
        *descriptionLength =
            readdb_sequenceData[readdb_sequenceCount].descriptionLength;
        *sequenceLength = readdb_sequenceData[readdb_sequenceCount].sequenceLength;
        *encodedLength = readdb_sequenceData[readdb_sequenceCount].encodedLength;
        *sequence = readdb_sequenceData[readdb_sequenceCount].sequence;

        readdb_sequenceCount++;
        return 1;
    } else {
        return 0;
    }
}

// Get the children
struct child *readdb_getChildren(unsigned char *sequence, uint4 sequenceLength,
        uint4 encodedLength, uint4 descriptionLocation,
        uint4 *numChildren) {
    unsigned char *edits, *editsEnd;
    struct child *child, *children;
    uint4 editNum, position;

    editsEnd = sequence + encodedLength - 1;

    // Advance to start of child information
    edits = sequence + sequenceLength + 1;

    *numChildren = 0;
    while (edits < editsEnd) {
        // Increase size of children buffer if required
        if (*numChildren >= readdb_sizeChildBuffer) {
            readdb_sizeChildBuffer = (readdb_sizeChildBuffer + 1) * 2;
            readdb_childBuffer = (struct child *)global_realloc(
                    readdb_childBuffer, sizeof(struct child) * readdb_sizeChildBuffer);
        }

        child = readdb_childBuffer + *numChildren;

        // Read child details
        vbyte_getVbyte(edits, &(child->descriptionLength));
        vbyte_getVbyte(edits, &(child->regionStart));
        vbyte_getVbyte(edits, &(child->length));
        vbyte_getVbyte(edits, &(child->numEdits));

        // Start with copy of parent sequence
        child->sequence = (unsigned char *)global_malloc(child->length + 2);
        child->sequence++;
        memcpy(child->sequence, sequence + child->regionStart, child->length);

        // Add sentinal codes to either end
        child->sequence[-1] = encoding_sentinalCode;
        child->sequence[child->length] = encoding_sentinalCode;

        child->edits =
            (struct edit *)global_malloc(sizeof(struct edit) * child->numEdits);

        // Read edits and update sequence
        position = 0;
        editNum = 0;
        while (editNum < child->numEdits) {
            // Locate the next wildcard in the child
            while (child->sequence[position] < encoding_aaStartWildcards)
                position++;

            // Read position and character
            child->edits[editNum].position = position;
            child->edits[editNum].code = *edits;
            edits++;

            // Change character in child
            child->sequence[position] = child->edits[editNum].code;

            editNum++;
        }

        child->descriptionLocation = descriptionLocation;
        descriptionLocation += child->descriptionLength;
        (*numChildren)++;
    }

    // Copy children from buffer into new memory block and return
    children = (struct child *)global_malloc(sizeof(struct child) * *numChildren);
    memcpy(children, readdb_childBuffer, sizeof(struct child) * *numChildren);

    return children;
}

// Close the database for reading
void readdb_close() {
    free(readdb_sequenceData);
    free(readdb_sequenceFilename);
    free(readdb_descriptionsFilename);
    free(readdb_dataFilename);
    readFile_close(readdb_readSequences);
    readFile_close(readdb_readData);
    descriptions_close();
    free(readdb_childBuffer);
    readdb_childBuffer = NULL;
    readdb_sizeChildBuffer = 0;
}

void readdb_close_mem() {
    free(readdb_sequenceData);
    free(readdb_sequenceFilename);
    free(readdb_descriptionsFilename);
    free(readdb_dataFilename);
    readFile_close_mem(readdb_readSequences_mem);
    readFile_close_mem(readdb_readData_mem);
    descriptions_close_mem();

    free(readdb_childBuffer);
    readdb_childBuffer = NULL;
    readdb_sizeChildBuffer = 0;
}
