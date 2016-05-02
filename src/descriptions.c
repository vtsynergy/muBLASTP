// descriptions.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license
// agreement,
// provided that this statement is retained.
//
// Code for reading descriptions from FASTA format file

#include "blast.h"
#include <stdio.h>
#include <errno.h>

char *descriptions_filename;
FILE *descriptions_file;


struct readFile descriptions_file_mem;

// Open text file containing descriptions
void descriptions_open(char *filename) {
    int ii;
    descriptions_file = fopen(filename, "r");

    if (descriptions_file == NULL) {
        fprintf(stderr, "%s\n", strerror(errno));
        fprintf(stderr, "Error opening file %s for reading\n", filename);
        exit(-1);
    }

    descriptions_filename = filename;
}

void descriptions_open_mem(char *filename) {
    //
    descriptions_file_mem = readFile_open(filename);
}

// Get the description located at the given position in the file
char *descriptions_getDescription(uint4 descriptionLocation,
        uint4 descriptionLength) {

    char *description;

    // Declare memory for the description
    description = (char *)global_malloc(sizeof(char) * (descriptionLength + 1));

    fseek(descriptions_file, descriptionLocation, SEEK_SET);

    // Read the description into the new buffer
    if (fgets(description, descriptionLength + 1, descriptions_file) == NULL) {
        fprintf(stderr, "%s\n", strerror(errno));
        fprintf(stderr, "Error reading from file %s\n", descriptions_filename);
        exit(-1);
    }

    return description;
}

// Close the file
void descriptions_close() 
{ 
    fclose(descriptions_file); 
}

void descriptions_close_mem() 
{
    readFile_close(descriptions_file_mem);
}

char *descriptions_getDescription_mem(uint4 descriptionLocation,
        uint4 descriptionLength) {

    char *description;

    // Declare memory for the description
    description = (char *)global_malloc(sizeof(char) * (descriptionLength + 1));

    memcpy(description, descriptions_file_mem.address + descriptionLocation, descriptionLength);

    description[descriptionLength] = '\0';

    return description;
}
