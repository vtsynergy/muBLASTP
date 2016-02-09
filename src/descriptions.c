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

FILE *descriptions_file_multi[BATCH_SIZE];

char *descriptions_file_mem;
int descriptions_size;

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

void descriptions_open_multi(char *filename) {
  int ii;
  for (ii = 0; ii < BATCH_SIZE; ii++) {
    descriptions_file = fopen(filename, "r");

    if (descriptions_file == NULL) {
      fprintf(stderr, "%s\n", strerror(errno));
      fprintf(stderr, "Error opening file %s for reading\n", filename);
      exit(-1);
    }

    descriptions_file_multi[ii] = descriptions_file;
  }

  descriptions_filename = filename;
}

void descriptions_open_load(char *filename) {
    int ii;
    descriptions_file = fopen(filename, "r");

    if (descriptions_file == NULL) {
        fprintf(stderr, "%s\n", strerror(errno));
        fprintf(stderr, "Error opening file %s for reading\n", filename);
        exit(-1);
    }

    fseek(descriptions_file, 0L, SEEK_END);
    descriptions_size = ftell(descriptions_file);

    descriptions_file_mem = (char *)global_malloc(descriptions_size);

    rewind(descriptions_file);

    fread(descriptions_file_mem,1,descriptions_size,descriptions_file);

    //descriptions_file_multi[ii] = descriptions_file;
    //descriptions_filename = filename;
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

  //    printf("(%d,%d)%s\n", descriptionLocation, descriptionLength,
  // description);

  return description;
}

// Close the file
void descriptions_close() { fclose(descriptions_file); }
void descriptions_close_free() { fclose(descriptions_file); free(descriptions_file_mem); }

char *descriptions_getDescription_mem(uint4 descriptionLocation,
                                        uint4 descriptionLength) {
  char *description;

  // Declare memory for the description
  //description = (char *)global_malloc(sizeof(char) * (descriptionLength + 1));

  //fseek(descriptions_file_multi[queryNum], descriptionLocation, SEEK_SET);

  // Read the description into the new buffer
  //if (fgets(description, descriptionLength + 1,
            //descriptions_file_multi[queryNum]) == NULL) {
    //fprintf(stderr, "%s\n", strerror(errno));
    //fprintf(stderr, "Error reading from file %s\n", descriptions_filename);
    //exit(-1);
  //}

  description = descriptions_file_mem + descriptionLocation; 

  ASSERT(descriptionLocation + descriptionLength < descriptions_size)
  //    printf("(%d,%d)%s\n", descriptionLocation, descriptionLength,
  // description);

  return description;
}

// Get the description located at the given position in the file
char *descriptions_getDescription_multi(uint4 descriptionLocation,
                                        uint4 descriptionLength, int queryNum) {
  char *description;

  // Declare memory for the description
  description = (char *)global_malloc(sizeof(char) * (descriptionLength + 1));

  fseek(descriptions_file_multi[queryNum], descriptionLocation, SEEK_SET);

  // Read the description into the new buffer
  if (fgets(description, descriptionLength + 1,
            descriptions_file_multi[queryNum]) == NULL) {
    fprintf(stderr, "%s\n", strerror(errno));
    fprintf(stderr, "Error reading from file %s\n", descriptions_filename);
    exit(-1);
  }

  //    printf("(%d,%d)%s\n", descriptionLocation, descriptionLength,
  // description);

  return description;
}
