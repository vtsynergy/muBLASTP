// readFile.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license
// agreement,
// provided that this statement is retained.
//
// Code for reading the contents of a file using mmap

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>

#include "blast.h"
#include "readFile.h"

// On file for reading and map to memory, then return mapped memory address
struct readFile readFile_open(char *filename) {
  struct stat fileStats;
  struct readFile readFile;

  // Open file for reading
  if ((readFile.fileDescriptor = open(filename, O_RDONLY)) == -1) {
    fprintf(stderr, "%s\n", strerror(errno));
    fprintf(stderr, "Error opening file %s for reading: 1\n", filename);
    exit(-1);
  }

  // Get length of file
  if (fstat(readFile.fileDescriptor, &fileStats) == -1) {
    fprintf(stderr, "%s\n", strerror(errno));
    fprintf(stderr, "Error opening file %s for reading: 2\n", filename);
    exit(-1);
  }

  readFile.fileSize = fileStats.st_size;

  // Map address to fileSize bytes of application address space
  readFile.address = mmap(0, readFile.fileSize, PROT_READ, MAP_SHARED,
                          readFile.fileDescriptor, 0);

  // Check for error in mapping
  if (readFile.address == (void *)MAP_FAILED) {
    fprintf(stderr, "%s\n", strerror(errno));
    fprintf(stderr, "Error opening file %s for reading: 3\n", filename);
    exit(-1);
  }

  return readFile;
}

// On file for reading and map to memory, then return mapped memory address
struct readFile_mem readFile_open_mem(char *filename) {
    //struct stat fileStats;
  struct readFile_mem readFile;

  // Open file for reading
  if ((readFile.fileDescriptor = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "%s\n", strerror(errno));
    fprintf(stderr, "Error opening file %s for reading: 1\n", filename);
    exit(-1);
  }

 
  fseek(readFile.fileDescriptor, 0L, SEEK_END);
  readFile.fileSize = ftell(readFile.fileDescriptor);
  rewind(readFile.fileDescriptor);

  readFile.address = global_malloc(readFile.fileSize);

  if(fread(readFile.address, 1, readFile.fileSize, readFile.fileDescriptor) != readFile.fileSize)
  {
      fprintf(stderr, "readFile_open faialed\n");
      exit(0);
  }

  return readFile;

}

// On file for reading and map to memory, then return mapped memory address
struct readFile_mem readFile_open_mem_offset(char *filename, int8 offset, int8 size) {
    //struct stat fileStats;
  struct readFile_mem readFile;

  // Open file for reading
  if ((readFile.fileDescriptor = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "%s\n", strerror(errno));
    fprintf(stderr, "Error opening file %s for reading: 1\n", filename);
    exit(-1);
  }

 
  fseek(readFile.fileDescriptor, 0L, SEEK_END);
  readFile.fileSize = ftell(readFile.fileDescriptor);
  //rewind(readFile.fileDescriptor);

  ASSERT(offset + size <= readFile.fileSize);

  fseek(readFile.fileDescriptor, offset, SEEK_SET);

  readFile.address = global_malloc(size);

  if(fread(readFile.address, 1, size, readFile.fileDescriptor) != size)
  {
      fprintf(stderr, "readFile_open faialed\n");
      exit(0);
  }

  return readFile;

}



// Check file exists
int readFile_checkOpen(char *filename) {
  FILE *file;

  if ((file = fopen(filename, "r")) != NULL) {
    fclose(file);
    return 1;
  } else
    return 0;
}

// Unmap then close the file
void readFile_close(struct readFile readFile) {
  if (munmap(readFile.address, readFile.fileSize) < 0) {
    fprintf(stderr, "%s\n", strerror(errno));
    fprintf(stderr, "Error unmapping file\n");
    exit(-1);
  }
  close(readFile.fileDescriptor);
}

void readFile_close_mem(struct readFile_mem readFile) {
    free(readFile.address);
    fclose(readFile.fileDescriptor);
}
