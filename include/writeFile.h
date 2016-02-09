#ifndef _writeFile_
#define _writeFile_

// On file for reading AND writing and map to memory, then return mapped memory address
void* writeFile_open(char* filename)

// Get file size/length
int4 writeFile_getSize()

// Unmap then close the file
void writeFile_close()

#endif
