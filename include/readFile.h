#ifndef _readFile_
#define _readFile_

struct readFile
{
	int4 fileSize;
	int4 fileDescriptor;
	void* address;
};

// On file for reading and map to memory, then return mapped memory address
struct readFile readFile_open(char* filename);

// Unmap then close the file
void readFile_close(struct readFile readFile);

// Check file exists
int readFile_checkOpen(char* filename);

#endif
