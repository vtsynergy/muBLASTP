#ifndef _readFile_
#define _readFile_

struct readFile
{
	size_t fileSize;
	int4 fileDescriptor;
	void* address;
};

struct readFile_mem
{
	size_t fileSize;
    FILE *fileDescriptor;
	void* address;
};

// On file for reading and map to memory, then return mapped memory address
struct readFile readFile_open(char* filename);

struct readFile_mem readFile_open_mem_offset(char *filename, int8 offset, int8 size);

// Unmap then close the file
void readFile_close(struct readFile readFile);

// Check file exists
int readFile_checkOpen(char* filename);


struct readFile_mem readFile_open_mem(char *filename);

void readFile_close_mem(struct readFile_mem readFile);

#endif
