#ifndef _descriptions_
#define _descriptions_

// Open text file containing descriptions
void descriptions_open(char* filename);
void descriptions_open_multi(char* filename);
void descriptions_open_load(char* filename);

// Get the description located at the given position in the file
char* descriptions_getDescription(uint4 descriptionLocation, uint4 descriptionLength);


char* descriptions_getDescription_multi(uint4 descriptionLocation, uint4 descriptionLength, int queryNum);
char* descriptions_getDescription_mem(uint4 descriptionLocation, uint4 descriptionLength);
// Close the file
void descriptions_close();
void descriptions_close_free();

extern char *descriptions_file_mem;

#endif

