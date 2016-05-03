#ifndef _descriptions_
#define _descriptions_

// Open text file containing descriptions
void descriptions_open(char* filename);

void descriptions_open_mem(char* filename, uint8 offset, uint8 size);

// Get the description located at the given position in the file
char* descriptions_getDescription(uint4 descriptionLocation, uint4 descriptionLength);

char* descriptions_getDescription_mem(uint4 descriptionLocation, uint4 descriptionLength);

// Close the file
void descriptions_close();

void descriptions_close_mem();

#endif

