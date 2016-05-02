#ifndef _wordLookupDFA_
#define _wordLookupDFA_

struct neighbour
{
	int4 codeword;
    int4 score;
    int4 position;
};

struct group
{
	unsigned char* nextWords;
    struct group* nextGroups;
};

extern unsigned char *wordLookupDFA;
extern unsigned char **wordLookupDFA_blockAddresses;

extern uint2* wordLookupDFA_additionalQueryPositions;
extern int4 wordLookupDFA_numAdditionalQueryPositions;
extern struct group *wordLookupDFA_groups;
extern int4 wordLookupDFA_numCodes, wordLookupDFA_wordLength, wordLookupDFA_numBlocks;

// Build the word-lookup structure
void wordLookupDFA_build(struct PSSMatrix PSSMatrix, int4 numCodes, int4 wordLength);

// Print4 the contents of the word lookup table
void wordLookupDFA_print();

// Free memory used by the word lookup table
void wordLookupDFA_free();



unsigned char* wordLookupDFA_getCodes(int4 codeword, int4 wordLength);
int4 wordLookupDFA_getCodeword(unsigned char* codes, int4 wordLength);


void wordLookupDFA_getNeighbours(struct PSSMatrix PSSMatrix, int4 queryPosition,
                                 int4* numNeighbours, struct neighbour* neighbours);

void wordLookupSM_getNeighbours(char* codes, 
        struct scoreMatrix scoreMatrix, int4 queryPosition,
        int4* numNeighbours, struct neighbour* neighbours);

#endif

