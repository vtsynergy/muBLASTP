
#ifndef _ungappedExtension_
#define _ungappedExtension_

unsigned char* ungappedExtension_subjectEndReached;
int4 ungappedExtension_bestScore;
struct memBlocks* ungappedExtension_extensions;

// Scoring status constants
#define ungappedExtension_DELETED 0
#define ungappedExtension_UNGAPPED 1
#define ungappedExtension_SEMIGAPPED 2
#define ungappedExtension_GAPPED 3
#define ungappedExtension_JOINED 4

struct coordinate ungappedExtension_findProteinSeed(struct ungappedExtension* ungappedExtension,
                                                    struct PSSMatrix PSSMatrix, unsigned char* subject);


// Initialize the creation of ungapped extensions
void ungappedExtension_initialize();

// Perform an ungapped extension between point4s queryStart,subjectStart and queryEnd,subjectEnd
// and extend in each direction until score drops below best score yet minus a dropoff parameter
struct ungappedExtension* ungappedExtension_extend(int2** queryHit, unsigned char* subjectHit,
                        unsigned char* lastHit, struct PSSMatrix PSSMatrix, unsigned char* subject);

// Perform an ungapped extension when the seed is only a single hit on the diagonal, rather
// than a pair of hits.
struct ungappedExtension* ungappedExtension_oneHitExtend(int2** queryHit,
		unsigned char* subjectHit, struct PSSMatrix PSSMatrix, unsigned char* subject);

// Perform one-hit seeded ungapped extension for nucleotide, 1 packed-byte at a time
struct ungappedExtension* ungappedExtension_nucleotideExtend(int4 queryHitOffset,
	int4 subjectHitOffset, struct PSSMatrix PSSMatrix, unsigned char* subject,
    uint4 subjectLength);

// Find seed point4 for an ungapped extension
extern inline void ungappedExtension_findSeed(struct ungappedExtension* ungappedExtension,
		                                      struct PSSMatrix PSSMatrix, unsigned char* subject);

void ungappedExtension_print(struct ungappedExtension* extension);

#endif
