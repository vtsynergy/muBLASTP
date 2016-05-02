#ifndef _alignments_
#define _alignments_

// All alignments with an ungapped extension above gapping threshold
struct memBlocks* alignments_alignments;
// Good alignments which contain a semi-gapped extension above cutoff
struct memSingleBlock* alignments_goodAlignments;
// Final alignments which contain a gapped extension above cutoff
struct memSingleBlock* alignments_finalAlignments;
struct alignment* alignments_currentAlignment;
// Number of clusters with high-scoring alignments
extern uint4 alignments_numClusters;

// Initialize array storing point4ers to alignments
void alignments_initialize();

// Use the next available slot in the alignment object array
void alignments_createNew(uint4 descriptionLocation, uint4 descriptionLength,
                          unsigned char* subject, int4 subjectLength,
                          int4 encodedLength);

// Add an ungapped extension to an alignment's list of ungapped extensions, which
// are ordered highest nominal score first, lowest score last
void alignments_addUngappedExtension(struct ungappedExtension* newExtension);

// Perform the gapped-alignment scoring stage of BLAST
void alignments_performScoring(struct PSSMatrix PSSMatrix);

// Get the subject sequence descriptions for all of the final alignments
void alignments_getFinalAlignmentDescriptions();

// Copy subject sequences associated with good alignments into memory before
// the volume is closed
void alignments_loadSubjectsIntoMemory(struct PSSMatrix PSSMatrix);

// Unpack the subject
void alignment_unpackSubject(struct alignment* alignment);

// Perform initial scoring of all ungapped extensions to find "good" alignments that may
// score above the cutoff
void alignments_findGoodAlignments(struct PSSMatrix PSSMatrix);

// Given a collection of good alignments, find the final top N alignments above cutoff
void alignments_findFinalAlignments(struct PSSMatrix PSSMatrix);

// Get the tracebacks for all of the final alignments
void alignments_getTracebacks(struct PSSMatrix PSSMatrix);

// Returns true if the given alignment score is high enough for it to be
// one of the N displayed final alignments
int alignments_isFinalAlignment(uint4 score);

// Add a high-scoring gapped extension to this alignment
void alignments_addGappedExtension(struct alignment* alignment, struct gappedExtension* newExtension);

// Add the current alignment (which contains at least one gapped extension
// scoring above cutoff) to to-be-sorted list of final alignments
void alignments_addFinalAlignment(int4 highestNominalScore, struct alignment* alignment);

// Sort the array of final alignments in order of score
void alignments_sortFinalAlignments();

// Free memory used to store alignments
void alignments_free();


// Compare the two alignments' highest scores. Return -1 if alignment1 < alignment2
// 1 if alignment1 > alignment2 and 0 if they are equal
int4 alignments_compareFinalAlignments(const void* alignment1,
                                      const void* alignment2);

int4 alignments_compareAlignmentDescriptionLocations(const void* alignment1, const void* alignment2);


unsigned char* alignments_selectRegion(struct alignment* alignment,
               struct ungappedExtension* ungappedExtension);

struct trace fasterGappedExtension_traceBeforeSeed(struct dpResults beforeDpResults, struct coordinate seed);
struct trace fasterGappedExtension_traceAfterSeed(struct dpResults beforeDpResults, int4 queryLength);

struct trace fasterGappedExtension_joinTraces(struct trace beforeTrace, struct trace afterTrace);


void alignment_freeAlignment(struct alignment *alignment);


int alignments_contains(struct ungappedExtension *extension1,
                        struct ungappedExtension *extension2);

void alignments_pruneOverlappingExtensions(struct alignment *alignment);


int4 alignments_compareFinalAlignments2(const void *alignment1,
                                       const void *alignment2);

void
alignments_addUngappedExtensionAtEnd(struct alignment *alignment,
                                     struct ungappedExtension *newExtension);
#endif
