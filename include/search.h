#ifndef _search_
#define _search_

// Search a protein database using 1-hit extension mode
void search_protein1hit(struct PSSMatrix PSSMatrix, struct sequenceData* sequenceData,
                       uint4 numSequences, uint4 tickFrequency);

// Search a protein database using 2-hit extension mode
void search_protein2hit(struct PSSMatrix PSSMatrix, struct sequenceData* sequenceData,
                       uint4 numSequences, uint4 tickFrequency);

void search_protein2hit_lookup(struct PSSMatrix PSSMatrix,
                        struct sequenceData *sequenceData, uint4 numSequences,
                        uint4 tickFrequency);
#endif
