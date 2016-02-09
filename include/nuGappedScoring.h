#ifndef _nuGappedScoring_
#define _nuGappedScoring_

// Perform gapped alinment with restricted insertion
int4 nuGappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                          int4 subjectSize, unsigned char* subject, int4 dropoff);

void nuGappedScoring_free();

#endif

