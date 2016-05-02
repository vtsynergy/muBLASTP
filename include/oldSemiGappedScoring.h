#ifndef _oldSemiGappedScoring_
#define _oldSemiGappedScoring_

// Perform semi-gapped alignment without restricted insertion
int4 oldSemiGappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                        int4 subjectSize, unsigned char* subject, int4 dropoff);

void oldSemiGappedScoring_free();

#endif


