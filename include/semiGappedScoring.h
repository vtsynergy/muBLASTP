#ifndef _semiGappedScoring_
#define _semiGappedScoring_

// Perform semi-gapped alignment with restricted insertion
int4 semiGappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                        int4 subjectSize, unsigned char* subject, int4 dropoff);

void semiGappedScoring_free();

#endif
