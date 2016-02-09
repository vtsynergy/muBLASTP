#ifndef _tableGappedScoring_
#define _tableGappedScoring_

// Perform gapped alignment using byte-packed alignment technique
int4 tableGappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                             int4 subjectSize, unsigned char* packedSubject, int4 dropoff);

#endif

