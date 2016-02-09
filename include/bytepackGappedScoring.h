#ifndef _bytepackGappedScoring_
#define _bytepackGappedScoring_

// Perform gapped alignment using byte-packed alignment technique
int4 bytepackGappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                                int4 subjectSize, unsigned char* packedSubject, int4 dropoff);

void bytepackGappedScoring_free();

#endif

