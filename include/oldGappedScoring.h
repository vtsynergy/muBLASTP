#ifndef _oldGappedScoring_
#define _oldGappedScoring_

int4 oldGappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                        int4 subjectSize, unsigned char* subject, int4 dropoff);

void oldGappedScoring_free();

#endif

                        
