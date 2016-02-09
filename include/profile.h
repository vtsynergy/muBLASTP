#ifndef _profile_
#define _profile_
typedef struct profile_info
{
    long long blast_pseudoCycle;
    long long blast_hitDetectCycle;
    long long blast_sortCycle;
    long long blast_ungappedExtCycle;
    long long blast_findGoodCycle;
    long long blast_findFinalCycle;
    long long blast_getFinalCycle;
    long long blast_getTrackCycle;
}Profile;

#endif
