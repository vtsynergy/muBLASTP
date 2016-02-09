/* ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 */




//#include "ncbi_stat.h"
#include <limits.h>
#include <math.h>



#include "blast.h"

#define HAVE_ERF


//#ifndef MAX
///** returns larger of a and b. */
//#define MAX(a, b) ((a) >= (b) ? (a) : (b))
//#endif

//#ifndef MIN
///** returns smaller of a and b. */
//#define MIN(a, b) ((a) > (b) ? (b) : (a))
//#endif

typedef struct BLAST_LetterProb {
    char ch;  /**< residue */
    double p; /**< probability of residue. */
} BLAST_LetterProb;

#define BLASTAA_SEQ_CODE 11
//#define INT2_MIN SHRT_MIN
//#define INT2_MAX SHRT_MAX

#define BLAST_SCORE_MIN INT2_MIN
#define BLAST_SCORE_MAX INT2_MAX

#define BLAST_SCORE_RANGE_MAX (BLAST_SCORE_MAX - BLAST_SCORE_MIN)

#define BLAST_KARLIN_K_SUMLIMIT_DEFAULT                                        \
    0.0001 /**< K_SUMLIMIT_DEFAULT == sumlimit used in BlastKarlinLHtoK() */

#define BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT                                   \
    (1.e-5) /**< LAMBDA_ACCURACY_DEFAULT == accuracy to which Lambda should be   \
              calc'd */

#define BLAST_KARLIN_LAMBDA_ITER_DEFAULT                                       \
    17 /**< LAMBDA_ITER_DEFAULT == no. of iterations in LambdaBis =              \
         ln(accuracy)/ln(2)*/

#define BLAST_KARLIN_LAMBDA0_DEFAULT                                           \
    0.5 /**< Initial guess for the value of Lambda in BlastKarlinLambdaNR */

#define BLAST_KARLIN_K_ITER_MAX                                                \
    100 /**< upper limit on iterations for BlastKarlinLHtoK */

/** Number of statistical parameters in each row of the precomputed tables. */
#define BLAST_NUM_STAT_VALUES                                                  \
    11 /**< originally 8, now 11 to support Spouge's FSC. see notes below */

#define NCBIMATH_LN2 0.69314718055994530941723212145818

/* amino acid background frequencies from Robinson and Robinson */
static BLAST_LetterProb Robinson_prob_ncbi[] = {
    { 'A', 78.05 },
    { 'C', 19.25 },
    { 'D', 53.64 },
    { 'E', 62.95 },
    { 'F', 38.56 },
    { 'G', 73.77 },
    { 'H', 21.99 },
    { 'I', 51.42 },
    { 'K', 57.44 },
    { 'L', 90.19 },
    { 'M', 22.43 },
    { 'N', 44.87 },
    { 'P', 52.03 },
    { 'Q', 42.64 },
    { 'R', 51.29 },
    { 'S', 71.20 },
    { 'T', 58.41 },
    { 'V', 64.41 },
    { 'W', 13.30 },
    { 'Y', 32.16 }
}; /**< amino acid background frequencies from Robinson and Robinson */

const Uint1 AMINOACID_TO_NCBISTDAA[128] = {
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,  0,
    0,  0,  0,  0,  25, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  1,  2,  3,  4,  5,  6,  7,  8, 9, 27, 10,
    11, 12, 13, 26, 14, 15, 16, 17, 18, 24, 19, 20, 21, 22, 23, 0, 0, 0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

static BLAST_LetterProb nt_prob[] = {
    { 'A', 25.00 }, { 'C', 25.00 }, { 'G', 25.00 }, { 'T', 25.00 }
}; /**< nucleotide probabilities (25% each letter) */

#define STD_AMINO_ACID_FREQS                                                   \
    Robinson_prob_ncbi /**< points to the standard amino acid frequencies to     \
                         use. */

#define DIM(A) (sizeof(A) / sizeof((A)[0]))

typedef struct Blast_ResFreq {
    Uint1 alphabet_code; /**< indicates alphabet. */
    double *prob;        /**< letter probs, (possible) non-zero offset. */
    double *prob0;       /**< probs, zero offset. */
} Blast_ResFreq;

Blast_ResFreq *Blast_ResFreqFree(Blast_ResFreq *rfp) {
    if (rfp == NULL)
        return NULL;

    if (rfp->prob0 != NULL)
        free(rfp->prob0);

    free(rfp);

    return rfp;
}

#define ABS(a) ((a) >= 0 ? (a) : -(a))

Int4 BLAST_Gcd(Int4 a, Int4 b) {
    Int4 c;

    b = ABS(b);
    if (b > a)
        c = a, a = b, b = c;

    while (b != 0) {
        c = a % b;
        a = b;
        b = c;
    }
    return a;
}

/*
   Allocates the Blast_ResFreq* and then fills in the frequencies
   in the probabilities.
   */

Blast_ResFreq *Blast_ResFreqNew(const BlastScoreBlk *sbp) {
    Blast_ResFreq *rfp;

    if (sbp == NULL) {
        return NULL;
    }

    rfp = (Blast_ResFreq *)calloc(1, sizeof(Blast_ResFreq));
    if (rfp == NULL)
        return NULL;

    rfp->alphabet_code = sbp->alphabet_code;

    rfp->prob0 = (double *)calloc(sbp->alphabet_size, sizeof(double));
    if (rfp->prob0 == NULL) {
        rfp = Blast_ResFreqFree(rfp);
        return rfp;
    }
    rfp->prob = rfp->prob0 - sbp->alphabet_start;

    return rfp;
}

Int2 Blast_GetStdAlphabet(Uint1 alphabet_code, Uint1 *residues,
        Uint4 residues_size) {
    Int2 index;

    if (residues_size < DIM(STD_AMINO_ACID_FREQS))
        return -2;

    for (index = 0; index < (int)DIM(STD_AMINO_ACID_FREQS); index++) {
        if (alphabet_code == BLASTAA_SEQ_CODE) {
            residues[index] =
                // AMINOACID_TO_NCBISTDAA[toupper((unsigned char)
                // STD_AMINO_ACID_FREQS[index].ch)];
                encoding_codesArray
                [toupper((unsigned char)STD_AMINO_ACID_FREQS[index].ch)];
        } else {
            residues[index] = STD_AMINO_ACID_FREQS[index].ch;
        }
    }

    return index;
}

static Int2 Blast_ResFreqNormalize(const BlastScoreBlk *sbp, Blast_ResFreq *rfp,
        double norm) {
    Int2 alphabet_stop, index;
    double sum = 0., p;

    if (norm == 0.)
        return 1;

    alphabet_stop = sbp->alphabet_start + sbp->alphabet_size;
    for (index = sbp->alphabet_start; index < alphabet_stop; index++) {
        p = rfp->prob[index];
        if (p < 0.)
            return 1;
        sum += p;
    }
    if (sum <= 0.)
        return 0;

    for (index = sbp->alphabet_start; index < alphabet_stop; index++) {
        rfp->prob[index] /= sum;
        rfp->prob[index] *= norm;
    }
    return 0;
}

Int2 Blast_ResFreqStdComp(const BlastScoreBlk *sbp, Blast_ResFreq *rfp) {
    Uint4 index;

    if (sbp->protein_alphabet == TRUE) {
        Int2 retval;
        Uint1 *residues = (Uint1 *)calloc(DIM(STD_AMINO_ACID_FREQS), sizeof(Uint1));
        retval = Blast_GetStdAlphabet(sbp->alphabet_code, residues,
                DIM(STD_AMINO_ACID_FREQS));
        if (retval < 1)
            return retval;

        for (index = 0; index < DIM(STD_AMINO_ACID_FREQS); index++) {
            rfp->prob[residues[index]] = STD_AMINO_ACID_FREQS[index].p;
            // fprintf(stderr, "rfp->prob : %f residues: %d \n",
            // rfp->prob[residues[index]], residues[index]);
        }
        free(residues);
    } else { /* beginning of blastna and ncbi2na are the same. */
        /* Only blastna used  for nucleotides. */
        for (index = 0; index < DIM(nt_prob); index++) {
            rfp->prob[index] = nt_prob[index].p;
        }
    }

    Blast_ResFreqNormalize(sbp, rfp, 1.0);

    return 0;
}

static Int2 BlastScoreChk(Int4 lo, Int4 hi) {
    if (lo >= 0 || hi <= 0 || lo < BLAST_SCORE_MIN || hi > BLAST_SCORE_MAX)
        return 1;

    if (hi - lo > BLAST_SCORE_RANGE_MAX)
        return 1;

    return 0;
}

Blast_ScoreFreq *Blast_ScoreFreqFree(Blast_ScoreFreq *sfp) {
    if (sfp == NULL)
        return NULL;

    if (sfp->sprob0 != NULL)
        free(sfp->sprob0);
    free(sfp);
    return sfp;
}

Blast_ScoreFreq *Blast_ScoreFreqNew(Int4 score_min, Int4 score_max) {
    Blast_ScoreFreq *sfp;
    Int4 range;

    if (BlastScoreChk(score_min, score_max) != 0)
        return NULL;

    sfp = (Blast_ScoreFreq *)calloc(1, sizeof(Blast_ScoreFreq));
    if (sfp == NULL)
        return NULL;

    range = score_max - score_min + 1;
    sfp->sprob = (double *)calloc(range, sizeof(double));
    if (sfp->sprob == NULL) {
        Blast_ScoreFreqFree(sfp);
        return NULL;
    }

    sfp->sprob0 = sfp->sprob;
    sfp->sprob -= score_min; /* center around 0 */
    sfp->score_min = score_min;
    sfp->score_max = score_max;
    sfp->obs_min = sfp->obs_max = 0;
    sfp->score_avg = 0.0;
    return sfp;
}

static Int2 BlastScoreFreqCalc(const BlastScoreBlk *sbp, Blast_ScoreFreq *sfp,
        Blast_ResFreq *rfp1, Blast_ResFreq *rfp2) {
    Int2 **matrix;
    Int4 score, obs_min, obs_max;
    double score_sum, score_avg;
    Int2 alphabet_start, alphabet_end, index1, index2;

    if (sbp == NULL || sfp == NULL)
        return 1;

    if (sbp->loscore < sfp->score_min || sbp->hiscore > sfp->score_max)
        return 1;

    for (score = sfp->score_min; score <= sfp->score_max; score++)
        sfp->sprob[score] = 0.0;

    matrix = sbp->matrix;

    alphabet_start = sbp->alphabet_start;
    alphabet_end = alphabet_start + sbp->alphabet_size;
    for (index1 = alphabet_start; index1 < alphabet_end; index1++) {
        for (index2 = alphabet_start; index2 < alphabet_end; index2++) {
            score = matrix[index1][index2];
            if (score >= sbp->loscore) {
                sfp->sprob[score] += rfp1->prob[index1] * rfp2->prob[index2];
            }
        }
    }

    score_sum = 0.;
    obs_min = obs_max = BLAST_SCORE_MIN;
    for (score = sfp->score_min; score <= sfp->score_max; score++) {
        if (sfp->sprob[score] > 0.) {
            score_sum += sfp->sprob[score];
            obs_max = score;
            if (obs_min == BLAST_SCORE_MIN)
                obs_min = score;
        }
    }
    sfp->obs_min = obs_min;
    sfp->obs_max = obs_max;

    score_avg = 0.0;
    if (score_sum > 0.0001 || score_sum < -0.0001) {
        for (score = obs_min; score <= obs_max; score++) {
            sfp->sprob[score] /= score_sum;
            score_avg += score * sfp->sprob[score];
        }
    }
    sfp->score_avg = score_avg;

    return 0;
}

Blast_KarlinBlk *Blast_KarlinBlkNew(void) {
    Blast_KarlinBlk *kbp;

    kbp = (Blast_KarlinBlk *)calloc(1, sizeof(Blast_KarlinBlk));

    return kbp;
}

static double NlmKarlinLambdaNR(double *probs, Int4 d, Int4 low, Int4 high,
        double lambda0, double tolx, Int4 itmax,
        Int4 maxNewton, Int4 *itn) {
    Int4 k;
    double x0, x, a = 0, b = 1;
    double f = 4;      /* Larger than any possible value of the poly in [0,1] */
    Int4 isNewton = 0; /* we haven't yet taken a Newton step. */

    // assert( d > 0 );

    x0 = exp(-lambda0);
    x = (0 < x0 && x0 < 1) ? x0 : .5;

    for (k = 0; k < itmax; k++) { /* all iteration indices k */
        Int4 i;
        double g, fold = f;
        Int4 wasNewton = isNewton; /* If true, then the previous step was a */
        /* Newton step */
        isNewton = 0;              /* Assume that this step is not */

        /* Horner's rule for evaluating a polynomial and its derivative */
        g = 0;
        f = probs[low];
        for (i = low + d; i < 0; i += d) {
            g = x * g + f;
            f = f * x + probs[i];
        }
        g = x * g + f;
        f = f * x + probs[0] - 1;
        for (i = d; i <= high; i += d) {
            g = x * g + f;
            f = f * x + probs[i];
        }
        /* End Horner's rule */

        if (f > 0) {
            a = x; /* move the left endpoint */
        } else if (f < 0) {
            b = x; /* move the right endpoint */
        } else { /* f == 0 */
            break; /* x is an exact solution */
        }
        if (b - a < 2 * a * (1 - b) * tolx) {
            /* The midpoint of the interval converged */
            x = (a + b) / 2;
            break;
        }

        if (k >= maxNewton ||
                /* If convergence of Newton's method appears to be failing; or */
                (wasNewton && fabs(f) > .9 * fabs(fold)) ||
                /* if the previous iteration was a Newton step but didn't decrease
                 * f sufficiently; or */
                g >= 0
                /* if a Newton step will move us away from the desired solution */
           ) { /* then */
            /* bisect */
            x = (a + b) / 2;
        } else {
            /* try a Newton step */
            double p = -f / g;
            double y = x + p;
            if (y <= a || y >= b) { /* The proposed iterate is not in (a,b) */
                x = (a + b) / 2;
            } else { /* The proposed iterate is in (a,b). Accept it. */
                isNewton = 1;
                x = y;
                if (fabs(p) < tolx * x * (1 - x))
                    break; /* Converged */
            }          /* else the proposed iterate is in (a,b) */
        }            /* else try a Newton step. */
    }              /* end for all iteration indices k */
    *itn = k;
    return -log(x) / d;
}

double Blast_KarlinLambdaNR(Blast_ScoreFreq *sfp, double initialLambdaGuess) {
    Int4 low;  /* Lowest score (must be negative)  */
    Int4 high; /* Highest score (must be positive) */
    Int4 itn;
    Int4 i, d;
    double *sprob;
    double returnValue;

    low = sfp->obs_min;
    high = sfp->obs_max;
    if (sfp->score_avg >= 0.) { /* Expected score must be negative */
        return -1.0;
    }
    if (BlastScoreChk(low, high) != 0)
        return -1.;

    sprob = sfp->sprob;
    /* Find greatest common divisor of all scores */
    for (i = 1, d = -low; i <= high - low && d > 1; ++i) {
        if (sprob[i + low] != 0.0) {
            d = BLAST_Gcd(d, i);
        }
    }
    returnValue = NlmKarlinLambdaNR(sprob, d, low, high, initialLambdaGuess,
            BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT, 20,
            20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT, &itn);

    return returnValue;
}

double BLAST_Powi(double x, Int4 n) {
    double y;

    if (n == 0)
        return 1.;

    if (x == 0.) {
        if (n < 0) {
            return HUGE_VAL;
        }
        return 0.;
    }

    if (n < 0) {
        x = 1. / x;
        n = -n;
    }

    y = 1.;
    while (n > 0) {
        if (n & 1)
            y *= x;
        n /= 2;
        x *= x;
    }
    return y;
}

static double BlastKarlinLtoH(Blast_ScoreFreq *sfp, double lambda) {
    Int4 score;
    double H, etonlam, sum, scale;

    double *probs = sfp->sprob;
    Int4 low = sfp->obs_min, high = sfp->obs_max;

    if (lambda < 0.) {
        return -1.;
    }
    if (BlastScoreChk(low, high) != 0)
        return -1.;

    etonlam = exp(-lambda);
    sum = low * probs[low];
    for (score = low + 1; score <= high; score++) {
        sum = score * probs[score] + etonlam * sum;
    }

    scale = BLAST_Powi(etonlam, high);
    if (scale > 0.0) {
        H = lambda * sum / scale;
    } else { /* Underflow of exp( -lambda * high ) */
        H = lambda * exp(lambda * high + log(sum));
    }
    return H;
}

double BLAST_Expm1(double x) {
    double absx = ABS(x);

    if (absx > .33)
        return exp(x) - 1.;

    if (absx < 1.e-16)
        return x;

    return x *
        (1. +
         x * (1. / 2. +
             x * (1. / 6. +
                 x * (1. / 24. +
                     x * (1. / 120. +
                         x * (1. / 720. +
                             x * (1. / 5040. +
                                 x * (1. / 40320. +
                                     x * (1. / 362880. +
                                         x * (1. / 3628800. +
                                             x * (1. / 39916800. +
                                                 x * (1. /
                                                     479001600. +
                                                     x / 6227020800.))))))))))));
}

static double BlastKarlinLHtoK(Blast_ScoreFreq *sfp, double lambda, double H) {
    /*The next array stores the probabilities of getting each possible
      score in an alignment of fixed length; the array is shifted
      during part of the computation, so that
      entry 0 is for score 0.  */
    double *alignmentScoreProbabilities = NULL;
    Int4 low;        /* Lowest score (must be negative) */
    Int4 high;       /* Highest score (must be positive) */
    Int4 range;      /* range of scores, computed as high - low*/
    double K;        /* local copy of K  to return*/
    int i;           /*loop index*/
    int iterCounter; /*counter on iterations*/
    Int4 divisor;    /*candidate divisor of all scores with
                       non-zero probabilities*/
    /*highest and lowest possible alignment scores for current length*/
    Int4 lowAlignmentScore, highAlignmentScore;
    Int4 first, last; /*loop indices for dynamic program*/
    register double innerSum;
    double oldsum, oldsum2; /* values of innerSum on previous
                               iterations*/
    double outerSum; /* holds sum over j of (innerSum
                        for iteration j/j)*/

    double score_avg; /*average score*/
    /*first term to use in the closed form for the case where
      high == 1 or low == -1, but not both*/
    double firstTermClosedForm; /*usually store H/lambda*/
    int iterlimit;              /*upper limit on iterations*/
    double sumlimit;            /*lower limit on contributions
                                  to sum over scores*/

    /*array of score probabilities reindexed so that low is at index 0*/
    double *probArrayStartLow;

    /*pointers used in dynamic program*/
    double *ptrP, *ptr1, *ptr2, *ptr1e;
    double expMinusLambda; /*e^^(-Lambda) */

    if (lambda <= 0. || H <= 0.) {
        /* Theory dictates that H and lambda must be positive, so
         * return -1 to indicate an error */
        return -1.;
    }

    /*Karlin-Altschul theory works only if the expected score
      is negative*/
    if (sfp->score_avg >= 0.0) {
        return -1.;
    }

    low = sfp->obs_min;
    high = sfp->obs_max;
    range = high - low;

    probArrayStartLow = &sfp->sprob[low];
    /* Look for the greatest common divisor ("delta" in Appendix of PNAS 87 of
       Karlin&Altschul (1990) */
    for (i = 1, divisor = -low; i <= range && divisor > 1; ++i) {
        if (probArrayStartLow[i] != 0.0)
            divisor = BLAST_Gcd(divisor, i);
    }

    high /= divisor;
    low /= divisor;
    lambda *= divisor;

    range = high - low;

    firstTermClosedForm = H / lambda;
    expMinusLambda = exp((double)-lambda);

    if (low == -1 && high == 1) {
        K = (sfp->sprob[low * divisor] - sfp->sprob[high * divisor]) *
            (sfp->sprob[low * divisor] - sfp->sprob[high * divisor]) /
            sfp->sprob[low * divisor];
        return (K);
    }

    if (low == -1 || high == 1) {
        if (high != 1) {
            score_avg = sfp->score_avg / divisor;
            firstTermClosedForm = (score_avg * score_avg) / firstTermClosedForm;
        }
        return firstTermClosedForm * (1.0 - expMinusLambda);
    }

    sumlimit = BLAST_KARLIN_K_SUMLIMIT_DEFAULT;
    iterlimit = BLAST_KARLIN_K_ITER_MAX;

    alignmentScoreProbabilities = (double *)calloc(
            (iterlimit * range + 1), sizeof(*alignmentScoreProbabilities));
    if (alignmentScoreProbabilities == NULL)
        return -1.;

    outerSum = 0.;
    lowAlignmentScore = highAlignmentScore = 0;
    alignmentScoreProbabilities[0] = innerSum = oldsum = oldsum2 = 1.;

    for (iterCounter = 0; ((iterCounter < iterlimit) && (innerSum > sumlimit));
            outerSum += innerSum /= ++iterCounter) {
        first = last = range;
        lowAlignmentScore += low;
        highAlignmentScore += high;
        /*dynamic program to compute P(i,j)*/
        for (ptrP = alignmentScoreProbabilities +
                (highAlignmentScore - lowAlignmentScore);
                ptrP >= alignmentScoreProbabilities; *ptrP-- = innerSum) {
            ptr1 = ptrP - first;
            ptr1e = ptrP - last;
            ptr2 = probArrayStartLow + first;
            for (innerSum = 0.; ptr1 >= ptr1e;) {
                innerSum += *ptr1 * *ptr2;
                ptr1--;
                ptr2++;
            }
            if (first)
                --first;
            if (ptrP - alignmentScoreProbabilities <= range)
                --last;
        }
        /* Horner's rule */
        innerSum = *++ptrP;
        for (i = lowAlignmentScore + 1; i < 0; i++) {
            innerSum = *++ptrP + innerSum * expMinusLambda;
        }
        innerSum *= expMinusLambda;

        for (; i <= highAlignmentScore; ++i)
            innerSum += *++ptrP;
        oldsum2 = oldsum;
        oldsum = innerSum;
    }

#ifdef ADD_GEOMETRIC_TERMS_TO_K
    /*old code assumed that the later terms in sum were
      asymptotically comparable to those of a geometric
      progression, and tried to speed up convergence by
      guessing the estimated ratio between sucessive terms
      and using the explicit terms of a geometric progression
      to speed up convergence. However, the assumption does not
      always hold, and convergenece of the above code is fast
      enough in practice*/
    /* Terms of geometric progression added for correction */
    {
        double ratio; /* fraction used to generate the
                         geometric progression */

        ratio = oldsum / oldsum2;
        if (ratio >= (1.0 - sumlimit * 0.001)) {
            K = -1.;
            if (alignmentScoreProbabilities != NULL)
                sfree(alignmentScoreProbabilities);
            return K;
        }
        sumlimit *= 0.01;
        while (innerSum > sumlimit) {
            oldsum *= ratio;
            outerSum += innerSum = oldsum / ++iterCounter;
        }
    }
#endif

    K = -exp((double)-2.0 * outerSum) /
        (firstTermClosedForm * BLAST_Expm1(-(double)lambda));

    if (alignmentScoreProbabilities != NULL)
        free(alignmentScoreProbabilities);

    return K;
}

Int2 Blast_KarlinBlkUngappedCalc(Blast_KarlinBlk *kbp, Blast_ScoreFreq *sfp) {

    if (kbp == NULL || sfp == NULL)
        return 1;

    /* Calculate the parameter Lambda */

    kbp->Lambda = Blast_KarlinLambdaNR(sfp, BLAST_KARLIN_LAMBDA0_DEFAULT);
    if (kbp->Lambda < 0.)
        goto ErrExit;

    /* Calculate H */

    kbp->H = BlastKarlinLtoH(sfp, kbp->Lambda);
    if (kbp->H < 0.)
        goto ErrExit;

    /* Calculate K and log(K) */

    kbp->K = BlastKarlinLHtoK(sfp, kbp->Lambda, kbp->H);
    if (kbp->K < 0.)
        goto ErrExit;
    kbp->logK = log(kbp->K);

    /* Normal return */
    return 0;

ErrExit:
    kbp->Lambda = kbp->H = kbp->K = -1.;
    kbp->logK = HUGE_VAL;
    return 1;
}

Int2 Blast_ScoreBlkKbpIdealCalc(BlastScoreBlk *sbp) {
    Blast_ResFreq *stdrfp = NULL;
    Blast_ScoreFreq *sfp = NULL;
    Int2 status = 0;

    if (!sbp) {
        return (status = 1);
    }

    stdrfp = Blast_ResFreqNew(sbp);
    Blast_ResFreqStdComp(sbp, stdrfp);
    sfp = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
    BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
    sbp->kbp_ideal = Blast_KarlinBlkNew();
    Blast_KarlinBlkUngappedCalc(sbp->kbp_ideal, sfp);

    stdrfp = Blast_ResFreqFree(stdrfp);
    sfp = Blast_ScoreFreqFree(sfp);

    return status;
}

typedef struct Blast_ResComp {
    Uint1 alphabet_code; /**< indicates alphabet. */
    Int4 *comp;          /**< store composition of a string. */
    Int4 *comp0;         /**< Same array as above, starts at zero. */
} Blast_ResComp;

static Blast_ResComp *BlastResCompDestruct(Blast_ResComp *rcp) {
    if (rcp == NULL)
        return NULL;

    if (rcp->comp0 != NULL)
        free(rcp->comp0);

    free(rcp);
    return NULL;
}

static Blast_ResComp *BlastResCompNew(const BlastScoreBlk *sbp) {
    Blast_ResComp *rcp;

    rcp = (Blast_ResComp *)calloc(1, sizeof(Blast_ResComp));
    if (rcp == NULL)
        return NULL;

    rcp->alphabet_code = sbp->alphabet_code;

    /* comp0 has zero offset, comp starts at 0, only one
       array is allocated.  */
    rcp->comp0 = (Int4 *)calloc(sbp->alphabet_size, sizeof(Int4));
    if (rcp->comp0 == NULL) {
        rcp = BlastResCompDestruct(rcp);
        return rcp;
    }

    rcp->comp = rcp->comp0 - sbp->alphabet_start;

    return rcp;
}

#define BLASTNA_SEQ_CODE 99

const Uint1 IUPACNA_TO_BLASTNA[128] = {
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 0,  10, 1,  11, 15, 15, 2,  12, 15, 15, 7,
    15, 6,  14, 15, 15, 15, 4,  9,  3,  15, 13, 8,  15, 5,  15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15
};

const Uint1 IUPACNA_TO_NCBI4NA[128] = {
    0,  0,  0,  0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,  0,  0,  0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,  0,  0,  0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 1,
    14, 2,  13, 0, 0, 4, 11, 0, 0, 12, 0, 3, 15, 0, 0, 0, 5, 6, 8, 0, 7, 9,
    0,  10, 0,  0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,  0,  0,  0, 0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0
};

#define NCBI4NA_SEQ_CODE 4

Int2 BLAST_ScoreSetAmbigRes(BlastScoreBlk *sbp, char ambiguous_res) {
    Int2 index;
    Uint1 *ambig_buffer;

    if (sbp == NULL)
        return 1;

    if (sbp->ambig_occupy >= sbp->ambig_size) {
        sbp->ambig_size += 5;
        ambig_buffer = (Uint1 *)calloc(sbp->ambig_size, sizeof(Uint1));
        for (index = 0; index < sbp->ambig_occupy; index++) {
            ambig_buffer[index] = sbp->ambiguous_res[index];
        }
        free(sbp->ambiguous_res);
        sbp->ambiguous_res = ambig_buffer;
    }

    if (sbp->alphabet_code == BLASTAA_SEQ_CODE) {
        sbp->ambiguous_res[sbp->ambig_occupy] =
            AMINOACID_TO_NCBISTDAA[toupper((unsigned char)ambiguous_res)];
    } else {
        if (sbp->alphabet_code == BLASTNA_SEQ_CODE)
            sbp->ambiguous_res[sbp->ambig_occupy] =
                IUPACNA_TO_BLASTNA[toupper((unsigned char)ambiguous_res)];
        else if (sbp->alphabet_code == NCBI4NA_SEQ_CODE)
            sbp->ambiguous_res[sbp->ambig_occupy] =
                IUPACNA_TO_NCBI4NA[toupper((unsigned char)ambiguous_res)];
    }
    (sbp->ambig_occupy)++;

    return 0;
}

static Int2 BlastResCompStr(const BlastScoreBlk *sbp, Blast_ResComp *rcp,
        char *str, Int4 length) {
    char *lp, *lpmax;
    Int2 index;
    Uint1 mask;

    if (sbp == NULL || rcp == NULL || str == NULL)
        return 1;

    if (rcp->alphabet_code != sbp->alphabet_code)
        return 1;

    /* For megablast, check only the first 4 bits of the sequence values */
    mask = (sbp->protein_alphabet ? 0xff : 0x0f);

    /* comp0 starts at zero and extends for "num", comp is the same array, but
       "start_at" offset. */
    for (index = 0; index < (sbp->alphabet_size); index++)
        rcp->comp0[index] = 0;

    for (lp = str, lpmax = lp + length; lp < lpmax; lp++) {
        // fprintf(stderr, "%c", encoding_getLetter(*lp));
        ++rcp->comp[(int)(*lp & mask)];
    }

    // for(index = 0; index < (sbp->alphabet_size); index++)
    //{
    // fprintf(stderr, "%d : %d\n", index, rcp->comp[index]);
    //}

    // fprintf(stderr, "\n");

    /* Don't count ambig. residues. */
    for (index = 0; index < sbp->ambig_occupy; index++) {
        rcp->comp[sbp->ambiguous_res[index]] = 0;
    }

    return 0;
}

static Int2 Blast_ResFreqClr(const BlastScoreBlk *sbp, Blast_ResFreq *rfp) {
    Int2 alphabet_max, index;

    alphabet_max = sbp->alphabet_start + sbp->alphabet_size;
    for (index = sbp->alphabet_start; index < alphabet_max; index++)
        rfp->prob[index] = 0.0;

    return 0;
}

static Int2 Blast_ResFreqResComp(const BlastScoreBlk *sbp, Blast_ResFreq *rfp,
        const Blast_ResComp *rcp) {
    Int2 alphabet_max, index;
    double sum = 0.;

    if (rfp == NULL || rcp == NULL)
        return 1;

    if (rfp->alphabet_code != rcp->alphabet_code)
        return 1;

    alphabet_max = sbp->alphabet_start + sbp->alphabet_size;
    for (index = sbp->alphabet_start; index < alphabet_max; index++) {
        sum += rcp->comp[index];
    }

    if (sum == 0.) {
        Blast_ResFreqClr(sbp, rfp);
        return 0;
    }

    for (index = sbp->alphabet_start; index < alphabet_max; index++)
        rfp->prob[index] = rcp->comp[index] / sum;

    return 0;
}

static Int2 Blast_ResFreqString(const BlastScoreBlk *sbp, Blast_ResFreq *rfp,
        char *string, Int4 length) {
    Blast_ResComp *rcp;

    rcp = BlastResCompNew(sbp);

    BlastResCompStr(sbp, rcp, string, length);

    Blast_ResFreqResComp(sbp, rfp, rcp);

    rcp = BlastResCompDestruct(rcp);

    return 0;
}

Blast_KarlinBlk *Blast_KarlinBlkFree(Blast_KarlinBlk *kbp) {
    free(kbp);

    return kbp;
}

Int2 Blast_ScoreBlkKbpUngappedCalc(BlastScoreBlk *sbp,
        // int num_query,
        struct PSSMatrix PSSMatrix,
        Int4 *blast_ungappedNominalTrigger, Blast_KarlinBlk *kbp) {
    
    Int2 status = 0;
    Int4 context; /* loop variable. */
    Blast_ResFreq *rfp, *stdrfp;

    /* Ideal Karlin block is filled unconditionally. */
    status = Blast_ScoreBlkKbpIdealCalc(sbp);
    if (status)
        return status;

    stdrfp = Blast_ResFreqNew(sbp);
    Blast_ResFreqStdComp(sbp, stdrfp);
    rfp = Blast_ResFreqNew(sbp);

    Int4 context_offset;
    Int4 query_length;
    Uint1 *buffer; /* holds sequence */
    Int2 loop_status; /* status flag for functions in this loop. */

    query_length = PSSMatrix.length;
    buffer = PSSMatrix.queryCodes;

    Blast_ResFreqString(sbp, rfp, (char *)buffer, query_length);
    Blast_ScoreFreq *sfp = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
    BlastScoreFreqCalc(sbp, sfp, rfp, stdrfp);
    loop_status = Blast_KarlinBlkUngappedCalc(kbp, sfp);
    *blast_ungappedNominalTrigger =
        (Int4)((22 * NCBIMATH_LN2 + kbp->logK) / kbp->Lambda) * sbp->scale_factor;

    Blast_ScoreFreqFree(sfp);
    Blast_ResFreqFree(stdrfp);
    Blast_ResFreqFree(rfp);
    // printf("blast_ungappedNominalTrigger: %d\n",
    // *blast_ungappedNominalTrigger);


    return status;
}

static const double one = 1.0, halF[2] = { 0.5, -0.5, }, huge = 1.0e+300,
             twom1000 =
             9.33263618503218878990e-302, /* 2**-1000=0x01700000,0*/
             o_threshold = 7.09782712893383973096e+02,    /* 0x40862E42, 0xFEFA39EF */
             u_threshold = -7.45133219101941108420e+02,   /* 0xc0874910, 0xD52D3051 */
             ln2HI[2] = { 6.93147180369123816490e-01,     /* 0x3fe62e42, 0xfee00000 */
                 -6.93147180369123816490e-01, }, /* 0xbfe62e42, 0xfee00000 */
                 ln2LO[2] = { 1.90821492927058770002e-10,     /* 0x3dea39ef, 0x35793c76 */
                     -1.90821492927058770002e-10, }, /* 0xbdea39ef, 0x35793c76 */
                     invln2 = 1.44269504088896338700e+00,         /* 0x3ff71547, 0x652b82fe */
                     P1 = 1.66666666666666019037e-01,             /* 0x3FC55555, 0x5555553E */
                     P2 = -2.77777777770155933842e-03,            /* 0xBF66C16C, 0x16BEBD93 */
                     P3 = 6.61375632143793436117e-05,             /* 0x3F11566A, 0xAF25DE2C */
                     P4 = -1.65339022054652515390e-06,            /* 0xBEBBBD41, 0xC5D26BF1 */
                     P5 = 4.13813679705723846039e-08;             /* 0x3E663769, 0x72BEA4D0 */

#ifndef HAVE_ERF
static const double tiny = 1e-300,
             half =
             5.00000000000000000000e-01, /* 0x3FE00000, 0x00000000 */
             /* one =  1.00000000000000000000e+00, */        /* 0x3FF00000, 0x00000000 */
             two = 2.00000000000000000000e+00,               /* 0x40000000, 0x00000000 */
             /* c = (float)0.84506291151 */
             erx = 8.45062911510467529297e-01, /* 0x3FEB0AC1, 0x60000000 */
             /*
              * Coefficients for approximation to  erf on [0,0.84375]
              */
             efx = 1.28379167095512586316e-01,  /* 0x3FC06EBA, 0x8214DB69 */
             efx8 = 1.02703333676410069053e+00, /* 0x3FF06EBA, 0x8214DB69 */
             pp0 = 1.28379167095512558561e-01,  /* 0x3FC06EBA, 0x8214DB68 */
             pp1 = -3.25042107247001499370e-01, /* 0xBFD4CD7D, 0x691CB913 */
             pp2 = -2.84817495755985104766e-02, /* 0xBF9D2A51, 0xDBD7194F */
             pp3 = -5.77027029648944159157e-03, /* 0xBF77A291, 0x236668E4 */
             pp4 = -2.37630166566501626084e-05, /* 0xBEF8EAD6, 0x120016AC */
             qq1 = 3.97917223959155352819e-01,  /* 0x3FD97779, 0xCDDADC09 */
             qq2 = 6.50222499887672944485e-02,  /* 0x3FB0A54C, 0x5536CEBA */
             qq3 = 5.08130628187576562776e-03,  /* 0x3F74D022, 0xC4D36B0F */
             qq4 = 1.32494738004321644526e-04,  /* 0x3F215DC9, 0x221C1A10 */
             qq5 = -3.96022827877536812320e-06, /* 0xBED09C43, 0x42A26120 */
             /*
              * Coefficients for approximation to  erf  in [0.84375,1.25]
              */
             pa0 = -2.36211856075265944077e-03, /* 0xBF6359B8, 0xBEF77538 */
             pa1 = 4.14856118683748331666e-01,  /* 0x3FDA8D00, 0xAD92B34D */
             pa2 = -3.72207876035701323847e-01, /* 0xBFD7D240, 0xFBB8C3F1 */
             pa3 = 3.18346619901161753674e-01,  /* 0x3FD45FCA, 0x805120E4 */
             pa4 = -1.10894694282396677476e-01, /* 0xBFBC6398, 0x3D3E28EC */
             pa5 = 3.54783043256182359371e-02,  /* 0x3FA22A36, 0x599795EB */
             pa6 = -2.16637559486879084300e-03, /* 0xBF61BF38, 0x0A96073F */
             qa1 = 1.06420880400844228286e-01,  /* 0x3FBB3E66, 0x18EEE323 */
             qa2 = 5.40397917702171048937e-01,  /* 0x3FE14AF0, 0x92EB6F33 */
             qa3 = 7.18286544141962662868e-02,  /* 0x3FB2635C, 0xD99FE9A7 */
             qa4 = 1.26171219808761642112e-01,  /* 0x3FC02660, 0xE763351F */
             qa5 = 1.36370839120290507362e-02,  /* 0x3F8BEDC2, 0x6B51DD1C */
             qa6 = 1.19844998467991074170e-02,  /* 0x3F888B54, 0x5735151D */
             /*
              * Coefficients for approximation to  erfc in [1.25,1/0.35]
              */
             ra0 = -9.86494403484714822705e-03, /* 0xBF843412, 0x600D6435 */
             ra1 = -6.93858572707181764372e-01, /* 0xBFE63416, 0xE4BA7360 */
             ra2 = -1.05586262253232909814e+01, /* 0xC0251E04, 0x41B0E726 */
             ra3 = -6.23753324503260060396e+01, /* 0xC04F300A, 0xE4CBA38D */
             ra4 = -1.62396669462573470355e+02, /* 0xC0644CB1, 0x84282266 */
             ra5 = -1.84605092906711035994e+02, /* 0xC067135C, 0xEBCCABB2 */
             ra6 = -8.12874355063065934246e+01, /* 0xC0545265, 0x57E4D2F2 */
             ra7 = -9.81432934416914548592e+00, /* 0xC023A0EF, 0xC69AC25C */
             sa1 = 1.96512716674392571292e+01,  /* 0x4033A6B9, 0xBD707687 */
             sa2 = 1.37657754143519042600e+02,  /* 0x4061350C, 0x526AE721 */
             sa3 = 4.34565877475229228821e+02,  /* 0x407B290D, 0xD58A1A71 */
             sa4 = 6.45387271733267880336e+02,  /* 0x40842B19, 0x21EC2868 */
             sa5 = 4.29008140027567833386e+02,  /* 0x407AD021, 0x57700314 */
             sa6 = 1.08635005541779435134e+02,  /* 0x405B28A3, 0xEE48AE2C */
             sa7 = 6.57024977031928170135e+00,  /* 0x401A47EF, 0x8E484A93 */
             sa8 = -6.04244152148580987438e-02, /* 0xBFAEEFF2, 0xEE749A62 */
             /*
              * Coefficients for approximation to  erfc in [1/.35,28]
              */
             rb0 = -9.86494292470009928597e-03, /* 0xBF843412, 0x39E86F4A */
             rb1 = -7.99283237680523006574e-01, /* 0xBFE993BA, 0x70C285DE */
             rb2 = -1.77579549177547519889e+01, /* 0xC031C209, 0x555F995A */
             rb3 = -1.60636384855821916062e+02, /* 0xC064145D, 0x43C5ED98 */
             rb4 = -6.37566443368389627722e+02, /* 0xC083EC88, 0x1375F228 */
             rb5 = -1.02509513161107724954e+03, /* 0xC0900461, 0x6A2E5992 */
             rb6 = -4.83519191608651397019e+02, /* 0xC07E384E, 0x9BDC383F */
             sb1 = 3.03380607434824582924e+01,  /* 0x403E568B, 0x261D5190 */
             sb2 = 3.25792512996573918826e+02,  /* 0x40745CAE, 0x221B9F0A */
             sb3 = 1.53672958608443695994e+03,  /* 0x409802EB, 0x189D5118 */
             sb4 = 3.19985821950859553908e+03,  /* 0x40A8FFB7, 0x688C246A */
             sb5 = 2.55305040643316442583e+03,  /* 0x40A3F219, 0xCEDF3BE6 */
             sb6 = 4.74528541206955367215e+02,  /* 0x407DA874, 0xE79FE763 */
             sb7 = -2.24409524465858183362e+01; /* 0xC03670E2, 0x42712D62 */
#endif

#ifdef WORDS_BIGENDIAN
#define __HI(x) *(int *)&x
#define __LO(x) *(1 + (int *)&x)
#define __HIp(x) *(int *)x
#define __LOp(x) *(1 + (int *)x)
#else
#define __HI(x) *(1 + (int *)&x)
#define __LO(x) *(int *)&x
#define __HIp(x) *(1 + (int *)x)
#define __LOp(x) *(int *)x
#endif

static double s_IEEE754_Exp(double x) /* default IEEE double exp */
{
    double y, hi, lo, c, t;
    int k, xsb;
    unsigned hx;

    hx = __HI(x);         /* high word of x */
    xsb = (hx >> 31) & 1; /* sign bit of x */
    hx &= 0x7fffffff;     /* high word of |x| */

    /* filter out non-finite argument */
    if (hx >= 0x40862E42) { /* if |x|>=709.78... */
        if (hx >= 0x7ff00000) {
            if (((hx & 0xfffff) | __LO(x)) != 0)
                return x + x; /* NaN */
            else
                return (xsb == 0) ? x : 0.0; /* exp(+-inf)={inf,0} */
        }
        if (x > o_threshold)
            return huge * huge; /* NCBI_FAKE_WARNING [deliberate overflow] */
        if (x < u_threshold)
            return twom1000 * twom1000; /* underflow */
    }

    /* argument reduction */
    if (hx > 0x3fd62e42) {   /* if  |x| > 0.5 ln2 */
        if (hx < 0x3FF0A2B2) { /* and |x| < 1.5 ln2 */
            hi = x - ln2HI[xsb];
            lo = ln2LO[xsb];
            k = 1 - xsb - xsb;
        } else {
            k = (int)(invln2 * x + halF[xsb]);
            t = k;
            hi = x - t * ln2HI[0]; /* t*ln2HI is exact here */
            lo = t * ln2LO[0];
        }
        x = hi - lo;
    } else if (hx < 0x3e300000) { /* when |x|<2**-28 */
        if (huge + x > one)
            return one + x; /* trigger inexact */
    } else
        k = 0;

    /* x is now in primary range */
    t = x * x;
    c = x - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
    if (k == 0)
        return one - ((x * c) / (c - 2.0) - x);
    else
        y = one - ((lo - (x * c) / (2.0 - c)) - hi);
    if (k >= -1021) {
        __HI(y) += (k << 20); /* add k to y's exponent */
        return y;
    } else {
        __HI(y) += ((k + 1000) << 20); /* add k to y's exponent */
        return y * twom1000;
    }
}

double NCBI_Erf(double x) {
#ifdef HAVE_ERF
    return erf(x);
#else
    int hx, ix, i;
    double R, S, P, Q, s, y, z, r;
    hx = __HI(x);
    ix = hx & 0x7fffffff;
    if (ix >= 0x7ff00000) { /* erf(nan)=nan */
        i = ((unsigned)hx >> 31) << 1;
        return (double)(1 - i) + one / x; /* erf(+-inf)=+-1 */
    }

    if (ix < 0x3feb0000) {   /* |x|<0.84375 */
        if (ix < 0x3e300000) { /* |x|<2**-28 */
            if (ix < 0x00800000)
                return 0.125 * (8.0 * x + efx8 * x); /*avoid underflow */
            return x + efx * x;
        }
        z = x * x;
        r = pp0 + z * (pp1 + z * (pp2 + z * (pp3 + z * pp4)));
        s = one + z * (qq1 + z * (qq2 + z * (qq3 + z * (qq4 + z * qq5))));
        y = r / s;
        return x + x * y;
    }
    if (ix < 0x3ff40000) { /* 0.84375 <= |x| < 1.25 */
        s = fabs(x) - one;
        P = pa0 +
            s * (pa1 + s * (pa2 + s * (pa3 + s * (pa4 + s * (pa5 + s * pa6)))));
        Q = one +
            s * (qa1 + s * (qa2 + s * (qa3 + s * (qa4 + s * (qa5 + s * qa6)))));
        if (hx >= 0)
            return erx + P / Q;
        else
            return -erx - P / Q;
    }
    if (ix >= 0x40180000) { /* inf>|x|>=6 */
        if (hx >= 0)
            return one - tiny;
        else
            return tiny - one;
    }
    x = fabs(x);
    s = one / (x * x);
    if (ix < 0x4006DB6E) { /* |x| < 1/0.35 */
        R = ra0 +
            s * (ra1 +
                    s * (ra2 +
                        s * (ra3 + s * (ra4 + s * (ra5 + s * (ra6 + s * ra7))))));
        S = one +
            s * (sa1 +
                    s * (sa2 +
                        s * (sa3 +
                            s * (sa4 +
                                s * (sa5 + s * (sa6 + s * (sa7 + s * sa8)))))));
    } else { /* |x| >= 1/0.35 */
        R = rb0 +
            s * (rb1 + s * (rb2 + s * (rb3 + s * (rb4 + s * (rb5 + s * rb6)))));
        S = one +
            s * (sb1 +
                    s * (sb2 +
                        s * (sb3 + s * (sb4 + s * (sb5 + s * (sb6 + s * sb7))))));
    }
    z = x;
    __LO(z) = 0;
    r = s_IEEE754_Exp(-z * z - 0.5625) * s_IEEE754_Exp((z - x) * (z + x) + R / S);
    if (hx >= 0)
        return one - r / x;
    else
        return r / x - one;
#endif
}

#define BLAST_Erf NCBI_Erf

double BLAST_SpougeStoE(Int4 y_, Blast_KarlinBlk *kbp, Blast_GumbelBlk *gbp,
        Int4 m_, Int4 n_) {
    // the score and lambda may have been rescaled.  We will compute the scaling
    // factor
    // and use it to scale a, alpha, and Sigma similarly.
    double scale_factor = kbp->Lambda / gbp->Lambda;

    // the pair-wise e-value must be scaled back to db-wise e-value
    double db_scale_factor =
        (gbp->db_length) ? (double)gbp->db_length / (double)n_ : 1.0;

    double lambda_ = kbp->Lambda;
    double k_ = kbp->K;
    double ai_hat_ = gbp->a * scale_factor;
    double bi_hat_ = gbp->b;
    double alphai_hat_ = gbp->Alpha * scale_factor;
    double betai_hat_ = gbp->Beta;
    double sigma_hat_ = gbp->Sigma * scale_factor;
    double tau_hat_ = gbp->Tau;

    /* here we consider symmetric matrix only */
    double aj_hat_ = ai_hat_;
    double bj_hat_ = bi_hat_;
    double alphaj_hat_ = alphai_hat_;
    double betaj_hat_ = betai_hat_;

    /* this is 1/sqrt(2.0*PI) */
    static double const_val = 0.39894228040143267793994605993438;

    double m_li_y, vi_y, sqrt_vi_y, m_F, P_m_F;
    double n_lj_y, vj_y, sqrt_vj_y, n_F, P_n_F;
    double c_y, p1, p2, area;

    m_li_y = m_ - (ai_hat_ * y_ + bi_hat_);
    vi_y = MAX(2.0 * alphai_hat_ / lambda_, alphai_hat_ * y_ + betai_hat_);
    sqrt_vi_y = sqrt(vi_y);
    m_F = m_li_y / sqrt_vi_y;
    P_m_F = 0.5 + 0.5 * BLAST_Erf(m_F);
    p1 = m_li_y * P_m_F + sqrt_vi_y * const_val * exp(-0.5 * m_F * m_F);

    n_lj_y = n_ - (aj_hat_ * y_ + bj_hat_);
    vj_y = MAX(2.0 * alphaj_hat_ / lambda_, alphaj_hat_ * y_ + betaj_hat_);
    sqrt_vj_y = sqrt(vj_y);
    n_F = n_lj_y / sqrt_vj_y;
    P_n_F = 0.5 + 0.5 * BLAST_Erf(n_F);
    p2 = n_lj_y * P_n_F + sqrt_vj_y * const_val * exp(-0.5 * n_F * n_F);

    c_y = MAX(2.0 * sigma_hat_ / lambda_, sigma_hat_ * y_ + tau_hat_);
    area = p1 * p2 + c_y * P_m_F * P_n_F;

    return area * k_ * exp(-lambda_ * y_) * db_scale_factor;
}

Int4 BLAST_SpougeEtoS(double e0, Blast_KarlinBlk *kbp, Blast_GumbelBlk *gbp,
        Int4 m, Int4 n) {
    Int4 a, b, c;
    double e;
    double db_scale_factor = (gbp->db_length) ? (double)gbp->db_length : 1.0;

    b = MAX((int)(log(db_scale_factor / e0) / kbp->Lambda), 2);

    e = BLAST_SpougeStoE(b, kbp, gbp, m, n);

    //printf("b: %d e: %f\n", b, e);

    if (e > e0) {
        while (e > e0) {
            a = b;
            b *= 2;
            e = BLAST_SpougeStoE(b, kbp, gbp, m, n);
        }
    } else {
        a = 0;
    }
    while (b - a > 1) {
        c = (a + b) / 2;
        e = BLAST_SpougeStoE(c, kbp, gbp, m, n);
        if (e > e0) {
            a = c;
        } else {
            b = c;
        }
    }
    return a;
}


double
BLAST_GapDecayDivisor(double decayrate, unsigned nsegs )
{
    return (1. - decayrate) * BLAST_Powi(decayrate, nsegs - 1);
}

static double
s_BlastGetBestEvalue(struct ungappedExtension **hsp_array, Int4 hspcnt)
{
    int index = 0;
    double best_evalue = (double) INT4_MAX;

    for (index=0; index<hspcnt; index++)
       best_evalue = MIN(hsp_array[index]->eValue, best_evalue);

    return best_evalue;
}

Int2 Blast_HSPListGetEvalues(
        //EBlastProgramType program_number,
        //const BlastQueryInfo* query_info,
        int query_length,
        Int4 subject_length,
        Blast_KarlinBlk* kbp,
        Blast_GumbelBlk* gbp,
        struct ungappedExtension **ungappedExtensions,
        //struct alignment *alignment,
        Int4 hsp_cnt,
        //BlastHSPList* hsp_list, 
        Boolean gapped_calculation, 
        Boolean RPS_prelim,
        double gap_decay_rate,
        double scaling_factor,
        double *best_eValue)
{
    //BlastHSP* hsp;
    //BlastHSP** hsp_array;
    //Int4 hsp_cnt;
    Int4 index;
    Int4 kbp_context;
    double gap_decay_divisor = 1.;

    //struct ungappedExtension** ungappedExtensions;
    //Boolean isRPS = Blast_ProgramIsRpsBlast(program_number);

    //if (hsp_list == NULL || hsp_list->hspcnt == 0)
    //return 0;

    //kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
    //hsp_cnt = hsp_list->hspcnt;
    //hsp_array = hsp_list->hsp_array;

    if (gap_decay_rate != 0.)
        gap_decay_divisor = BLAST_GapDecayDivisor(gap_decay_rate, 1);

    for (index=0; index<hsp_cnt; index++) {
        //hsp = hsp_array[index];

        //ASSERT(hsp != NULL);
        //ASSERT(scaling_factor != 0.0);
        //ASSERT(sbp->round_down == FALSE || (hsp->score & 1) == 0);

        /* Divide Lambda by the scaling factor, so e-value is 
           calculated correctly from a scaled score. This is needed only
           for RPS BLAST, where scores are scaled, but Lambda is not. */
        //kbp_context = hsp->context;
        //if (RPS_prelim) {
        ///* All kbp in preliminary stage are equivalent.  However, some
        //may be invalid.  Search for the first populated kbp */
        //int i;
        //for (i=0; i<6; ++i) {
        //if (kbp[i]) break;
        //}
        //kbp_context = i;
        //}
        if(ungappedExtensions[index]->status == ungappedExtension_DELETED)
            continue;

        kbp->Lambda /= scaling_factor;

        if (gbp) {
            /* Only try Spouge's method if gumbel parameters are available */
            //if (ungappedExtensions[index]->status == ungappedExtension_SEMIGAPPED || ungappedExtensions[index]->status == ungappedExtension_GAPPED)
            {
                ungappedExtensions[index]->eValue =
                    BLAST_SpougeStoE(ungappedExtensions[index]->nominalScore, kbp, gbp, 
                            query_length, 
                            subject_length);
                //fprintf(stderr, "evalue: %f score: %d subjectlength: %d query_length: %d\n", ungappedExtensions[index]->eValue, ungappedExtensions[index]->nominalScore, subject_length, query_length);
            }
            //else {
            ///* for RPS blast, query and subject is swapped */
            //hsp->evalue =
            //BLAST_SpougeStoE(hsp->score, kbp[kbp_context], sbp->gbp, 
            //subject_length,
            //query_info->contexts[hsp->context].query_length);
            //}
        } 
        //else {
        ///* Get effective search space from the query information block */
        //hsp->evalue =
        //BLAST_KarlinStoE_simple(hsp->score, kbp[kbp_context], 
        //query_info->contexts[hsp->context].eff_searchsp);
        //}

        ungappedExtensions[index]->eValue /= gap_decay_divisor;
        /* Put back the unscaled value of Lambda. */
        kbp->Lambda *= scaling_factor;
    }

    /* Assign the best e-value field. Here the best e-value will always be
       attained for the first HSP in the list. Check that the incoming
       HSP list is properly sorted by score. */
    //ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
    *best_eValue = s_BlastGetBestEvalue(ungappedExtensions, hsp_cnt);
    return 0;
}

Int2 
s_Blast_HSPListReapByPrelimEvalue(struct ungappedExtension **ungappedExtensions, int4 numExts, double prelim_evalue)
{
    //BlastHSP* hsp;
    //BlastHSP** hsp_array;
    //Int4 hsp_cnt = 0;
    Int4 index;
    double cutoff;

    //if (hsp_list == NULL)
    //return 0;

    cutoff = prelim_evalue;

    //hsp_array = hsp_list->hsp_array;
    for (index = 0; index < numExts; index++) {
        //hsp = hsp_array[index];

        //ASSERT(hsp != NULL);

        if(ungappedExtensions[index]->status == ungappedExtension_SEMIGAPPED)
        {
            //fprintf(stderr, "q_start: %d s_start: %d evalue: %f\n", ungappedExtensions[index]->start.queryOffset, ungappedExtensions[index]->start.subjectOffset, ungappedExtensions[index]->eValue);
            if (ungappedExtensions[index]->eValue > cutoff) {
                ungappedExtensions[index]->status = ungappedExtension_DELETED;
            }

        }

    }

    //hsp_list->hspcnt = hsp_cnt;

    return 0;
}


/** Callback for sorting HSPs by ending offset in query. Sorting is by
 * increasing context, then increasing query end offset, then increasing
 * subject end offset, then decreasing score, then decreasing query start
 * offset, then decreasing subject start offset. Null HSPs are moved to the 
 * end of the array.
 * @param v1 pointer to first HSP [in]
 * @param v2 pointer to second HSP [in]
 * @return Result of comparison.
 */
static int
s_QueryEndCompareHSPs(const void* v1, const void* v2)
{
   struct ungappedExtension *h1, *h2;
   struct ungappedExtension **hp1, **hp2;

   hp1 = (struct ungappedExtension**) v1;
   hp2 = (struct ungappedExtension**) v2;
   h1 = *hp1;
   h2 = *hp2;

   if (!h1 && !h2)
      return 0;
   else if (!h1) 
      return 1;
   else if (!h2)
      return -1;

   /* If these are from different contexts, don't compare offsets */
   //if (h1->context < h2->context) 
      //return -1;
   //if (h1->context > h2->context)
      //return 1;

   if (h1->end.queryOffset < h2->end.queryOffset)
      return -1;
   if (h1->end.queryOffset > h2->end.queryOffset)
      return 1;

   if (h1->end.subjectOffset < h2->end.subjectOffset)
      return -1;
   if (h1->end.subjectOffset > h2->end.subjectOffset)
      return 1;

   /* tie breakers: sort by decreasing score, then 
      by increasing size of query range, then by
      increasing size of subject range. The shortest range 
      means the *largest* sequence offset must come 
      first */
   if (h1->nominalScore < h2->nominalScore)
      return 1;
   if (h1->nominalScore > h2->nominalScore)
      return -1;

   if (h1->start.queryOffset < h2->start.queryOffset)
      return 1;
   if (h1->start.queryOffset > h2->start.queryOffset)
      return -1;

   if (h1->start.subjectOffset < h2->start.subjectOffset)
      return 1;
   if (h1->start.subjectOffset > h2->start.subjectOffset)
      return -1;

   return 0;
}


static int
s_QueryOffsetCompareHSPs(const void* v1, const void* v2)
{
   struct ungappedExtension *h1, *h2;
   struct ungappedExtension **hp1, **hp2;

   hp1 = (struct ungappedExtension**) v1;
   hp2 = (struct ungappedExtension**) v2;
   h1 = *hp1;
   h2 = *hp2;

   if (!h1 && !h2)
      return 0;
   else if (!h1) 
      return 1;
   else if (!h2)
      return -1;

   /* If these are from different contexts, don't compare offsets */
   //if (h1->context < h2->context) 
      //return -1;
   //if (h1->context > h2->context)
      //return 1;

   if (h1->start.queryOffset < h2->start.queryOffset)
      return -1;
   if (h1->start.queryOffset > h2->start.queryOffset)
      return 1;

   if (h1->start.subjectOffset < h2->start.subjectOffset)
      return -1;
   if (h1->start.subjectOffset > h2->start.subjectOffset)
      return 1;

   /* tie breakers: sort by decreasing nominalScore, then 
      by increasing size of query range, then by
      increasing subject range. */

   if (h1->nominalScore < h2->nominalScore)
      return 1;
   if (h1->nominalScore > h2->nominalScore)
      return -1;

   if (h1->end.queryOffset < h2->end.queryOffset)
      return 1;
   if (h1->end.queryOffset > h2->end.queryOffset)
      return -1;

   if (h1->end.subjectOffset < h2->end.subjectOffset)
      return 1;
   if (h1->end.subjectOffset > h2->end.subjectOffset)
      return -1;

   return 0;
}

Int4
Blast_HSPListPurgeHSPsWithCommonEndpoints( 
        struct ungappedExtension **hsp_array,
        Int4 hsp_count,
        Boolean purge)

{
    //BlastHSP** hsp_array;  /* hsp_array to purge. */
    struct ungappedExtension* hsp;
    Int4 i, j, k;
    //Int4 hsp_count;
    //purge |= (program != eBlastTypeBlastn);

    /* If HSP list is empty, return immediately. */
    //if (hsp_list == NULL || hsp_list->hspcnt == 0)
    //return 0;

    /* Do nothing for PHI BLAST, because HSPs corresponding to different pattern
       occurrences may have common end points, but should all be kept. */
    //if (Blast_ProgramIsPhiBlast(program))
    //return hsp_list->hspcnt;

    //hsp_array = hsp_list->hsp_array;
    //hsp_count = hsp_list->hspcnt;

    //qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
    qsort(hsp_array, hsp_count, sizeof(struct ungappedExtension *), s_QueryOffsetCompareHSPs);
    i = 0;
    while (i < hsp_count) {
        j = 1;
        while (i+j < hsp_count &&
                hsp_array[i] && hsp_array[i+j] &&
                hsp_array[i]->start.queryOffset == hsp_array[i+j]->start.queryOffset &&
                hsp_array[i]->start.subjectOffset == hsp_array[i+j]->start.subjectOffset) {
            hsp_count--;
            //hsp_array[i+j]->status = ungappedExtension_DELETED;
            hsp = hsp_array[i+j];
            //if (!purge && (hsp->query.end > hsp_array[i]->query.end)) {
            //s_CutOffGapEditScript(hsp, hsp_array[i]->query.end,
            //hsp_array[i]->subject.end, TRUE);
            //} else {
            //hsp = Blast_HSPFree(hsp);
            hsp->status = ungappedExtension_DELETED;
            //hsp = NULL;
            //}
            for (k=i+j; k<hsp_count; k++) {
                hsp_array[k] = hsp_array[k+1];
            }
            hsp_array[hsp_count] = hsp;
        }
        i += j;
    }

    qsort(hsp_array, hsp_count, sizeof(struct ungappedExtension *), s_QueryEndCompareHSPs);
    i = 0;
    while (i < hsp_count) {
        j = 1;
        while (i+j < hsp_count &&
                hsp_array[i] && hsp_array[i+j] &&
                hsp_array[i]->end.queryOffset == hsp_array[i+j]->end.queryOffset &&
                hsp_array[i]->end.subjectOffset == hsp_array[i+j]->end.subjectOffset) {
            hsp_count--;
            //hsp_array[i+j]->status = ungappedExtension_DELETED;
            hsp = hsp_array[i+j];
            //if (!purge && (hsp->query.offset < hsp_array[i]->query.offset)) {
            //s_CutOffGapEditScript(hsp, hsp_array[i]->query.offset,
            //hsp_array[i]->subject.offset, FALSE);
            //} else {
            //hsp = Blast_HSPFree(hsp);
            hsp->status = ungappedExtension_DELETED;
            //hsp = NULL;
            //}
            for (k=i+j; k<hsp_count; k++) {
                hsp_array[k] = hsp_array[k+1];
            }
            hsp_array[hsp_count] = hsp;
        }
        i += j;
    }

    //if (purge) {
        //Blast_HSPListPurgeNullHSPs(hsp_list);
    //}

    return hsp_count;
}

Int4 get_ncbi_dropoff_score(struct PSSMatrix PSSMatrix)
{
    int compositionBasedStats = 2;
    int cbs_stretch = (compositionBasedStats > 1) ? 5 : 1;
    int evalue = 10; 
    Blast_KarlinBlk kbp;
    Blast_GumbelBlk gbp;

    gbp.Lambda = kbp.Lambda = 0.26700000000000002;
    kbp.K = 0.041000000000000002;
    kbp.logK = -3.1941832122778293;
    kbp.H = 0.14000000000000001;
    kbp.paramC = 0; 

    gbp.C = 0.66971999999999998;

    gbp.G = 12;
    gbp.a = 1.8999999999999999;
    gbp.Alpha = 42.602800000000002;
    gbp.Sigma = 43.636200000000002;
    gbp.a_un = 0.79159999999999997;
    gbp.Alpha_un = 4.9646600000000003; 
    gbp.b = -26.601600000000001;
    gbp.Beta = -903.31536000000006; 
    gbp.Tau = -928.11696000000006; 
    gbp.db_length = readdb_numberOfLetters; 
    gbp.filled = 1;

    int avg_subject_length = 10; 
    int new_cutoff = BLAST_SpougeEtoS(cbs_stretch * evalue,
            &kbp, &gbp,
            PSSMatrix.length, avg_subject_length);

    return new_cutoff;
}

#define DBL_MAX 1.79769e+308
double min_lambda = DBL_MAX;


Int4 get_ncbi_ungappedNominalTrigger(struct PSSMatrix PSSMatrix, struct scoreMatrix scoreMatrix)
{
    int blast_ungappedNominalTrigger_t = 0;
    BlastScoreBlk sbp;
    sbp.protein_alphabet = 1;
    sbp.alphabet_code = 11;
    sbp.alphabet_size = 28;
    sbp.alphabet_start = 0;
    sbp.matrix = scoreMatrix.matrix;
    sbp.loscore = -4;
    sbp.hiscore = 11;
    sbp.penalty = 0;
    sbp.scale_factor = 1;
    sbp.reward = 0;

    sbp.ambiguous_res = NULL;
    sbp.ambig_occupy = 0;
    sbp.ambig_size = 0;

    //BLAST_ScoreSetAmbigRes(&sbp, 'X');
    //BLAST_ScoreSetAmbigRes(&sbp, '*');
    Blast_KarlinBlk kbp;

    //kbp = Blast_KarlinBlkNew();
    Blast_ScoreBlkKbpUngappedCalc(&sbp, PSSMatrix, &blast_ungappedNominalTrigger_t, &kbp);

    Blast_KarlinBlkFree(sbp.kbp_ideal);


    min_lambda = MIN(kbp.Lambda, min_lambda);
    //printf("Lambda: %f K: %f logK: %f H: %f paramC: %f min_lambda: %f\n", kbp.Lambda, kbp.K,
    //kbp.logK, kbp.H, kbp.paramC, min_lambda);
    //printf("blast_ungappedNominalTrigger_t: %d\n", blast_ungappedNominalTrigger_t);

    //Blast_KarlinBlkFree(kbp);
    return blast_ungappedNominalTrigger_t;
}
