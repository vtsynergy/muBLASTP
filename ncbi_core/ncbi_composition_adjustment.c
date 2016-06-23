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



#include "blast.h"
#include <limits.h>


#ifdef __cplusplus
extern "C" {
#endif


/* Some characters in the NCBIstdaa alphabet, including ambiguity
   characters, selenocysteine and the stop character. */
enum { eGapChar = 0, eBchar = 2,  eDchar = 4,  eEchar = 5, eIchar = 9,
    eLchar = 11,  eNchar = 13, eQchar = 15, eXchar = 21,
    eZchar = 23,  eSelenocysteine = 24, eStopChar = 25,
    eOchar = 26,  eJchar = 27};

typedef struct Compo_FrequencyData {
    const char * name;   /**< name of the matrix */
    const double (*joint_probs)[COMPO_NUM_TRUE_AA]; /**< joint probabilities */
    const double * background; /**< background frequencies */
} Compo_FrequencyData;

#define NUM_SUPPORTED_MATRICES 8

static const double
BLOSUM62_JOINT_PROBS[COMPO_NUM_TRUE_AA][COMPO_NUM_TRUE_AA] =
{
  {2.1497573378347484e-02, 2.3470224274721213e-03, 1.9493235258876179e-03,
   2.1674844853066858e-03, 1.5903351423026848e-03, 1.9242657898716525e-03,
   2.9879059292799641e-03, 5.8158526388051033e-03, 1.1076584657559144e-03,
   3.1880644746334580e-03, 4.4186245468471547e-03, 3.3466571942021082e-03,
   1.3412107617355408e-03, 1.6360627863999076e-03, 2.1568959784943114e-03,
   6.2524987419815400e-03, 3.7180506975672363e-03, 4.0281679108936688e-04,
   1.2999956675626666e-03, 5.0679056444508912e-03},
  {2.3470224274721213e-03, 1.7757465118386322e-02, 1.9786027128591904e-03,
   1.5865480081162602e-03, 3.9365984789376245e-04, 2.4858611089731411e-03,
   2.6933867548771758e-03, 1.7221140903704937e-03, 1.2407382229440791e-03,
   1.2435878276496955e-03, 2.4193952633248727e-03, 6.2339060289407083e-03,
   8.0309461712520876e-04, 9.3181986323789834e-04, 9.5783034332718700e-04,
   2.2660898636037261e-03, 1.7802796534180537e-03, 2.6571979312581875e-04,
   9.2634607111251918e-04, 1.5810185245264004e-03},
  {1.9493235258876179e-03, 1.9786027128591904e-03, 1.4140291972553610e-02,
   3.7201973506001745e-03, 4.3845466068066216e-04, 1.5304436972610567e-03,
   2.2097156829738759e-03, 2.8591871815612977e-03, 1.4301072616183181e-03,
   9.9437221166923172e-04, 1.3690958423974782e-03, 2.4402105140841090e-03,
   5.2943633069226512e-04, 7.5004227978192801e-04, 8.6016459857770028e-04,
   3.1466019144814608e-03, 2.2360795375444384e-03, 1.6159545671597605e-04,
   7.0048422794024819e-04, 1.2014015528772706e-03},
  {2.1674844853066858e-03, 1.5865480081162602e-03, 3.7201973506001745e-03,
   2.1274574617480089e-02, 3.9909227141697264e-04, 1.6481246723433428e-03,
   4.9158017471929655e-03, 2.5221102126636373e-03, 9.5384849402143984e-04,
   1.2347404942429857e-03, 1.5202051791453383e-03, 2.4453087721980561e-03,
   4.6429229320514104e-04, 7.6023722413111566e-04, 1.2373315413524663e-03,
   2.8035127901697272e-03, 1.8961512776990257e-03, 1.6218020183662784e-04,
   5.9842263937853702e-04, 1.3158365660538270e-03},
  {1.5903351423026848e-03, 3.9365984789376245e-04, 4.3845466068066216e-04,
   3.9909227141697264e-04, 1.1931352277704348e-02, 3.0937204045913537e-04,
   3.8338775043186374e-04, 7.6951976030099293e-04, 2.2976387481074697e-04,
   1.0956590131781735e-03, 1.5682982157153873e-03, 5.0124929379033781e-04,
   3.7717165634097634e-04, 5.1389991547056834e-04, 3.6111795849154795e-04,
   1.0432626586831986e-03, 9.3041313726939057e-04, 1.4474923964368156e-04,
   3.4603772624580643e-04, 1.3606607271146112e-03},
  {1.9242657898716525e-03, 2.4858611089731411e-03, 1.5304436972610567e-03,
   1.6481246723433428e-03, 3.0937204045913537e-04, 7.3292255467189687e-03,
   3.5385780499965817e-03, 1.3683038039160171e-03, 1.0489026828741754e-03,
   8.9102936026571569e-04, 1.6174411456311808e-03, 3.0968229715707327e-03,
   7.3993258722701268e-04, 5.4255147972143906e-04, 8.4668181752066874e-04,
   1.8931125300036275e-03, 1.3796838284921874e-03, 2.2737931366728891e-04,
   6.7584155312457842e-04, 1.1660966117775285e-03},
  {2.9879059292799641e-03, 2.6933867548771758e-03, 2.2097156829738759e-03,
   4.9158017471929655e-03, 3.8338775043186374e-04, 3.5385780499965817e-03,
   1.6133927472163669e-02, 1.9380952488713059e-03, 1.3667885452189439e-03,
   1.2192061706431622e-03, 2.0030316026648431e-03, 4.1322603720305197e-03,
   6.7909745467514783e-04, 8.5179405867513139e-04, 1.4216207127018586e-03,
   2.9539180653600089e-03, 2.0493063257644955e-03, 2.6488552587183780e-04,
   8.7044186256788659e-04, 1.6987763526262680e-03},
  {5.8158526388051033e-03, 1.7221140903704937e-03, 2.8591871815612977e-03,
   2.5221102126636373e-03, 7.6951976030099293e-04, 1.3683038039160171e-03,
   1.9380952488713059e-03, 3.7804346453413303e-02, 9.5813607255887238e-04,
   1.3849118546156933e-03, 2.0864716056392773e-03, 2.5392537741810947e-03,
   7.3281559749652399e-04, 1.1976708695723554e-03, 1.3641171883713547e-03,
   3.8342830901664762e-03, 2.1858459940987062e-03, 4.0740829083805248e-04,
   8.3467413018106177e-04, 1.8218235950233687e-03},
  {1.1076584657559144e-03, 1.2407382229440791e-03, 1.4301072616183181e-03,
   9.5384849402143984e-04, 2.2976387481074697e-04, 1.0489026828741754e-03,
   1.3667885452189439e-03, 9.5813607255887238e-04, 9.2802502369336622e-03,
   5.8089627083019206e-04, 9.8696608463236094e-04, 1.1873625842258938e-03,
   3.8264639620910225e-04, 8.1041076335565583e-04, 4.7770135861914477e-04,
   1.1052034635193162e-03, 7.4371746073077327e-04, 1.5168037757411286e-04,
   1.5213771111755425e-03, 6.4882907765797669e-04},
  {3.1880644746334580e-03, 1.2435878276496955e-03, 9.9437221166923172e-04,
   1.2347404942429857e-03, 1.0956590131781735e-03, 8.9102936026571569e-04,
   1.2192061706431622e-03, 1.3849118546156933e-03, 5.8089627083019206e-04,
   1.8441526588740136e-02, 1.1382470627796603e-02, 1.5655862274689192e-03,
   2.5081290988482057e-03, 3.0458868657559346e-03, 1.0068164685944146e-03,
   1.7225081689171561e-03, 2.6953622613315018e-03, 3.6183761166072852e-04,
   1.3821121844492116e-03, 1.1972663837662637e-02},
  {4.4186245468471547e-03, 2.4193952633248727e-03, 1.3690958423974782e-03,
   1.5202051791453383e-03, 1.5682982157153873e-03, 1.6174411456311808e-03,
   2.0030316026648431e-03, 2.0864716056392773e-03, 9.8696608463236094e-04,
   1.1382470627796603e-02, 3.7141460156350926e-02, 2.4634345023228079e-03,
   4.9293545515183088e-03, 5.4151301166464015e-03, 1.4146090399381900e-03,
   2.4277107072013821e-03, 3.3238031308707055e-03, 7.3206640617832933e-04,
   2.2096734692836624e-03, 9.4786263030457313e-03},
  {3.3466571942021082e-03, 6.2339060289407083e-03, 2.4402105140841090e-03,
   2.4453087721980561e-03, 5.0124929379033781e-04, 3.0968229715707327e-03,
   4.1322603720305197e-03, 2.5392537741810947e-03, 1.1873625842258938e-03,
   1.5655862274689192e-03, 2.4634345023228079e-03, 1.6113385590544604e-02,
   9.0876633395557617e-04, 9.4875149773685364e-04, 1.5773020912564391e-03,
   3.1016069999481111e-03, 2.3467014804084987e-03, 2.7198500003555514e-04,
   9.9908866586876396e-04, 1.9360424083099779e-03},
  {1.3412107617355408e-03, 8.0309461712520876e-04, 5.2943633069226512e-04,
   4.6429229320514104e-04, 3.7717165634097634e-04, 7.3993258722701268e-04,
   6.7909745467514783e-04, 7.3281559749652399e-04, 3.8264639620910225e-04,
   2.5081290988482057e-03, 4.9293545515183088e-03, 9.0876633395557617e-04,
   4.0477309321969848e-03, 1.1901770463553603e-03, 4.0824445213456919e-04,
   8.5603787638552766e-04, 1.0095451907679563e-03, 1.9872537223131380e-04,
   5.7145288352831449e-04, 2.3123361470140736e-03},
  {1.6360627863999076e-03, 9.3181986323789834e-04, 7.5004227978192801e-04,
   7.6023722413111566e-04, 5.1389991547056834e-04, 5.4255147972143906e-04,
   8.5179405867513139e-04, 1.1976708695723554e-03, 8.1041076335565583e-04,
   3.0458868657559346e-03, 5.4151301166464015e-03, 9.4875149773685364e-04,
   1.1901770463553603e-03, 1.8277684015431908e-02, 5.2528021756783813e-04,
   1.1939618185901600e-03, 1.1624184369750680e-03, 8.4917468952377874e-04,
   4.2392005745634370e-03, 2.5763052227920180e-03},
  {2.1568959784943114e-03, 9.5783034332718700e-04, 8.6016459857770028e-04,
   1.2373315413524663e-03, 3.6111795849154795e-04, 8.4668181752066874e-04,
   1.4216207127018586e-03, 1.3641171883713547e-03, 4.7770135861914477e-04,
   1.0068164685944146e-03, 1.4146090399381900e-03, 1.5773020912564391e-03,
   4.0824445213456919e-04, 5.2528021756783813e-04, 1.9066033679132538e-02,
   1.6662567934883051e-03, 1.3511005665728870e-03, 1.4152209821874487e-04,
   4.5224391125285910e-04, 1.2451325046931832e-03},
  {6.2524987419815400e-03, 2.2660898636037261e-03, 3.1466019144814608e-03,
   2.8035127901697272e-03, 1.0432626586831986e-03, 1.8931125300036275e-03,
   2.9539180653600089e-03, 3.8342830901664762e-03, 1.1052034635193162e-03,
   1.7225081689171561e-03, 2.4277107072013821e-03, 3.1016069999481111e-03,
   8.5603787638552766e-04, 1.1939618185901600e-03, 1.6662567934883051e-03,
   1.2585947097159817e-02, 4.7004857686835334e-03, 2.8731729176487776e-04,
   1.0299846310599138e-03, 2.3587292053265561e-03},
  {3.7180506975672363e-03, 1.7802796534180537e-03, 2.2360795375444384e-03,
   1.8961512776990257e-03, 9.3041313726939057e-04, 1.3796838284921874e-03,
   2.0493063257644955e-03, 2.1858459940987062e-03, 7.4371746073077327e-04,
   2.6953622613315018e-03, 3.3238031308707055e-03, 2.3467014804084987e-03,
   1.0095451907679563e-03, 1.1624184369750680e-03, 1.3511005665728870e-03,
   4.7004857686835334e-03, 1.2514818886617953e-02, 2.8575770858467209e-04,
   9.4161039895612720e-04, 3.6402328079338207e-03},
  {4.0281679108936688e-04, 2.6571979312581875e-04, 1.6159545671597605e-04,
   1.6218020183662784e-04, 1.4474923964368156e-04, 2.2737931366728891e-04,
   2.6488552587183780e-04, 4.0740829083805248e-04, 1.5168037757411286e-04,
   3.6183761166072852e-04, 7.3206640617832933e-04, 2.7198500003555514e-04,
   1.9872537223131380e-04, 8.4917468952377874e-04, 1.4152209821874487e-04,
   2.8731729176487776e-04, 2.8575770858467209e-04, 6.4699301717154852e-03,
   8.8744160259272527e-04, 3.5578318710317554e-04},
  {1.2999956675626666e-03, 9.2634607111251918e-04, 7.0048422794024819e-04,
   5.9842263937853702e-04, 3.4603772624580643e-04, 6.7584155312457842e-04,
   8.7044186256788659e-04, 8.3467413018106177e-04, 1.5213771111755425e-03,
   1.3821121844492116e-03, 2.2096734692836624e-03, 9.9908866586876396e-04,
   5.7145288352831449e-04, 4.2392005745634370e-03, 4.5224391125285910e-04,
   1.0299846310599138e-03, 9.4161039895612720e-04, 8.8744160259272527e-04,
   1.0246100213822419e-02, 1.5489827890922993e-03},
  {5.0679056444508912e-03, 1.5810185245264004e-03, 1.2014015528772706e-03,
   1.3158365660538270e-03, 1.3606607271146112e-03, 1.1660966117775285e-03,
   1.6987763526262680e-03, 1.8218235950233687e-03, 6.4882907765797669e-04,
   1.1972663837662637e-02, 9.4786263030457313e-03, 1.9360424083099779e-03,
   2.3123361470140736e-03, 2.5763052227920180e-03, 1.2451325046931832e-03,
   2.3587292053265561e-03, 3.6402328079338207e-03, 3.5578318710317554e-04,
   1.5489827890922993e-03, 1.9631915140537640e-02}};


static const double
BLOSUM62_bg[COMPO_NUM_TRUE_AA] =
{7.4216205067993410e-02, 5.1614486141284638e-02, 4.4645808512757915e-02,
 5.3626000838554413e-02, 2.4687457167944848e-02, 3.4259650591416023e-02,
 5.4311925684587502e-02, 7.4146941452644999e-02, 2.6212984805266227e-02,
 6.7917367618953756e-02, 9.8907868497150955e-02, 5.8155682303079680e-02,
 2.4990197579643110e-02, 4.7418459742284751e-02, 3.8538003320306206e-02,
 5.7229029476494421e-02, 5.0891364550287033e-02, 1.3029956129972148e-02,
 3.2281512313758580e-02, 7.2919098205619245e-02};


#define MATRIX_DATA(MAT) \
    { #MAT, MAT##_JOINT_PROBS, MAT##_bg }

static Compo_FrequencyData
s_FrequencyData[NUM_SUPPORTED_MATRICES] =
{
    MATRIX_DATA(BLOSUM62),
};

static Compo_FrequencyData *
s_LocateFrequencyData(const char * matrix) 
{
    return &s_FrequencyData[0];
}


const double *
Blast_GetMatrixBackgroundFreq(const char *matrix_name)
{
    Compo_FrequencyData * data = s_LocateFrequencyData(matrix_name);
    if (NULL != data) {
        return data->background;
    } else {
        fprintf(stderr, "matrix %s is not supported "
                "for RE based adjustment\n", matrix_name);
        return NULL;
    }
}

double
Blast_GetRelativeEntropy(const double A[], const double B[])
{
    int i;                 /* loop index over letters */
    double temp;           /* intermediate term */
    double value = 0.0;    /* square of relative entropy */

    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        temp = (A[i] + B[i]) / 2;
        if (temp > 0) {
            if (A[i] > 0) {
                value += A[i] * log(A[i] / temp) / 2;
            }
            if (B[i] > 0) {
                value += B[i] * log(B[i] / temp) / 2;
            }
        }
    }
    if (value < 0) {             /* must be numerical rounding error */
        value = 0;
    }
    return sqrt(value);
}

#define HALF_CIRCLE_DEGREES 180
#define PI 3.1415926543
#define QUERY_MATCH_DISTANCE_THRESHOLD 0.16
#define LENGTH_RATIO_THRESHOLD 3.0
#define ANGLE_DEGREE_THRESHOLD 70.0
#define HIGH_PAIR_THRESHOLD 0.4
#define LENGTH_LOWER_THRESHOLD 50

static int
s_HighPairFrequencies(const double * letterProbs, int length)
{
    int i; /*index*/
    double max, second; /*two highest letter probabilities*/

    if (length <= LENGTH_LOWER_THRESHOLD) {
        return FALSE;
    }
    max = 0;
    second = 0;
    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        if (letterProbs[i] > second) {
            second = letterProbs[i];
            if (letterProbs[i] > max) {
                second = max;
                max = letterProbs[i];
            }
        }
    }
    return (max + second) > HIGH_PAIR_THRESHOLD;
}


static int
s_HighPairEitherSeq(const double * P_query, int length1,
                    const double * P_match, int length2)
{
    int result1, result2;

    result1 = s_HighPairFrequencies(P_query, length1);
    result2 = s_HighPairFrequencies(P_match, length2);

    return result1 || result2;
}

static EMatrixAdjustRule
s_TestToApplyREAdjustmentConditional(int Len_query,
                                     int Len_match,
                                     const double * P_query,
                                     const double * P_match,
                                     const char *matrix_name)
{
    EMatrixAdjustRule which_rule; /* which relative entropy mode to
                                     return */
    int i;                       /* loop indices */
    double p_query[COMPO_NUM_TRUE_AA];
    double p_match[COMPO_NUM_TRUE_AA]; /*letter probabilities
                                                for query and match*/
    const double *p_matrix;       /* letter probabilities used in
                                     constructing matrix name*/
    double D_m_mat, D_q_mat, D_m_q;  /* distances between match and
                                        original between query and
                                        original between match and
                                        query*/
    double corr_factor = 0.0;     /* correlation between how p_query
                                     and p_match deviate from p_matrix
                                     */
    double len_q, len_m;          /* lengths of query and matching
                                     sequence in floating point */
    double len_large, len_small;  /* store the larger and smaller of
                                     len_q and len_m */
    double angle;                 /* angle between query and match
                                     probabilities */

    p_matrix = Blast_GetMatrixBackgroundFreq(matrix_name);

    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        p_query[i] = P_query[i];
        p_match[i] = P_match[i];
        corr_factor +=
            (p_query[i] - p_matrix[i]) * (p_match[i] - p_matrix[i]);
    }
    D_m_mat = Blast_GetRelativeEntropy(p_match, p_matrix);
    D_q_mat = Blast_GetRelativeEntropy(p_query, p_matrix);
    D_m_q   = Blast_GetRelativeEntropy(p_match, p_query);

    angle =
        acos((D_m_mat * D_m_mat + D_q_mat * D_q_mat -
              D_m_q * D_m_q) / 2.0 / D_m_mat / D_q_mat);
    /* convert from radians to degrees */
    angle = angle * HALF_CIRCLE_DEGREES / PI;

    len_q = 1.0 * Len_query;
    len_m = 1.0 * Len_match;
    if (len_q > len_m) {
        len_large = len_q;
        len_small = len_m;
    } else {
        len_large = len_m;
        len_small = len_q;
    }
    if (s_HighPairEitherSeq(P_query, Len_query, P_match, Len_match)) {
        which_rule = eUserSpecifiedRelEntropy;
    } else {
      if ((D_m_q > QUERY_MATCH_DISTANCE_THRESHOLD) &&
        (len_large / len_small > LENGTH_RATIO_THRESHOLD) &&
        (angle > ANGLE_DEGREE_THRESHOLD)) {
        which_rule = eCompoScaleOldMatrix;
      } else {
        which_rule = eUserSpecifiedRelEntropy;
      }
    }
    return which_rule;
}


EMatrixAdjustRule
Blast_ChooseMatrixAdjustRule(int length1,
                             int length2,
                             const double * probArray1,
                             const double * probArray2,
                             const char *matrixName,
                             ECompoAdjustModes composition_adjust_mode)
{
    int testFunctionIndex = (int) composition_adjust_mode;

    return
        s_TestToApplyREAdjustmentConditional(length1,    length2,
                                      probArray1, probArray2, matrixName);
}


static int alphaConvert[COMPO_LARGEST_ALPHABET] =
  {(-1), 0, (-1),  4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15,
   16, 19, 17, (-1), 18, (-1), (-1), (-1), (-1), (-1)};


static void
s_GatherLetterProbs(double * outputLetterProbs,
                    const double * inputLetterProbs, int alphsize)
{
    int c; /*index over characters*/

    for (c = 0;  c < alphsize;  c++) {
        if ((-1) != alphaConvert[c]) {
            outputLetterProbs[alphaConvert[c]] = inputLetterProbs[c];
        }
    }
}


void
Blast_CalcFreqRatios(double ** ratios, int alphsize,
                     double row_prob[], double col_prob[])
{
    int i, j;
    for (i = 0;  i < alphsize;  i++) {
        if (row_prob[i] > 0) {
            for (j = 0;  j < alphsize;  j++) {
                if (col_prob[j] > 0) {
                    ratios[i][j] /= (row_prob[i] * col_prob[j]);
                }
            }
        }
    }
}

#define COMPO_SCORE_MIN INT2_MIN

void
Blast_CalcLambdaFullPrecision(double * plambda, int *piterations,
                              double **score, int alphsize,
                              const double row_prob[],
                              const double col_prob[],
                              double lambda_tolerance,
                              double function_tolerance,
                              int max_iterations)
{
    double f = 4;               /* The current function value; initially
                                   set to a value greater than any possible
                                   real value of f */
    double left = 0, right = 1; /* (left, right) is an interval containing
                                   a solution */
    double x = 0.367879441171;  /* The current iterate; initially exp(-1) */
    int is_newton = 0;          /* true if the last iteration was a Newton
                                   step; initially false */
    int i, j, k;                /* iteration indices */
    /* maximum score that occurs with nonzero probability */
    double max_score = COMPO_SCORE_MIN;
    /* average score */
    double avg_score = 0.0;

    /* Find the maximum score with nonzero probability */
    for (i = 0;  i < alphsize;  i++) {
        if (row_prob[i] == 0.0) {
            continue;
        }
        for (j = 0;  j < alphsize;  j++) {
            if (col_prob[j] == 0.0) {
                continue;
            }
            if (max_score < score[i][j]) {
                max_score = score[i][j];
            }
            avg_score += row_prob[i] * col_prob[j] * score[i][j];
        }
    }
    if (max_score <= 0.0 || avg_score >= 0) { 
        /* The iteration cannot converge if max_score is nonpositive
         * or the average score is nonnegative; lambda doesn't exist */
        *piterations = max_iterations;
        *plambda = -1.0;
        return;
    }
    for (k = 0;  k < max_iterations;  k++) {
        double slope;               /* slope of f at x */
        double fold = f;            /* previous value of f */
        double x_pow_max_score;     /* x raised to the power max_score */
        double lambda = -log(x);    /* the iterate lambda, see above */
        int was_newton = is_newton; /* true if the previous iteration
                                       was a Newton step; instead of a
                                       bisection step */
        /* Evaluate the function and its derivative */
        x_pow_max_score = exp(-max_score * lambda);
        f = -x_pow_max_score;
        slope = max_score * f / x;
        for (i = 0;  i < alphsize;  i++) {
            if (row_prob[i] == 0.0) {
                continue;
            }
            for (j = 0;  j < alphsize;  j++) {
                double ff;  /* a term in the sum used to compute f */

               if (col_prob[j] == 0.0) {
                   continue;
               }
               if (max_score != score[i][j]) {
                   double diff_score = max_score - score[i][j];

                   ff = row_prob[i] * col_prob[j] * exp(-lambda * diff_score);
                   slope += diff_score * ff / x;
               } else {
                   ff = row_prob[i] * col_prob[j];
               }
               f += ff;
            }
        }
        /* Finished evaluating the function and its derivative */
        if (f > 0) {
            left = x; /* move the left endpoint */
        } else if (f < 0) {
            right = x; /* move the right endpoint */
        } else { /* f == 0 */
            break; /* x is an exact solution */
        }
        if (right - left <= 2 * left * (1 - right) * lambda_tolerance &&
            fabs(f/x_pow_max_score) <= function_tolerance) {
            /* The midpoint of the interval converged */
            x = (left + right) / 2;
            break;
        }
        if ((was_newton && fabs(f) > .9 * fabs(fold))
            /* if the previous iteration was a Newton step but didn't
             * decrease f sufficiently; or */
             || slope >= 0
             /* if a Newton step will move us away from the desired solution */
        ) {/* then */
            x = (left + right)/2;  /* bisect */
        } else {
            double p = -f/slope;   /* The Newton step */
            double y = x + p;      /* The proposed next iterate */
            if (y <= left || y >= right) { /* The proposed iterate is
                                              not in (left,right) */
                x = (left + right)/2;  /* bisect */
            } else {/* The proposed iterate is in (left,right). Accept it. */
                is_newton = 1;
                x = y;
                if (fabs(p) <= lambda_tolerance * x * (1-x) &&
                    fabs(f/x_pow_max_score) <= function_tolerance) break;
            }
        }
    }  /* End for all iterations k */
    *plambda = (k < max_iterations) ? -log(x) : -1.0;
    *piterations = k;
} 

/** iteration limit for Newton's method for computing Lambda */
static const int kLambdaIterationLimit = 100;
/** bound on error for estimating lambda */
static const double kLambdaErrorTolerance = 0.0000001;
/** bound on the difference between the expected value of
 *  exp(lambda x S) and 1 */
static const double kLambdaFunctionTolerance = 0.00001;

double
Blast_MatrixEntropy(double ** matrix, int alphsize, const double row_prob[],
                    const double col_prob[], double Lambda)
{
    int i, j;
    double entropy = 0.0;
    for (i = 0;  i < alphsize;  i++) {
        for (j = 0;  j < alphsize;  j++) {
            /* the score at (i,j), rescaled to nats */
            double nat_score = Lambda * matrix[i][j];
            entropy += nat_score * exp(nat_score) * row_prob[i] * col_prob[j];
        }
    }
    return entropy;
}

//void
//Nlm_DenseMatrixFree(double *** mat)
//{
    //if(*mat != NULL) {
        //free((*mat)[0]);
        //free(*mat);
    //}
    //*mat = NULL;
//}

int
Blast_EntropyOldFreqNewContext(double * entropy,
                               double * Lambda,
                               int * iter_count,
                               double ** target_freq,
                               const double row_prob[],
                               const double col_prob[])
{
    /* iteration indices */
    int i, j;
    /* Status flag; will be set to zero on success */
    int status = 1;
    /* A matrix of scores in the context consistent with the target 
     * frequencies */
    double  **scores;
    /* Row and column probabilities consistent with the target
     * frequencies; the old context */
    double old_col_prob[COMPO_NUM_TRUE_AA] = {0.0,};
    double old_row_prob[COMPO_NUM_TRUE_AA] = {0.0,};

    *entropy = 0;
    status = 1;

    /* Calculate the matrix "scores" from the target frequencies */
    scores = Nlm_DenseMatrixNew(COMPO_NUM_TRUE_AA, COMPO_NUM_TRUE_AA);
    if (scores == NULL) {
        return -1;
    }
    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        for (j = 0;  j < COMPO_NUM_TRUE_AA; j++) {
            old_row_prob[i] += target_freq[i][j];
            old_col_prob[j] += target_freq[i][j];
        }
    }
    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        memcpy(scores[i], target_freq[i], COMPO_NUM_TRUE_AA * sizeof(double));
    }
    Blast_CalcFreqRatios(scores, COMPO_NUM_TRUE_AA,
                         old_row_prob, old_col_prob);
    Blast_FreqRatioToScore(scores, COMPO_NUM_TRUE_AA, COMPO_NUM_TRUE_AA, 1.0);
    /* Finished calculating the matrix "scores" */

    Blast_CalcLambdaFullPrecision(Lambda, iter_count, scores,
                                  COMPO_NUM_TRUE_AA, row_prob,
                                  col_prob, kLambdaErrorTolerance,
                                  kLambdaFunctionTolerance,
                                  kLambdaIterationLimit);
    if (*iter_count <  kLambdaIterationLimit) {
        *entropy = Blast_MatrixEntropy(scores, COMPO_NUM_TRUE_AA,
                                       row_prob, col_prob, *Lambda);
        status = 0;
    }
    Nlm_DenseMatrixFree(&scores);
    return status;
}

double
Blast_TargetFreqEntropy(double ** target_freq)
{
    int i, j;          /* iteration indices */
    double entropy;    /* the entropy to be returned */
    /* Row probabilities consistent with the target frequencies */
    double row_prob[COMPO_NUM_TRUE_AA] = {0,};
    double col_prob[COMPO_NUM_TRUE_AA] = {0,};

    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        for (j = 0;  j < COMPO_NUM_TRUE_AA;  j++) {
            row_prob[i] += target_freq[i][j];
            col_prob[j] += target_freq[i][j];
        }
    }
    entropy = 0.0;
    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        for (j = 0;  j < COMPO_NUM_TRUE_AA;  j++) {
            double freq = target_freq[i][j];
            entropy += freq * log(freq / row_prob[i] / col_prob[j]);
        }
    }
    return entropy;
}

void
Blast_ApplyPseudocounts(double * probs20,
                    int number_of_observations,
                    const double * background_probs20,
                    int pseudocounts)
{
    int i;                 /* loop index */
    double weight;         /* weight assigned to pseudocounts */
    double sum;            /* sum of the observed frequencies */
    /* pseudocounts as a double */
    double dpseudocounts = pseudocounts;

    /* Normalize probabilities */
    sum = 0.0;
    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        sum += probs20[i];
    }
    if (sum == 0.0) {  /* Can't normalize a zero vector */
        sum = 1.0;
    }
    weight = dpseudocounts / (number_of_observations + dpseudocounts);
    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        probs20[i] = (1.0 - weight) * probs20[i] / sum
            + weight * background_probs20[i];
    }
}

double **
Nlm_LtriangMatrixNew(int n)
{
    int i;                      /* iteration index */
    double ** L;                /* the new, lower triangular matrix */
    size_t nelts;               /* the number of elements in
                                   the matrix */
    nelts = ((size_t) n * (n + 1))/2;

    L    = (double**) calloc(n, sizeof(double *));
    if (L != NULL) {
        L[0] = (double*) malloc(nelts * sizeof(double));
        if (L[0] != NULL) {
            for (i = 1;  i < n;  i++) {
                L[i] = L[i - 1] + i;
            }
        } else {
            free(L);
            L = NULL;
        }
    }
    return L;
}

void
ReNewtonSystemFree(ReNewtonSystem ** newton_system)
{
    if (*newton_system != NULL) {
        Nlm_DenseMatrixFree(&(*newton_system)->W);

        free((*newton_system)->Dinv);
        (*newton_system)->Dinv = NULL;

        free((*newton_system)->grad_re);
        (*newton_system)->grad_re = NULL;

        free(*newton_system);
        *newton_system = NULL;
    }
}

ReNewtonSystem * ReNewtonSystemNew(int alphsize)
{
    ReNewtonSystem * newton_system;  /* the new ReNewtonSystem */

    newton_system = (ReNewtonSystem *) malloc(sizeof(ReNewtonSystem));
    if (newton_system != NULL) {
        newton_system->alphsize = alphsize;
        newton_system->constrain_rel_entropy = 1;
        newton_system->W = NULL;
        newton_system->Dinv = NULL;
        newton_system->grad_re = NULL;

        newton_system->W = Nlm_LtriangMatrixNew(2 * alphsize);
        if (newton_system->W == NULL)
            goto error_return;
        newton_system->Dinv =
            (double *) malloc(alphsize * alphsize * sizeof(double));
        if (newton_system->Dinv == NULL)
            goto error_return;
        newton_system->grad_re =
            (double *) malloc(alphsize * alphsize * sizeof(double));
        if (newton_system->grad_re == NULL)
            goto error_return;
    }
    goto normal_return;
error_return:
    ReNewtonSystemFree(&newton_system);
normal_return:

    return newton_system;
}

static void
ComputeScoresFromProbs(double scores[],
                       int alphsize,
                       const double target_freqs[],
                       const double row_freqs[],
                       const double col_freqs[])
{
    int i, j;     /* iteration indices over characters in the alphabet */
    int k;        /* index into scores and target_freqs */

    for (i = 0;  i < alphsize;  i++) {
        for (j = 0;  j < alphsize;  j++) {
            k = i * alphsize + j;

            scores[k] = log(target_freqs[k] / (row_freqs[i] * col_freqs[j]));
        }
    }
}

static void
EvaluateReFunctions(double values[], double ** grads, int alphsize,
                    const double x[], const double q[],
                    const double scores[],
                    int constrain_rel_entropy)
{
    int k;         /* iteration index over elements of x, q and scores */
    double temp;   /* holds intermediate values in a computation */

    values[0] = 0.0; values[1] = 0.0;
    for (k = 0;  k < alphsize * alphsize;  k++) {
        temp = log(x[k] / q[k]);

        values[0]   += x[k] * temp;
        grads[0][k]  = temp + 1;

        if (constrain_rel_entropy) {
            temp += scores[k];

            values[1]   += x[k] * temp;
            grads[1][k]  = temp + 1;
        }
    }
}

double
Nlm_EuclideanNorm(const double v[], int n)
{
    double sum   = 1.0;   /* sum of squares of elements in v */
    double scale = 0.0;   /* a scale factor for the elements in v */
    int i;                /* iteration index */

    for (i = 0;  i < n;  i++) {
        if (v[i] != 0.0) {
            double absvi = fabs(v[i]);
            if (scale < absvi) {
                sum = 1.0 + sum * (scale/absvi) * (scale/absvi);
                scale = absvi;
            } else {
                sum += (absvi/scale) * (absvi/scale);
            }
        }
    }
    return scale * sqrt(sum);
}

static void
MultiplyByAtranspose(double beta, double y[], int alphsize,
                     double alpha, const double x[])
{
    int i, j;     /* iteration indices over characters in the alphabet */
    int k;        /* index of a row of A transpose (a column of A); also
                      an index into y */

    if (beta == 0.0) {
        /* Initialize y to zero, without reading any elements from y */
        for (k = 0;  k < alphsize * alphsize;  k++) {
            y[k] = 0.0;
        }
    } else if (beta != 1.0) {
        /* rescale y */
        for (k = 0;  k < alphsize * alphsize;  k++) {
            y[k] *= beta;
        }
    }
    for (i = 0;  i < alphsize;  i++) {
        for (j = 0;  j < alphsize;  j++) {
            k = i * alphsize + j;

            y[k] += alpha * x[j];
            if (i > 0) {
                y[k] += alpha * x[i + alphsize - 1];
            }
        }
    }
}


static void
DualResiduals(double resids_x[], int alphsize, double ** grads,
              const double z[], int constrain_rel_entropy)
{
    int i;                        /* iteration index */
    int n = alphsize * alphsize;  /* size of resids_x */

    if (constrain_rel_entropy) {
        double eta;                /* dual variable for the relative
                                      entropy constraint */
        eta = z[2 * alphsize - 1];
        for (i = 0;  i < n;  i++) {
            resids_x[i] = -grads[0][i] + eta * grads[1][i];
        }
    } else {
        for (i = 0;  i < n;  i++) {
            resids_x[i] = -grads[0][i];
        }
    }
    MultiplyByAtranspose(1.0, resids_x, alphsize, 1.0, z);
}

static void
MultiplyByA(double beta, double y[], int alphsize,
            double alpha, const double x[])
{
    int i, j;     /* iteration indices over characters in the alphabet */
    if (beta == 0.0) {
        /* Initialize y to zero, without reading any elements from y */
        for (i = 0;  i < 2 * alphsize - 1;  i++) {
            y[i] = 0.0;
        }
    } else if (beta != 1.0) {
        /* rescale y */
        for (i = 0;  i < 2 * alphsize - 1;  i++) {
            y[i] *= beta;
        }
    }
    for (i = 0;  i < alphsize;  i++) {
        for (j = 0;  j < alphsize;  j++) {
            y[j] += alpha * x[i * alphsize + j];
        }
    }
    for (i = 1;  i < alphsize;  i++) {
        for (j = 0;  j < alphsize;  j++) {
            y[i + alphsize - 1] += alpha * x[i * alphsize + j];
        }
    }
}

static void
ResidualsLinearConstraints(double rA[], int alphsize, const double x[],
                           const double row_sums[], const double col_sums[])
{
    int i;             /* iteration index */

    for (i = 0;  i < alphsize;  i++) {
        rA[i] = col_sums[i];
    }
    for (i = 1;  i < alphsize;  i++) {
        rA[i + alphsize - 1] = row_sums[i];
    }
    MultiplyByA(1.0, rA, alphsize, -1.0, x);
}

static void
CalculateResiduals(double * rnorm,
                   double resids_x[],
                   int alphsize,
                   double resids_z[],
                   const double values[],
                   double ** grads,
                   const double row_sums[],
                   const double col_sums[],
                   const double x[],
                   const double z[],
                   int constrain_rel_entropy,
                   double relative_entropy)
{
    /* Euclidean norms of the primal and dual residuals */
    double norm_resids_z, norm_resids_x;

    DualResiduals(resids_x, alphsize, grads, z, constrain_rel_entropy);
    norm_resids_x = Nlm_EuclideanNorm(resids_x, alphsize * alphsize);

    ResidualsLinearConstraints(resids_z, alphsize, x, row_sums, col_sums);

    if (constrain_rel_entropy) {
        resids_z[2 * alphsize - 1] = relative_entropy - values[1];

        norm_resids_z = Nlm_EuclideanNorm(resids_z, 2 * alphsize);
    } else {
        norm_resids_z = Nlm_EuclideanNorm(resids_z, 2 * alphsize - 1);
    }
    *rnorm =
        sqrt(norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z);
}

static void
ScaledSymmetricProductA(double ** W, const double diagonal[], int alphsize)
{
    int rowW, colW;   /* iteration indices over the rows and columns of W */
    int i, j;         /* iteration indices over characters in the alphabet */
    int m;            /* The number of rows in A; also the size of W */

    m = 2 * alphsize - 1;

    for (rowW = 0;  rowW < m;  rowW++) {
        for (colW = 0;  colW <= rowW;  colW++) {
            W[rowW][colW] = 0.0;
        }
    }
    for (i = 0;  i < alphsize;  i++) {
        for (j = 0;  j < alphsize;  j++) {
            double dd;     /* an individual diagonal element */

            dd = diagonal[i * alphsize + j];

            W[j][j] += dd;
            if (i > 0) {
                W[i + alphsize - 1][j] += dd;
                W[i + alphsize - 1][i + alphsize - 1] += dd;
            }
        }
    }
}

void
Nlm_FactorLtriangPosDef(double ** A, int n)
{
    int i, j, k;                /* iteration indices */
    double temp;                /* temporary variable for intermediate
                                   values in a computation */

    for (i = 0;  i < n;  i++) {
        for (j = 0;  j < i;  j++) {
            temp = A[i][j];
            for (k = 0;  k < j;  k++) {
                temp -= A[i][k] * A[j][k];
            }
            A[i][j] = temp/A[j][j];
        }
        temp = A[i][i];
        for (k = 0;  k < i;  k++) {
            temp -= A[i][k] * A[i][k];
        }
        A[i][i] = sqrt(temp);
    }
}


static void
FactorReNewtonSystem(ReNewtonSystem * newton_system,
                     const double x[],
                     const double z[],
                     double ** grads,
                     int constrain_rel_entropy,
                     double * workspace)
{
    int i;          /* iteration index */
    int n;          /* the length of x */
    int m;          /* the length of z */

    /* Pointers to fields in newton_systems; the names of the local
     * variables match the names of the fields. */
    double ** W       = newton_system->W;
    int alphsize     = newton_system->alphsize;
    double * Dinv     = newton_system->Dinv;
    double * grad_re  = newton_system->grad_re;

    n = alphsize * alphsize;
    m = constrain_rel_entropy ? 2 * alphsize : 2 * alphsize - 1;

    newton_system->constrain_rel_entropy = constrain_rel_entropy;

    /* The original system has the form
     *
     *     (D     J^T)
     *     (J     0  ).
     *
     * We block reduce the system to
     *
     *     (D    J^T          )
     *     (0    -J D^{-1} J^T).
     *
     * First we find the inverse of the diagonal matrix D. */

     if (constrain_rel_entropy) {
        double eta;             /* dual variable for the relative
                                   entropy constraint */
        eta = z[m - 1];
        for (i = 0;  i < n;  i++) {
            Dinv[i] = x[i] / (1 - eta);
        }
    } else {
        memcpy(Dinv, x, n * sizeof(double));
    }

    /* Then we compute J D^{-1} J^T; First fill in the part that corresponds
     * to the linear constraints */
    ScaledSymmetricProductA(W, Dinv, alphsize);

    if (constrain_rel_entropy) {
        /* Save the gradient of the relative entropy constraint. */
        memcpy(grad_re, grads[1], n * sizeof(double));

        /* Fill in the part of J D^{-1} J^T that corresponds to the relative
         * entropy constraint. */
        W[m - 1][m - 1] = 0.0;
        for (i = 0;  i < n;  i++) {
            workspace[i] = Dinv[i] * grad_re[i];

            W[m - 1][m - 1] += grad_re[i] * workspace[i];
        }
        MultiplyByA(0.0, &W[m - 1][0], alphsize, 1.0, workspace);
    }
    /* Factor J D^{-1} J^T and save the result in W. */
    Nlm_FactorLtriangPosDef(W, m);
}


/* Documented in nlm_linear_algebra.h. */
void Nlm_SolveLtriangPosDef(double x[], int n,
                            double ** L )
{
    int i, j;                   /* iteration indices */
    double temp;                /* temporary variable for intermediate
                                   values in a computation */

    /* At point x = b in the equation L L\T y = b */

    /* Forward solve; L z = b */
    for (i = 0;  i < n;  i++) {
        temp = x[i];
        for (j = 0;  j < i;  j++) {
            temp -= L[i][j] * x[j];
        }
        x[i] = temp/L[i][i];
    }
    /* Now x = z.  Back solve the system L\T y = z */
    for (j = n - 1;  j >= 0;  j--) {
        x[j] /= L[j][j];
        for (i = 0;  i < j;  i++) {
            x[i] -= L[j][i] * x[j];
        }
    }
    /* Now x = y, the solution to  L L\T y = b */
}

static void
SolveReNewtonSystem(double x[], double z[],
                    const ReNewtonSystem * newton_system, double workspace[])
{
    int i;                     /* iteration index */
    int n;                     /* the size of x */
    int mA;                    /* the number of linear constraints */
    int m;                     /* the size of z */

    /* Local variables that represent fields of newton_system */
    double ** W       = newton_system->W;
    double * Dinv     = newton_system->Dinv;
    double * grad_re  = newton_system->grad_re;
    int alphsize               = newton_system->alphsize;
    int constrain_rel_entropy  = newton_system->constrain_rel_entropy;

    n  = alphsize * alphsize;
    mA = 2 * alphsize - 1;
    m  = constrain_rel_entropy ? mA + 1 : mA;

    /* Apply the same block reduction to the right-hand side as was
     * applied to the matrix:
     *
     *     rzhat = rz - J D^{-1} rx
     */
    for (i = 0;  i < n;  i++) {
        workspace[i] = x[i] * Dinv[i];
    }
    MultiplyByA(1.0, z, alphsize, -1.0, workspace);

    if (constrain_rel_entropy) {
        for (i = 0;  i < n;  i++) {
            z[m - 1] -= grad_re[i] * workspace[i];
        }
    }

    /* Solve for step in z, using the inverse of J D^{-1} J^T */
    Nlm_SolveLtriangPosDef(z, m, W);

    /* Backsolve for the step in x, using the newly-computed step in z.
     *
     *     x = D^{-1) (rx + J\T z)
     */
    if (constrain_rel_entropy) {
        for(i = 0; i < n; i++) {
            x[i] += grad_re[i] * z[m - 1];
        }
    }
    MultiplyByAtranspose(1.0, x, alphsize, 1.0, z);

    for (i = 0;  i < n;  i++) {
        x[i] *= Dinv[i];
    }
}

double
Nlm_StepBound(const double x[], int n, const double step_x[], double max)
{
    int i;                 /* iteration index */
    double alpha = max;    /* current largest permitted step */

    for (i = 0; i < n; i++) {
        double alpha_i;    /* a step to the boundary for the current i */

        alpha_i = -x[i] / step_x[i];
        if (alpha_i >= 0 && alpha_i < alpha) {
            alpha = alpha_i;
        }
    }
    return alpha;
}

void Nlm_AddVectors(double y[], int n, double alpha, const double x[])
{
    int i;                     /* iteration index */

    for (i = 0; i < n; i++) {
        y[i] += alpha * x[i];
    }
}

int
Blast_OptimizeTargetFrequencies2(double x[],
                                int alphsize,
                                int *iterations,
                                const double q[],
                                const double row_sums[],
                                const double col_sums[],
                                int constrain_rel_entropy,
                                double relative_entropy,
                                double tol,
                                int maxits,
                                ReNewtonSystem *newton_system,
                                double *z, double *resids_x, double *resids_z, 
                                double *old_scores, double *workspace, double **grads)
{
    int its;       /* number of iterations that have been performed */
    int n;         /* number of target frequencies; the size of x */
    int mA;        /* number of linear constraints */
    int m;         /* total number of constraints */

    double         values[2];   /* values of the nonlinear functions
                                   at this iterate */
    //double ** grads = NULL;     /* gradients of the nonlinear
    //functions at this iterate */

    //ReNewtonSystem *
    //newton_system = NULL;   /* factored matrix of the linear
    //system to be solved at this
    //iteration */
    //double * z = NULL;          /* dual variables (Lagrange multipliers) */
    //double * resids_x = NULL;   /* dual residuals (gradient of Lagrangian) */
    //double * resids_z = NULL;   /* primal (constraint) residuals */
    double rnorm;               /* norm of the residuals for the
                                   current iterate */
    //double * old_scores = NULL; /* a scoring matrix, with lambda = 1,
    //generated from q, row_sums and
    //col_sums */
    //double * workspace = NULL;  /* A vector for intermediate computations */
    int converged;              /* true if Newton's method converged
                                   to a *minimizer* (strong
                                   second-order point) */
    int status;                 /* the return status */
    n  = alphsize * alphsize;
    mA = 2 * alphsize - 1;
    m  = constrain_rel_entropy ? mA + 1 : mA;

    //newton_system = ReNewtonSystemNew(alphsize);
    //if (newton_system == NULL) goto error_return;
    //resids_x = (double *) malloc(n * sizeof(double));
    //if (resids_x == NULL) goto error_return;
    //resids_z = (double *) malloc((mA + 1) * sizeof(double));
    //if (resids_z == NULL) goto error_return;
    ///* z must be initialized to zero */
    //z = (double *) calloc( mA + 1,   sizeof(double));
    //if (z == NULL) goto error_return;
    //old_scores = (double *) malloc(n * sizeof(double));
    //if (old_scores == NULL) goto error_return;
    //workspace = (double *) malloc(n * sizeof(double));
    //if (workspace == NULL) goto error_return;
    //grads = Nlm_DenseMatrixNew(2, n);
    //if (grads == NULL) goto error_return;

    ComputeScoresFromProbs(old_scores, alphsize, q, row_sums, col_sums);

    /* Use q as the initial value for x */
    memcpy(x, q, n * sizeof(double));
    its = 0;        /* Initialize the iteration count. Note that we may
                       converge in zero iterations if the initial x is
                       optimal. */
    while (its <= maxits) {
        /* Compute the residuals */
        EvaluateReFunctions(values, grads, alphsize, x, q, old_scores,
                            constrain_rel_entropy);
        CalculateResiduals(&rnorm, resids_x, alphsize, resids_z, values,
                           grads, row_sums, col_sums, x, z,
                           constrain_rel_entropy, relative_entropy);

        //fprintf(stderr, "it: %d rnorm: %lf\n", its, rnorm);
        /* and check convergence; the test correctly handles the case
           in which rnorm is NaN (not a number). */
        if ( !(rnorm > tol) ) {
            /* We converged at the current iterate */
            break;
        } else {
            /* we did not converge, so increment the iteration counter
               and start a new iteration */
            if (++its <= maxits) {
                /* We have not exceeded the maximum number of iterations;
                   take a Newton step. */
                double alpha;       /* a positive number used to scale the
                                       Newton step. */

                FactorReNewtonSystem(newton_system, x, z, grads,
                                     constrain_rel_entropy, workspace);
                SolveReNewtonSystem(resids_x, resids_z, newton_system,
                                    workspace);

                /* Calculate a value of alpha that ensure that x is
                   positive */
                alpha = Nlm_StepBound(x, n, resids_x, 1.0 / .95);
                alpha *= 0.95;

                Nlm_AddVectors(x, n, alpha, resids_x);
                Nlm_AddVectors(z, m, alpha, resids_z);
            }
        }
    }
    converged = 0;
    if (its <= maxits && rnorm <= tol) {
        /* Newton's iteration converged */
        if ( !constrain_rel_entropy || z[m - 1] < 1 ) {
            /* and the final iterate is a minimizer */
            converged = 1;
        }
    }
    status = converged ? 0 : 1;
    *iterations = its;
    goto normal_return;

error_return:
    status = -1;
    *iterations = 0;
normal_return:

    //Nlm_DenseMatrixFree(&grads);
    //free(workspace);
    //free(old_scores);
    //free(z);
    //free(resids_z);
    //free(resids_x);
    //ReNewtonSystemFree(&newton_system);

    return status;
}

int
Blast_OptimizeTargetFrequencies(double x[],
                                int alphsize,
                                int *iterations,
                                const double q[],
                                const double row_sums[],
                                const double col_sums[],
                                int constrain_rel_entropy,
                                double relative_entropy,
                                double tol,
                                int maxits)
{
    int its;       /* number of iterations that have been performed */
    int n;         /* number of target frequencies; the size of x */
    int mA;        /* number of linear constraints */
    int m;         /* total number of constraints */

    double         values[2];   /* values of the nonlinear functions
                                   at this iterate */
    double ** grads = NULL;     /* gradients of the nonlinear
                                   functions at this iterate */

    ReNewtonSystem *
        newton_system = NULL;   /* factored matrix of the linear
                                   system to be solved at this
                                   iteration */
    double * z = NULL;          /* dual variables (Lagrange multipliers) */
    double * resids_x = NULL;   /* dual residuals (gradient of Lagrangian) */
    double * resids_z = NULL;   /* primal (constraint) residuals */
    double rnorm;               /* norm of the residuals for the
                                   current iterate */
    double * old_scores = NULL; /* a scoring matrix, with lambda = 1,
                                   generated from q, row_sums and
                                   col_sums */
    double * workspace = NULL;  /* A vector for intermediate computations */
    int converged;              /* true if Newton's method converged
                                   to a *minimizer* (strong
                                   second-order point) */
    int status;                 /* the return status */
    n  = alphsize * alphsize;
    mA = 2 * alphsize - 1;
    m  = constrain_rel_entropy ? mA + 1 : mA;

    newton_system = ReNewtonSystemNew(alphsize);
    if (newton_system == NULL) goto error_return;
    resids_x = (double *) malloc(n * sizeof(double));
    if (resids_x == NULL) goto error_return;
    resids_z = (double *) malloc((mA + 1) * sizeof(double));
    if (resids_z == NULL) goto error_return;
    /* z must be initialized to zero */
    z = (double *) calloc( mA + 1,   sizeof(double));
    if (z == NULL) goto error_return;
    old_scores = (double *) malloc(n * sizeof(double));
    if (old_scores == NULL) goto error_return;
    workspace = (double *) malloc(n * sizeof(double));
    if (workspace == NULL) goto error_return;
    grads = Nlm_DenseMatrixNew(2, n);
    if (grads == NULL) goto error_return;

    ComputeScoresFromProbs(old_scores, alphsize, q, row_sums, col_sums);

    /* Use q as the initial value for x */
    memcpy(x, q, n * sizeof(double));
    its = 0;        /* Initialize the iteration count. Note that we may
                       converge in zero iterations if the initial x is
                       optimal. */
    while (its <= maxits) {
        /* Compute the residuals */
        EvaluateReFunctions(values, grads, alphsize, x, q, old_scores,
                            constrain_rel_entropy);
        CalculateResiduals(&rnorm, resids_x, alphsize, resids_z, values,
                           grads, row_sums, col_sums, x, z,
                           constrain_rel_entropy, relative_entropy);

        //fprintf(stderr, "it: %d rnorm: %lf\n", its, rnorm);
        /* and check convergence; the test correctly handles the case
           in which rnorm is NaN (not a number). */
        if ( !(rnorm > tol) ) {
            /* We converged at the current iterate */
            break;
        } else {
            /* we did not converge, so increment the iteration counter
               and start a new iteration */
            if (++its <= maxits) {
                /* We have not exceeded the maximum number of iterations;
                   take a Newton step. */
                double alpha;       /* a positive number used to scale the
                                       Newton step. */

                FactorReNewtonSystem(newton_system, x, z, grads,
                                     constrain_rel_entropy, workspace);
                SolveReNewtonSystem(resids_x, resids_z, newton_system,
                                    workspace);

                /* Calculate a value of alpha that ensure that x is
                   positive */
                alpha = Nlm_StepBound(x, n, resids_x, 1.0 / .95);
                alpha *= 0.95;

                Nlm_AddVectors(x, n, alpha, resids_x);
                Nlm_AddVectors(z, m, alpha, resids_z);
            }
        }
    }
    converged = 0;
    if (its <= maxits && rnorm <= tol) {
        /* Newton's iteration converged */
        if ( !constrain_rel_entropy || z[m - 1] < 1 ) {
            /* and the final iterate is a minimizer */
            converged = 1;
        }
    }
    status = converged ? 0 : 1;
    *iterations = its;
    goto normal_return;

error_return:
    status = -1;
    *iterations = 0;
normal_return:

    //Nlm_DenseMatrixFree(&grads);
    //free(workspace);
    //free(old_scores);
    //free(z);
    //free(resids_z);
    //free(resids_x);
    //ReNewtonSystemFree(&newton_system);

    return status;
}

static void
s_UnpackLetterProbs(double std_probs[], int alphsize, const double probs[])
{
    int c; /*index over characters*/

    for (c = 0;  c < alphsize;  c++) {
        if ((-1) != alphaConvert[c]) {
            std_probs[c] = probs[alphaConvert[c]];
        } else {
            std_probs[c] = 0.0;
        }
    }
}

static void
s_SetPairAmbigProbsToSum(double probs[], int alphsize)
{
    probs[eBchar] = probs[eDchar] + probs[eNchar];
    probs[eZchar] = probs[eEchar] + probs[eQchar];
    if (alphsize > eJchar) {
        probs[eJchar] = probs[eIchar] + probs[eLchar];
    }
}

void
Blast_TrueAaToStdTargetFreqs(double ** StdFreq, int StdAlphsize,
                             double ** freq)
{
    /* Note I'm using a rough convention for this routine that uppercase
     * letters refer to quantities in the standard (larger) alphabet
     * and lowercase letters refer to the true amino acid (smaller)
     * alphabet.
     */
    /* Shorter names for the sizes of the two alphabets */
    const int small_alphsize = COMPO_NUM_TRUE_AA;
    int A, B;          /* characters in the std (big) alphabet */
    int a, b;          /* characters in the small alphabet */
    double sum;        /* sum of values in target_freq; used to normalize */
    sum = 0.0;
    for (a = 0;  a < small_alphsize;  a++) {
        for (b = 0;  b < small_alphsize;  b++) {
            sum +=  freq[a][b];
        }
    }
    for (A = 0;  A < StdAlphsize;  A++) {
        /* for all rows */
        if (alphaConvert[A] < 0) {
            /* the row corresponds to a nonstandard reside */
            for (B = 0;  B < StdAlphsize;  B++) {
                StdFreq[A][B] = 0.0;
            }
        } else {
            /* the row corresponds to a standard reside */
            a = alphaConvert[A];

            for (B = 0;  B < StdAlphsize;  B++) {
                /* for all columns */
                if (alphaConvert[B] < 0) {
                    /* the column corresponds to a nonstandard reside */
                    StdFreq[A][B] = 0.0;
                } else {
                    /* the column corresponds to a standard reside */
                    b = alphaConvert[B];
                    StdFreq[A][B] = freq[a][b] / sum;
                }
            }
            /* Set values for two-character ambiguities */
            StdFreq[A][eBchar] = StdFreq[A][eDchar] + StdFreq[A][eNchar];
            StdFreq[A][eZchar] = StdFreq[A][eEchar] + StdFreq[A][eQchar];
            if (StdAlphsize > eJchar) {
                StdFreq[A][eJchar] = StdFreq[A][eIchar] + StdFreq[A][eLchar];
            }
        }
    }
    /* Add rows to set values for two-character ambiguities */
    memcpy(StdFreq[eBchar], StdFreq[eDchar], StdAlphsize * sizeof(double));
    Nlm_AddVectors(StdFreq[eBchar], StdAlphsize, 1.0, StdFreq[eNchar]);

    memcpy(StdFreq[eZchar], StdFreq[eEchar], StdAlphsize * sizeof(double));
    Nlm_AddVectors(StdFreq[eZchar], StdAlphsize, 1.0, StdFreq[eQchar]);

    if (StdAlphsize > eJchar) {
        memcpy(StdFreq[eJchar],StdFreq[eIchar], StdAlphsize * sizeof(double));
        Nlm_AddVectors(StdFreq[eJchar], StdAlphsize, 1.0, StdFreq[eLchar]);
    }

}

static double
s_CalcAvgScore(double * M, int alphsize, int incM, const double probs[])
{
    int j;                   /* iteration index */
    double score_iX = 0.0;   /* score of character i substituted by X */

    for (j = 0;  j < alphsize;  j++) {
        if (alphaConvert[j] >= 0) {
            /* If the column corresponds to a true amino acid */
            score_iX += M[j * incM] * probs[j];
        }
    }
    return score_iX;
}

static const double kMaximumXscore = -1.0;

static double
s_CalcXScore(double * M, int alphsize, int incM, const double probs[])
{
    return MIN( s_CalcAvgScore(M, alphsize, incM, probs), kMaximumXscore);
}


static void
s_SetXUOScores(double ** M, int alphsize, 
              const double row_probs[], const double col_probs[])
{
    int i;                      /* iteration index */
    double score_XX = 0.0;      /* score of matching an X to an X */
    /* the matrix has alphsize colums (this variable exists just to
       make things easier to read) */
    const int cols = alphsize;  

    for (i = 0;  i < alphsize;  i++) {
        if (alphaConvert[i] >= 0) {
            double avg_iX = s_CalcAvgScore(M[i], alphsize, 1, col_probs);
            M[i][eXchar] = MIN(avg_iX, kMaximumXscore);
            score_XX += avg_iX * row_probs[i];

            M[eXchar][i] = s_CalcXScore(&M[0][i], alphsize, cols, row_probs);
        }
    }
    M[eXchar][eXchar] = MIN(score_XX, kMaximumXscore);

    /* Set X scores for pairwise ambiguity characters */
    M[eBchar][eXchar] = s_CalcXScore(M[eBchar], alphsize, 1, col_probs);
    M[eXchar][eBchar] = s_CalcXScore(&M[0][eBchar], alphsize, cols, row_probs);

    M[eZchar][eXchar] = s_CalcXScore(M[eZchar], alphsize, 1, col_probs);
    M[eXchar][eZchar] = s_CalcXScore(&M[0][eZchar], alphsize, cols, row_probs);
    if( alphsize > eJchar) {
        M[eJchar][eXchar] = s_CalcXScore(M[eJchar], alphsize, 1, col_probs);
        M[eXchar][eJchar] =
            s_CalcXScore(&M[0][eJchar], alphsize, cols, row_probs);
    }
    /* Copy X scores to U and O */
    memcpy(M[eSelenocysteine], M[eXchar], alphsize * sizeof(double));
    for (i = 0;  i < alphsize;  i++) {
        M[i][eSelenocysteine] = M[i][eXchar];
    }
    if (alphsize > eOchar) {
        memcpy(M[eOchar], M[eXchar], alphsize * sizeof(double));
        for (i = 0;  i < alphsize;  i++) {
            M[i][eOchar] = M[i][eXchar];
        }
    }
}

static long Nint(double x)
{
    x += (x >= 0. ? 0.5 : -0.5);
    return (long)x;
}


static void
s_RoundScoreMatrix(int **matrix, int rows, int cols,    
                   double **floatScoreMatrix)
{
    int p, c; /*indices over positions and characters*/

    for (p = 0;  p < rows;  p++) {
        for (c = 0;  c < cols;  c++) {
            if (floatScoreMatrix[p][c] < INT_MIN) {
                matrix[p][c] = INT_MIN;
            } else {
                matrix[p][c] = Nint(floatScoreMatrix[p][c]);
            }
        }
    }
}

static int
s_ScoresStdAlphabet2(Int4 ** Matrix, int Alphsize,
                    double ** target_freq, int ** StartMatrix,
                    const double row_prob[], const double col_prob[],
                    double Lambda, double **Scores)
{
    /* Note: I'm using a rough convention for this routine that uppercase
     * letters refer to quantities in the standard (larger) alphabet
     * and lowercase letters refer to the true amino acid (smaller)
     * alphabet.
     */
    int i;
    /* row and column probabilities in the NCBIstdaa alphabet */
    double RowProb[COMPO_LARGEST_ALPHABET];
    double ColProb[COMPO_LARGEST_ALPHABET];
    /* A double precision score matrix */
    //double ** Scores = Nlm_DenseMatrixNew(Alphsize, Alphsize);
    //if (Scores == NULL) {
        //return -1;
    //}
    s_UnpackLetterProbs(RowProb, Alphsize, row_prob);
    s_SetPairAmbigProbsToSum(RowProb, Alphsize);

    s_UnpackLetterProbs(ColProb, Alphsize, col_prob);
    s_SetPairAmbigProbsToSum(ColProb, Alphsize);

    Blast_TrueAaToStdTargetFreqs(Scores, Alphsize, target_freq);
    Blast_CalcFreqRatios(Scores, Alphsize, RowProb, ColProb);
    Blast_FreqRatioToScore(Scores, Alphsize, Alphsize, Lambda);
    s_SetXUOScores(Scores, Alphsize, RowProb, ColProb);

    s_RoundScoreMatrix(Matrix, Alphsize, Alphsize, Scores);
    //Nlm_DenseMatrixFree(&Scores);

    for (i = 0;  i < Alphsize;  i++) {
        Matrix[i][eStopChar] = StartMatrix[i][eStopChar];
        Matrix[eStopChar][i] = StartMatrix[eStopChar][i];
    }
    return 0;
}

static int
s_ScoresStdAlphabet(Int4 ** Matrix, int Alphsize,
                    double ** target_freq, int ** StartMatrix,
                    const double row_prob[], const double col_prob[],
                    double Lambda)
{
    /* Note: I'm using a rough convention for this routine that uppercase
     * letters refer to quantities in the standard (larger) alphabet
     * and lowercase letters refer to the true amino acid (smaller)
     * alphabet.
     */
    int i;
    /* row and column probabilities in the NCBIstdaa alphabet */
    double RowProb[COMPO_LARGEST_ALPHABET];
    double ColProb[COMPO_LARGEST_ALPHABET];
    /* A double precision score matrix */
    double ** Scores = Nlm_DenseMatrixNew(Alphsize, Alphsize);
    if (Scores == NULL) {
        return -1;
    }
    s_UnpackLetterProbs(RowProb, Alphsize, row_prob);
    s_SetPairAmbigProbsToSum(RowProb, Alphsize);

    s_UnpackLetterProbs(ColProb, Alphsize, col_prob);
    s_SetPairAmbigProbsToSum(ColProb, Alphsize);

    Blast_TrueAaToStdTargetFreqs(Scores, Alphsize, target_freq);
    Blast_CalcFreqRatios(Scores, Alphsize, RowProb, ColProb);
    Blast_FreqRatioToScore(Scores, Alphsize, Alphsize, Lambda);
    s_SetXUOScores(Scores, Alphsize, RowProb, ColProb);

    s_RoundScoreMatrix(Matrix, Alphsize, Alphsize, Scores);
    Nlm_DenseMatrixFree(&Scores);

    for (i = 0;  i < Alphsize;  i++) {
        Matrix[i][eStopChar] = StartMatrix[i][eStopChar];
        Matrix[eStopChar][i] = StartMatrix[eStopChar][i];
    }
    return 0;
}

/** bound on error for Newton's method */
static const double kCompoAdjustErrTolerance = 0.00000001;
/** iteration limit for Newton's method */
static const int kCompoAdjustIterationLimit = 2000;
/** relative entropy of BLOSUM62 */
static const double kFixedReBlosum62 = 0.44;
/** largest permitted value of score(i,X), when scores are computed using
    composition adjustment */

int
Blast_CompositionMatrixAdj(Int4 ** matrix,
                           int alphsize,
                           EMatrixAdjustRule matrix_adjust_rule,
                           int length1,
                           int length2,
                           const double * stdaa_row_probs,
                           const double * stdaa_col_probs,
                           int pseudocounts,
                           double specifiedRE,
                           Blast_CompositionWorkspace * NRrecord,
                           const Blast_MatrixInfo * matrixInfo)
{
    int iteration_count, status;
    double row_probs[COMPO_NUM_TRUE_AA], col_probs[COMPO_NUM_TRUE_AA];
    /* Target RE when optimizing the matrix; zero if the relative
       entropy should not be constrained. */
    double dummy, desired_re = 0.0;

    s_GatherLetterProbs(row_probs, stdaa_row_probs, alphsize);
    s_GatherLetterProbs(col_probs, stdaa_col_probs, alphsize);

    switch (matrix_adjust_rule) {
    case eUnconstrainedRelEntropy:
        desired_re = 0.0;
        break;
    case eRelEntropyOldMatrixNewContext:
        /* Calculate the desired re using the new marginal probs with
           the old matrix */
        status = Blast_EntropyOldFreqNewContext(&desired_re, &dummy,
                                                &iteration_count,
                                                NRrecord->mat_b,
                                                row_probs, col_probs);
        if (status < 0)     /* Error, e.g. memory */
            return status;
        else if (status > 0) /* we could not calculate the desired re */
            desired_re = 0.0; /* so, leave the re unconstrained */

        break;
    case eRelEntropyOldMatrixOldContext:
        desired_re = Blast_TargetFreqEntropy(NRrecord->mat_b);
        break;
    case eUserSpecifiedRelEntropy:
        desired_re = specifiedRE;
        break;
    default:  /* I assert that we can't get here */
        fprintf(stderr, "Unknown flag for setting relative entropy"
                "in composition matrix adjustment");
        exit(1);
    }
    Blast_ApplyPseudocounts(row_probs, length1,
                            NRrecord->first_standard_freq, pseudocounts);
    Blast_ApplyPseudocounts(col_probs, length2,
                            NRrecord->second_standard_freq, pseudocounts);

    status =
        Blast_OptimizeTargetFrequencies(&NRrecord->mat_final[0][0],
                                        COMPO_NUM_TRUE_AA,
                                        &iteration_count,
                                        &NRrecord->mat_b[0][0],
                                        row_probs, col_probs,
                                        (desired_re > 0.0),
                                        desired_re,
                                        kCompoAdjustErrTolerance,
                                        kCompoAdjustIterationLimit);

    if (status != 0)            /* Did not compute the target freqs */
        return status;

    return
        s_ScoresStdAlphabet(matrix, alphsize, NRrecord->mat_final,
                            matrixInfo->startMatrix,
                            row_probs, col_probs,
                            matrixInfo->ungappedLambda);
}

int
Blast_CompositionMatrixAdj2(Int4 ** matrix,
                           int alphsize,
                           EMatrixAdjustRule matrix_adjust_rule,
                           int length1,
                           int length2,
                           const double * stdaa_row_probs,
                           const double * stdaa_col_probs,
                           int pseudocounts,
                           double specifiedRE,
                           Blast_CompositionWorkspace * NRrecord,
                           const Blast_MatrixInfo * matrixInfo,
                           ReNewtonSystem *newton_system,
                           double *z, double *resids_x, double *resids_z, 
                           double *old_scores, double *workspace, double **grads, double **Scores)
{
    int iteration_count, status;
    double row_probs[COMPO_NUM_TRUE_AA], col_probs[COMPO_NUM_TRUE_AA];
    /* Target RE when optimizing the matrix; zero if the relative
       entropy should not be constrained. */
    double dummy, desired_re = 0.0;

    s_GatherLetterProbs(row_probs, stdaa_row_probs, alphsize);
    s_GatherLetterProbs(col_probs, stdaa_col_probs, alphsize);

    switch (matrix_adjust_rule) {
    case eUnconstrainedRelEntropy:
        desired_re = 0.0;
        break;
    case eRelEntropyOldMatrixNewContext:
        /* Calculate the desired re using the new marginal probs with
           the old matrix */
        status = Blast_EntropyOldFreqNewContext(&desired_re, &dummy,
                                                &iteration_count,
                                                NRrecord->mat_b,
                                                row_probs, col_probs);
        if (status < 0)     /* Error, e.g. memory */
            return status;
        else if (status > 0) /* we could not calculate the desired re */
            desired_re = 0.0; /* so, leave the re unconstrained */

        break;
    case eRelEntropyOldMatrixOldContext:
        desired_re = Blast_TargetFreqEntropy(NRrecord->mat_b);
        break;
    case eUserSpecifiedRelEntropy:
        desired_re = specifiedRE;
        break;
    default:  /* I assert that we can't get here */
        fprintf(stderr, "Unknown flag for setting relative entropy"
                "in composition matrix adjustment");
        exit(1);
    }
    Blast_ApplyPseudocounts(row_probs, length1,
                            NRrecord->first_standard_freq, pseudocounts);
    Blast_ApplyPseudocounts(col_probs, length2,
                            NRrecord->second_standard_freq, pseudocounts);

    status =
        Blast_OptimizeTargetFrequencies2(&NRrecord->mat_final[0][0],
                                        COMPO_NUM_TRUE_AA,
                                        &iteration_count,
                                        &NRrecord->mat_b[0][0],
                                        row_probs, col_probs,
                                        (desired_re > 0.0),
                                        desired_re,
                                        kCompoAdjustErrTolerance,
                                        kCompoAdjustIterationLimit,
                                        newton_system,
                                        z, resids_x, resids_z, old_scores, workspace, grads
                                        );

    if (status != 0)            /* Did not compute the target freqs */
        return status;

    return
        s_ScoresStdAlphabet2(matrix, alphsize, NRrecord->mat_final,
                            matrixInfo->startMatrix,
                            row_probs, col_probs,
                            matrixInfo->ungappedLambda, Scores);
}

static int trueCharPositions[COMPO_NUM_TRUE_AA] =
{1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22};


static void s_GetScoreRange(int * obs_min, int * obs_max,
                            int ** matrix, int rows)
{
    int aa;                    /* index of an amino-acid in the 20
                                  letter alphabet */
    int irow, jcol;            /* matrix row and column indices */
    int minScore, maxScore;    /* largest and smallest observed scores */

    minScore = maxScore = 0;
    for (irow = 0;  irow < rows;  irow++) {
        for (aa = 0;  aa < COMPO_NUM_TRUE_AA;  aa++) {
            jcol = trueCharPositions[aa];
            if (matrix[irow][jcol] < minScore &&
                matrix[irow][jcol] > COMPO_SCORE_MIN)
                minScore = matrix[irow][jcol];
            if (matrix[irow][jcol] > maxScore)
                maxScore = matrix[irow][jcol];
        }
    }
    *obs_min = minScore;
    *obs_max = maxScore;
}


static int
s_GetMatrixScoreProbs(double **scoreProb, int * obs_min, int * obs_max,
                      int **matrix, int alphsize, 
                      const double *subjectProbArray,
                      const double *queryProbArray)
{
    int aa;          /* index of an amino-acid in the 20 letter
                        alphabet */
    int irow, jcol;  /* matrix row and column indices */
    double * sprob;  /* a pointer to the element of the score
                        probabilities array that represents the
                        probability of the score 0*/
    int minScore;    /* smallest score in matrix; the same value as
                        (*obs_min). */
    int range;       /* the range of scores in the matrix */

    s_GetScoreRange(obs_min, obs_max, matrix, alphsize);
    minScore = *obs_min;
    range = *obs_max - *obs_min + 1;
    *scoreProb = calloc(range, sizeof(double));
    if (*scoreProb == NULL) {
        return -1;
    }
    sprob = &((*scoreProb)[-(*obs_min)]); /*center around 0*/
    for (irow = 0;  irow < alphsize;  irow++) {
        for (aa = 0;  aa < COMPO_NUM_TRUE_AA;  aa++) {
            jcol = trueCharPositions[aa];
            if (matrix[irow][jcol] >= minScore) {
                sprob[matrix[irow][jcol]] +=
                    (queryProbArray[irow] * subjectProbArray[jcol]);
            }
        }
    }
    return 0;
}

#define LambdaRatioLowerBound 0.5

static int
s_ScaleSquareMatrix(int **matrix, int alphsize,
                    double ** freq_ratios, int ** start_matrix,
                    const double row_prob[], const double col_prob[],
                    double Lambda)
{
    double ** scores;     /* a double precision matrix of scores */
    int i;                /* iteration index */

    scores = Nlm_DenseMatrixNew(alphsize, alphsize);
    if (scores == 0) return -1;

    for (i = 0;  i < alphsize;  i++) {
        memcpy(scores[i], freq_ratios[i], alphsize * sizeof(double));
    }
    Blast_FreqRatioToScore(scores, alphsize, alphsize, Lambda);
    s_SetXUOScores(scores, alphsize, row_prob, col_prob);
    s_RoundScoreMatrix(matrix, alphsize, alphsize, scores);
    for (i = 0;  i < alphsize;  i++) {
        matrix[i][eStopChar] = start_matrix[i][eStopChar];
        matrix[eStopChar][i] = start_matrix[eStopChar][i];
    }
    Nlm_DenseMatrixFree(&scores);

    return 0;
}


int
Blast_CompositionBasedStats(int ** matrix, double * LambdaRatio,
                            const Blast_MatrixInfo * ss,
                            const double queryProb[], const double resProb[],
                            //double (*calc_lambda)(double*,int,int,double),
                            int pValueAdjustment)
{
    double correctUngappedLambda; /* new value of ungapped lambda */
    int obs_min, obs_max;         /* smallest and largest score in the
                                     unscaled matrix */
    double *scoreArray;           /* an array of score probabilities */
    int out_of_memory;            /* status flag to indicate out of memory */

    //if (ss->positionBased) {
        //out_of_memory =
            //s_GetPssmScoreProbs(&scoreArray, &obs_min, &obs_max,
                                //ss->startMatrix, ss->rows, resProb);
    //} else 
    {
        out_of_memory =
            s_GetMatrixScoreProbs(&scoreArray, &obs_min, &obs_max,
                                  ss->startMatrix, ss->cols,
                                  resProb, queryProb);
    }

    //fprintf(stderr, "obs_min: %d obs_max: %d\n", obs_min, obs_max);
    if (out_of_memory)
        return -1;
    correctUngappedLambda =
        s_CalcLambda(scoreArray, obs_min, obs_max, ss->ungappedLambda);


    /* calc_lambda will return -1 in the case where the
     * expected score is >=0; however, because of the MAX statement 3
     * lines below, LambdaRatio should always be > 0; the succeeding
     * test is retained as a vestige, in case one wishes to remove the
     * MAX statement and allow LambdaRatio to take on the error value
     * -1 */
    *LambdaRatio = correctUngappedLambda / ss->ungappedLambda;
    //fprintf(stderr, "correctUngappedLambda: %f LambdaRatio: %f\n", correctUngappedLambda, *LambdaRatio);
    if (0 == pValueAdjustment)
      *LambdaRatio = MIN(1, *LambdaRatio);
    *LambdaRatio = MAX(*LambdaRatio, LambdaRatioLowerBound);

    if (*LambdaRatio > 0) {
        double scaledLambda = ss->ungappedLambda/(*LambdaRatio);
        //if (ss->positionBased) {
            //s_ScalePSSM(matrix, ss->rows, ss->cols, ss->startFreqRatios,
                        //ss->startMatrix, resProb, scaledLambda);
        //} else 
        {
            s_ScaleSquareMatrix(matrix, ss->cols,
                                ss->startFreqRatios, ss->startMatrix,
                                queryProb, resProb, scaledLambda);
        }
    }
    free(scoreArray);

    return 0;
}

#define BLAST_SCORE_MIN INT2_MIN
#define BLAST_SCORE_MAX INT2_MAX

#define BLAST_SCORE_RANGE_MAX (BLAST_SCORE_MAX - BLAST_SCORE_MIN)

static Int2
BlastScoreChk(Int4 lo, Int4 hi)
{
   if (lo >= 0 || hi <= 0 ||
         lo < BLAST_SCORE_MIN || hi > BLAST_SCORE_MAX)
      return 1;

   if (hi - lo > BLAST_SCORE_RANGE_MAX)
      return 1;

   return 0;
}

//double 
//NlmKarlinLambdaNR(double* probs, Int4 d, Int4 low, Int4 high, double lambda0,
                  //double tolx, Int4 itmax, Int4 maxNewton, Int4 * itn ) 
//{
  //Int4 k;
  //double x0, x, a = 0, b = 1;
  //double f = 4;  /* Larger than any possible value of the poly in [0,1] */
  //Int4 isNewton = 0; /* we haven't yet taken a Newton step. */

  ////assert( d > 0 );

   //x0 = exp( -lambda0 );
  //x = ( 0 < x0 && x0 < 1 ) ? x0 : .5;
  
  //for( k = 0; k < itmax; k++ ) { /* all iteration indices k */
    //Int4 i;
    //double g, fold = f;
    //Int4 wasNewton = isNewton; /* If true, then the previous step was a */
                              ///* Newton step */
    //isNewton  = 0;            /* Assume that this step is not */
    
    ///* Horner's rule for evaluating a polynomial and its derivative */
    //g = 0;
    //f = probs[low];
    //for( i = low + d; i < 0; i += d ) {
      //g = x * g + f;
      //f = f * x + probs[i];
    //}
    //g = x * g + f;
    //f = f * x + probs[0] - 1;
    //for( i = d; i <= high; i += d ) {
      //g = x * g + f;
      //f = f * x + probs[i];
    //}
    ///* End Horner's rule */

    //if( f > 0 ) {
      //a = x; /* move the left endpoint */
    //} else if( f < 0 ) { 
      //b = x; /* move the right endpoint */
    //} else { /* f == 0 */
      //break; /* x is an exact solution */
    //}
    //if( b - a < 2 * a * ( 1 - b ) * tolx ) {
      ///* The midpoint of the interval converged */
      //x = (a + b) / 2; break;
    //}

    //if( k >= maxNewton ||
        ///* If convergence of Newton's method appears to be failing; or */
            //( wasNewton && fabs( f ) > .9 * fabs(fold) ) ||  
        ///* if the previous iteration was a Newton step but didn't decrease 
         //* f sufficiently; or */
        //g >= 0 
        ///* if a Newton step will move us away from the desired solution */
        //) { /* then */
      ///* bisect */
      //x = (a + b)/2;
    //} else {
      ///* try a Newton step */
      //double p = - f/g;
      //double y = x + p;
      //if( y <= a || y >= b ) { /* The proposed iterate is not in (a,b) */
        //x = (a + b)/2;
      //} else { /* The proposed iterate is in (a,b). Accept it. */
        //isNewton = 1;
        //x = y;
        //if( fabs( p ) < tolx * x * (1-x) ) break; /* Converged */
      //} /* else the proposed iterate is in (a,b) */
    //} /* else try a Newton step. */ 
  //} /* end for all iteration indices k */
   //*itn = k; 
  //return -log(x)/d;
//}

#define BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT                                   \
    (1.e-5) /**< LAMBDA_ACCURACY_DEFAULT == accuracy to which Lambda should be   \
              calc'd */

#define BLAST_KARLIN_LAMBDA_ITER_DEFAULT                                       \
    17 /**< LAMBDA_ITER_DEFAULT == no. of iterations in LambdaBis =              \
         ln(accuracy)/ln(2)*/



double
s_CalcLambda(double probs[], int min_score, int max_score, double lambda0)
{
   
    int i;                 /* loop index */      
    int score_range;       /* range of possible scores */
    double avg;            /* expected score of aligning two characters */
    Blast_ScoreFreq freq;  /* score frequency data */

    score_range = max_score - min_score + 1;
    avg = 0.0;
    for (i = 0;  i < score_range;  i++) {
        avg += (min_score + i) * probs[i];
        //fprintf(stderr, "%f ", probs[i]);
    }
    //fprintf(stderr, "\n");
    freq.score_min = min_score;
    freq.score_max = max_score;
    freq.obs_min = min_score;
    freq.obs_max = max_score;
    freq.sprob0 = probs;
    freq.sprob = &probs[-min_score];
    freq.score_avg = avg;

    return Blast_KarlinLambdaNR(&freq, lambda0);
}

int
Blast_AdjustScores2(Int4 ** matrix,
                   const Blast_AminoAcidComposition * query_composition,
                   int queryLength,
                   const Blast_AminoAcidComposition * subject_composition,
                   int subjectLength,
                   const Blast_MatrixInfo * matrixInfo,
                   ECompoAdjustModes composition_adjust_mode,
                   int RE_pseudocounts,
                   Blast_CompositionWorkspace *NRrecord,
                   EMatrixAdjustRule *matrix_adjust_rule,
                   double *pvalueForThisPair,
                   int compositionTestIndex,
                   double *ratioToPassBack,
                   ReNewtonSystem *newton_system,
                   double *z, double *resids_x, double *resids_z, 
                   double *old_scores, double *workspace, 
                   double **grads, double **Scores)
{
    const int alphsize = matrixInfo->cols;

    double lambdaForPair;     /*lambda for this pair of compositions*/
    int iter_count; /*used as argument to Blast_CalcLambdaFullPrecision*/

    /* The next two arrays are letter probabilities of query and
     * match in 20 letter ARND... alphabet. */
    double permutedQueryProbs[COMPO_NUM_TRUE_AA];
    double permutedMatchProbs[COMPO_NUM_TRUE_AA];

    if (query_composition->numTrueAminoAcids == 0 ||
        subject_composition->numTrueAminoAcids == 0) {
        /* Either the query or subject contains only ambiguity
           characters, most likely because the entire subject has been
           SEGed.  Compositional adjustment is meaningless. */
        return 1;
    }

    if ((compositionTestIndex > 0) ||
            ((!(matrixInfo->positionBased)) &&
             (composition_adjust_mode != eCompositionBasedStats))) {
        s_GatherLetterProbs(permutedQueryProbs,
                query_composition->prob, alphsize);
        s_GatherLetterProbs(permutedMatchProbs,
                subject_composition->prob, alphsize);
    }

    *matrix_adjust_rule =
        Blast_ChooseMatrixAdjustRule(queryLength, subjectLength,
                permutedQueryProbs,
                permutedMatchProbs,
                matrixInfo->matrixName,
                composition_adjust_mode);
   
    if (eCompoScaleOldMatrix != *matrix_adjust_rule) {
        /* Try matrix optimization, if it fails to converge, we
           fall back to traditional scaling below */
        int status =
            Blast_CompositionMatrixAdj2(matrix,
                                       alphsize,
                                       *matrix_adjust_rule,
                                       query_composition->
                                       numTrueAminoAcids,
                                       subject_composition->
                                       numTrueAminoAcids,
                                       query_composition->prob,
                                       subject_composition->prob,
                                       RE_pseudocounts,
                                       kFixedReBlosum62,
                                       NRrecord,
                                       matrixInfo,
                                       newton_system,
                                       z, resids_x, resids_z, old_scores, workspace, grads, Scores
                                       );

        


        *ratioToPassBack = 1.0;    /* meaningless for this mode */
        if (status <= 0)
        {

            //fprintf(stderr, "Blast_CompositionMatrixAdj <= 0\n");
            return status;      /* Success (=0) or fatal error (<0)*/
        }
            
    } /* End try matrix optimization */

    *matrix_adjust_rule = eCompoScaleOldMatrix; 
    return Blast_CompositionBasedStats(matrix, ratioToPassBack, matrixInfo,
                                       query_composition->prob,
                                       subject_composition->prob,
                                       //calc_lambda,
                                       (compositionTestIndex > 0));
}

int
Blast_AdjustScores(Int4 ** matrix,
                   const Blast_AminoAcidComposition * query_composition,
                   int queryLength,
                   const Blast_AminoAcidComposition * subject_composition,
                   int subjectLength,
                   const Blast_MatrixInfo * matrixInfo,
                   ECompoAdjustModes composition_adjust_mode,
                   int RE_pseudocounts,
                   Blast_CompositionWorkspace *NRrecord,
                   EMatrixAdjustRule *matrix_adjust_rule,
                   double *pvalueForThisPair,
                   int compositionTestIndex,
                   double *ratioToPassBack)
{
    const int alphsize = matrixInfo->cols;

    double lambdaForPair;     /*lambda for this pair of compositions*/
    int iter_count; /*used as argument to Blast_CalcLambdaFullPrecision*/

    /* The next two arrays are letter probabilities of query and
     * match in 20 letter ARND... alphabet. */
    double permutedQueryProbs[COMPO_NUM_TRUE_AA];
    double permutedMatchProbs[COMPO_NUM_TRUE_AA];

    if (query_composition->numTrueAminoAcids == 0 ||
        subject_composition->numTrueAminoAcids == 0) {
        /* Either the query or subject contains only ambiguity
           characters, most likely because the entire subject has been
           SEGed.  Compositional adjustment is meaningless. */
        return 1;
    }

    if ((compositionTestIndex > 0) ||
            ((!(matrixInfo->positionBased)) &&
             (composition_adjust_mode != eCompositionBasedStats))) {
        s_GatherLetterProbs(permutedQueryProbs,
                query_composition->prob, alphsize);
        s_GatherLetterProbs(permutedMatchProbs,
                subject_composition->prob, alphsize);
    }

    /* else call Yi-Kuo's code to choose mode for matrix adjustment. */
    *matrix_adjust_rule =
        Blast_ChooseMatrixAdjustRule(queryLength, subjectLength,
                permutedQueryProbs,
                permutedMatchProbs,
                matrixInfo->matrixName,
                composition_adjust_mode);
   //fprintf(stderr, "matrix_adjust_rule: %d\n", *matrix_adjust_rule);

   
    if (eCompoScaleOldMatrix != *matrix_adjust_rule) {
        /* Try matrix optimization, if it fails to converge, we
           fall back to traditional scaling below */
        int status =
            Blast_CompositionMatrixAdj(matrix,
                                       alphsize,
                                       *matrix_adjust_rule,
                                       query_composition->
                                       numTrueAminoAcids,
                                       subject_composition->
                                       numTrueAminoAcids,
                                       query_composition->prob,
                                       subject_composition->prob,
                                       RE_pseudocounts,
                                       kFixedReBlosum62,
                                       NRrecord,
                                       matrixInfo);

        


        *ratioToPassBack = 1.0;    /* meaningless for this mode */
        if (status <= 0)
        {

            //fprintf(stderr, "Blast_CompositionMatrixAdj <= 0\n");
            return status;      /* Success (=0) or fatal error (<0)*/
        }
            
    } /* End try matrix optimization */

    *matrix_adjust_rule = eCompoScaleOldMatrix; 
    return Blast_CompositionBasedStats(matrix, ratioToPassBack, matrixInfo,
                                       query_composition->prob,
                                       subject_composition->prob,
                                       (compositionTestIndex > 0));
}


static const int kCompositionMargin = 20;

void
Blast_GetCompositionRange(int * pleft, int * pright,
                          const Uint1 * subject_data, int length,
                          int start, int finish)
{
    int i;                /* iteration index */
    int left, right;

    left = start;
    /* Search leftward for a StopChar */
    for (i = left;  i > 0;  i--) {
        if (subject_data[i - 1] == eStopChar) {
            /* We have found a StopChar. Unless the StopChar is
             * too close to the start of the subject region of the
             * HSP, */
            if (i + kCompositionMargin < left) {
                /* reset the left endpoint. */
                left = i + kCompositionMargin;
            }
            break;
        }
    }
    if (i == 0) {
        /* No stop codon was found to the left. */
        left = 0;
    }
    right = finish;
    /* Search rightward for a StopChar */
    for (i = right;  i < length;  i++) {
        if (subject_data[i] == eStopChar) {
            /* We have found a StopChar. Unless the StopChar is
             * too close to the end of the subject region of the
             * HSP, */
            if (i - kCompositionMargin > right) {
                /* reset the right endpoint */
                right = i - kCompositionMargin;
            }
            break;
        }
    }
    if (i == length) {
        /* No stop codon was found to the right. */
        right = length;
    }
    *pleft = left; *pright = right;
}


static const char kNCBIstdaa[] = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";

void
Blast_ReadAaComposition_fsa(Blast_AminoAcidComposition * composition,
        int alphsize,
        const Uint1 * sequence, int length)
{
    int i; /* iteration index */

    /* fields of composition as local variables */
    int numTrueAminoAcids = 0;
    double * prob = composition->prob;

    for (i = 0;  i < alphsize;  i++) {
        prob[i] = 0.0;
    }
    for (i = 0;  i < length;  i++) {
        //fprintf(stderr, "%c", encoding_getLetter(sequence[i]));
        if (alphaConvert[to_ncbi[sequence[i]]] >= 0) {
            prob[to_ncbi[sequence[i]]]++;
            numTrueAminoAcids++;
        }
    }
    //fprintf(stderr, "\n");
    composition->numTrueAminoAcids = numTrueAminoAcids;
    if (numTrueAminoAcids > 0) {
        for (i = 0;  i < alphsize;  i++) {
            prob[i] /= numTrueAminoAcids;
        }
    }
}

void
Blast_ReadAaComposition(Blast_AminoAcidComposition * composition,
        int alphsize,
        const Uint1 * sequence, int length)
{
    int i; /* iteration index */

    /* fields of composition as local variables */
    int numTrueAminoAcids = 0;
    double * prob = composition->prob;

    for (i = 0;  i < alphsize;  i++) {
        prob[i] = 0.0;
    }
    for (i = 0;  i < length;  i++) {
        //fprintf(stderr, "%c", kNCBIstdaa[sequence[i]]);
        if (alphaConvert[sequence[i]] >= 0) {
            prob[sequence[i]]++;
            numTrueAminoAcids++;
        }
    }
    //fprintf(stderr, "\n");
    composition->numTrueAminoAcids = numTrueAminoAcids;
    if (numTrueAminoAcids > 0) {
        for (i = 0;  i < alphsize;  i++) {
            prob[i] /= numTrueAminoAcids;
        }
    }
}

void
s_GetComposition(Blast_AminoAcidComposition * composition,
        int alphsize,
        Uint1 *data,
        int length
        )
{
    int left, right;
    left = 0;
    right = length;
    Blast_ReadAaComposition(composition, alphsize, &data[left], right-left);
}

void
Blast_CompositionWorkspaceFree(Blast_CompositionWorkspace ** pNRrecord)
{
    Blast_CompositionWorkspace * NRrecord = *pNRrecord;

    if (NRrecord != NULL) {
        free(NRrecord->first_standard_freq);
        free(NRrecord->second_standard_freq);

        Nlm_DenseMatrixFree(&NRrecord->mat_final);
        Nlm_DenseMatrixFree(&NRrecord->mat_b);

        free(NRrecord);
    }
    pNRrecord = NULL;
}


Blast_CompositionWorkspace * Blast_CompositionWorkspaceNew(void)
{
    Blast_CompositionWorkspace * NRrecord;        /* record to allocate
                                                    and return */
    int i;                     /* loop index */

    NRrecord = (Blast_CompositionWorkspace *)
        malloc(sizeof(Blast_CompositionWorkspace));
    if (NRrecord == NULL) goto error_return;

    NRrecord->first_standard_freq      = NULL;
    NRrecord->second_standard_freq     = NULL;
    NRrecord->mat_final                = NULL;
    NRrecord->mat_b                    = NULL;

    NRrecord->first_standard_freq =
        (double *) malloc(COMPO_NUM_TRUE_AA * sizeof(double));
    if (NRrecord->first_standard_freq == NULL) goto error_return;

    NRrecord->second_standard_freq =
        (double *) malloc(COMPO_NUM_TRUE_AA * sizeof(double));
    if (NRrecord->second_standard_freq == NULL) goto error_return;

    NRrecord->mat_final   = Nlm_DenseMatrixNew(COMPO_NUM_TRUE_AA,
                                               COMPO_NUM_TRUE_AA);
    if (NRrecord->mat_final == NULL) goto error_return;

    NRrecord->mat_b       = Nlm_DenseMatrixNew(COMPO_NUM_TRUE_AA,
                                               COMPO_NUM_TRUE_AA);
    if (NRrecord->mat_b == NULL) goto error_return;

    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        NRrecord->first_standard_freq[i] =
            NRrecord->second_standard_freq[i] = 0.0;
    }

    goto normal_return;
error_return:
    Blast_CompositionWorkspaceFree(&NRrecord);
normal_return:
    return NRrecord;
}




/* Documented in matrix_frequency_data.h. */
int
Blast_GetJointProbsForMatrix(double ** probs, double row_sums[],
                             double col_sums[], const char *matrix_name)
{
    int i, j;
    Compo_FrequencyData * data = s_LocateFrequencyData(matrix_name);
    if (NULL == data) {
        fprintf(stderr, "matrix %s is not supported "
                "for RE based adjustment\n", matrix_name);
        return -1;
    }
    for (j = 0;  j < COMPO_NUM_TRUE_AA;  j++) {
        col_sums[j] = 0.0;
    }
    for (i = 0;  i < COMPO_NUM_TRUE_AA;  i++) {
        row_sums[i] = 0.0;
        for (j = 0;  j < COMPO_NUM_TRUE_AA;  j++) {
            probs[i][j] = data->joint_probs[i][j];
            row_sums[i] += probs[i][j];
            col_sums[j] += probs[i][j];
        }
    }
    return 0;
}



    int
Blast_CompositionWorkspaceInit(Blast_CompositionWorkspace * NRrecord,
        const char *matrixName)
{
    if (0 == Blast_GetJointProbsForMatrix(NRrecord->mat_b,
                NRrecord->first_standard_freq,
                NRrecord->second_standard_freq,
                matrixName)) {
        return 0;
    } else {
        fprintf(stderr,
                "Matrix %s not currently supported for RE based adjustment\n",
                matrixName);
        return -1;
    }
}

#define DBL_MAX 1.79769e+308

void alignments_get_eValue(struct ungappedExtension **ungappedExtensions, int4 numExts, Int4 queryLength, Int4 subjectLength, double *best_eValue)
{
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
    gbp.db_length = total_numberOfLetters; 
    gbp.filled = 1;

    Blast_HSPListGetEvalues(queryLength, subjectLength, &kbp, &gbp, ungappedExtensions, numExts, TRUE, FALSE, 0.0, 1.0, best_eValue);

}

void finalAlignments_get_eValue(struct ungappedExtension **ungappedExtensions, int4 numExts, Int4 queryLength, Int4 subjectLength, double *best_eValue)
{
    Blast_KarlinBlk kbp;
    Blast_GumbelBlk gbp;

    kbp.Lambda = 0.0083437500000000005; 
    kbp.K = 0.041000000000000002;
    kbp.logK = -3.1941832122778293;
    kbp.H = 0.14000000000000001;
    kbp.paramC = 0; 

    gbp.Lambda = 0.26700000000000002;
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
    gbp.db_length = total_numberOfLetters; 
    gbp.filled = 1;

    Blast_HSPListGetEvalues(queryLength, subjectLength, &kbp, &gbp, ungappedExtensions, numExts, TRUE, FALSE, 0.0, 1.0, best_eValue);

}

Int2 Blast_HSPListReapByEvalue(struct ungappedExtension ** hsp_array, Int4 hspcnt, Int4 *hspcnt_new, Int4 *bestScore, double expect_value)
{
    struct ungappedExtension* hsp;
    *bestScore = INT2_MIN;
    //BlastHSP** hsp_array;
    //Int4 hsp_cnt = 0;
   Int4 index;
   double cutoff;
   
   cutoff = expect_value;

   *hspcnt_new = 0;

   for (index = 0; index < hspcnt; index++) {
      hsp = hsp_array[index];
      *bestScore = MAX(*bestScore, hsp->nominalScore);

      ASSERT(hsp != NULL);
      
      if (hsp->eValue > cutoff) {
          hsp->status = ungappedExtension_DELETED;
      }
      else {
          if (index > *hspcnt_new)
              hsp_array[*hspcnt_new] = hsp_array[index];
          (*hspcnt_new)++;
      }
   }
      
   return 0;
}


int
s_HitlistEvaluateAndPurge(int * pbestScore, double *pbestEvalue,
                          struct ungappedExtension ** hsp_list,
                          int hspcnt,
                          int query_length,
                          int subject_length)
{
    int status = 0;
    *pbestEvalue = DBL_MAX;
    *pbestScore  = 0;

    finalAlignments_get_eValue(hsp_list, 
            hspcnt, query_length, 
            subject_length, pbestEvalue);
    int hspcnt_new = 0;
    double expect_value = 10;
    Blast_HSPListReapByEvalue(hsp_list, 
            hspcnt, &hspcnt_new, 
            pbestScore, expect_value);

    if(!hspcnt_new)
    {
        *pbestEvalue = 0;
        *pbestScore = 0;
    }

    return hspcnt_new;
}

