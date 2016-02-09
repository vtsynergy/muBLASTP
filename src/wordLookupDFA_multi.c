/*
* muBLASTP - A database indexed protein sequence search tool 
* Version 0.1 (beta)
*
* (c) 2010-2015 Virginia Tech
* This version of dbBLASTP is licensed for non-commercial use only,
*  as specified in LICENSE. For all other use contact vtiplicensing@vtip.org
* 
* Authors: Jing Zhang 
*
*/

#include "blast.h"

uint2 *wordLookupDFA_additionalQueryPositions_multi[BATCH_SIZE] = { NULL };
int4 wordLookupDFA_numAdditionalQueryPositions_multi[BATCH_SIZE];
struct group *wordLookupDFA_groups_multi[BATCH_SIZE] = { NULL };
int4 wordLookupDFA_numBlocks_multi[BATCH_SIZE] = { 0 };

unsigned char *wordLookupDFA_multi[BATCH_SIZE] = { NULL };
unsigned char **wordLookupDFA_blockAddresses_multi[BATCH_SIZE] = { NULL };

struct initialWord {
  int2 numQueryPositions;
  int2 allocQueryPositions;
  int2 *queryPositions;
};

struct aaFrequencyGroup {
  struct memSingleBlock *groups;
  float frequency;
};

void wordLookupDFA_assign(int queryNum) {
  wordLookupDFA_additionalQueryPositions_multi[queryNum] =
      wordLookupDFA_additionalQueryPositions;
  wordLookupDFA_numAdditionalQueryPositions_multi[queryNum] =
      wordLookupDFA_numAdditionalQueryPositions;
  wordLookupDFA_groups_multi[queryNum] = wordLookupDFA_groups;
  wordLookupDFA_numBlocks_multi[queryNum] = wordLookupDFA_numBlocks;
  wordLookupDFA_multi[queryNum] = wordLookupDFA;
  wordLookupDFA_blockAddresses_multi[queryNum] = wordLookupDFA_blockAddresses;
}

void wordLookupDFA_free_multi(int queryNum) {
  free(wordLookupDFA_groups_multi[queryNum]);
  free(wordLookupDFA_multi[queryNum]);
  free(wordLookupDFA_blockAddresses_multi[queryNum]);
}

// Print the contents of the word lookup table
void wordLookupDFA_print_multi(int queryNum) {
  unsigned char *currentBlock;
  unsigned char codes[wordLookupDFA_wordLength + 1], code;
  unsigned char word;
  int4 groupCodeword;
  uint2 *queryPositions;
  int4 totalQueryPositions = 0, totalEmptySlots = 0;
  int4 count;
  int4 numWords, numGroups;
  int4 numCodes, wordLength;
  uint4 codeCount;

  numCodes = wordLookupDFA_numCodes;
  wordLength = wordLookupDFA_wordLength;

  // Determine number of codewords and number of groups
  numGroups = ceil(pow(numCodes, wordLength - 1));
  numWords = ceil(pow(numCodes, wordLength));

  // Initialize word to first position
  codeCount = 0;
  while (codeCount <= wordLength - 1) {
    codes[codeCount] = 0;
    codeCount++;
  }

  // Iterate - For each possible group
  while (!codes[wordLength - 1]) {
    // Construct the codeword for array of codes
    groupCodeword = wordLookupDFA_getCodeword(codes, wordLength - 1);

    // Get current block
    currentBlock = wordLookupDFA_blockAddresses_multi[queryNum][groupCodeword];

    // For each word in the group
    code = 0;
    while (code < numCodes) {
      // Get word in the block
      word = currentBlock[code];

      if (encoding_alphabetType == encoding_protein || word != 0) {
        //  Print the word
        codeCount = 0;
        while (codeCount < wordLength - 1) {
          printf("%c", encoding_getLetter(codes[codeCount]));
          codeCount++;
        }
        printf("%c", encoding_getLetter(code));

        printf(" QueryPositions ");
        fflush(stdout);

        printf("[%d]: ", word);

        if (word == 0) {
          // No query positions
          totalEmptySlots++;
        } else {
          // At least one position at an external address
          queryPositions = ((uint2 *)currentBlock) - word;

          // If the zero flag is stored at the first query position
          if (!*queryPositions) {
            // Go to an outside address for additional positions
            queryPositions =
                wordLookupDFA_additionalQueryPositions_multi[queryNum] +
                *(queryPositions + 1);
          }

          count = 0;
          while (queryPositions[count] != 0) {
            printf("%d ", queryPositions[count] - 1);
            totalQueryPositions++;
            fflush(stdout);
            count++;
          }
        }

        printf("\n");
      }

      code++;
    }

    // Move to next word
    codes[0]++;
    codeCount = 0;
    while (codeCount < wordLength - 1) {
      if (codes[codeCount] >= numCodes) {
        codes[codeCount] = 0;
        codes[codeCount + 1]++;
      } else {
        break;
      }
      codeCount++;
    }
  }

  printf("Total query positions=%d\n", totalQueryPositions);
  printf("Empty slots=%d/%d (%d%%)\n", totalEmptySlots, numWords,
         totalEmptySlots * 100 / numWords);
  printf("Number of blocks/groups=%d/%d\n",
         wordLookupDFA_numBlocks_multi[queryNum], numGroups);
}
