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

// Get an unused entry from the block
void *memBlocks_newEntry_multi(struct memBlocks *memBlocks) {
  void *newEntry;

  // Check if we need to create a new block of memory
  if (memBlocks->numEntries[memBlocks->numBlocks - 1] >=
      memBlocks->blockSizes) {
    // Declare memory for the new block
    memBlocks->lastBlock =
        (void *)global_malloc(memBlocks->entrySize * memBlocks->blockSizes);

    // Check if we need more memory for block pointers
    if (memBlocks->numBlocks >= memBlocks->maxNumBlocks) {
      // Allocate more
      memBlocks->maxNumBlocks *= 2;
      memBlocks->blocks = (void **)global_realloc(
          memBlocks->blocks, sizeof(void *) * memBlocks->maxNumBlocks);
      memBlocks->numEntries = (int4 *)global_realloc(
          memBlocks->numEntries, sizeof(int4) * memBlocks->maxNumBlocks);
    }

    // Store the address of this new block
    memBlocks->blocks[memBlocks->numBlocks] = memBlocks->lastBlock;

    // Reset number of entries in this block
    memBlocks->numEntries[memBlocks->numBlocks] = 0;
    memBlocks->numBlocks++;
  }

  // Use the next available slot in the latest block
  newEntry =
      ((char *)(memBlocks->lastBlock)) +
      memBlocks->numEntries[memBlocks->numBlocks - 1] * memBlocks->entrySize;

  memBlocks->numEntries[memBlocks->numBlocks - 1]++;
  memBlocks->numTotalEntries++;
  return newEntry;
}

// Set the current position to the specific entry
void memBlocks_setCurrent(struct memBlocks *memBlocks, int entryNum) {
  memBlocks->currentBlock = entryNum / memBlocks->blockSizes;
  memBlocks->currentEntry = entryNum % memBlocks->blockSizes;
}
