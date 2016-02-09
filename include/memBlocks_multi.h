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



#ifndef _memBlocks_multi__
#define _memBlocks_multi__

void* memBlocks_newEntry_multi(struct memBlocks* memBlocks);

void memBlocks_setCurrent(struct memBlocks* memBlocks, int entryNum);
#endif
