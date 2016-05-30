
INC=include
SRC=src
OBJ=obj

NCBI_CORE=ncbi_core/

CC=gcc
CCPLUS=g++
CFLAGS=-fopenmp -g -I$(INC) -I$(NCBI_CORE) -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE #-DNEIGHBOR_INDEX
LDFLAGS=-g -lm -lpthread -lstdc++

OBJFILES=$(OBJ)/utils.o $(OBJ)/merge.o $(OBJ)/queryLookup.o $(OBJ)/dbLookup.o\
$(OBJ)/search_dbIdx.o $(OBJ)/alignments_multi.o $(OBJ)/alignments_dbIdx.o\
$(OBJ)/descriptions.o $(OBJ)/encoding.o  $(OBJ)/global.o $(OBJ)/karlin.o\
$(OBJ)/memBlocks.o $(OBJ)/memSingleBlock.o $(OBJ)/parameters.o $(OBJ)/print.o\
$(OBJ)/PSSMatrix.o $(OBJ)/qPosList.o $(OBJ)/readFasta.o $(OBJ)/readFile.o\
$(OBJ)/scoreMatrix.o  $(OBJ)/vbyte.o  $(OBJ)/statistics.o\
$(OBJ)/ungappedExtension_multi.o  $(OBJ)/ungappedExtension.o\
$(OBJ)/wordLookupDFA_multi.o $(OBJ)/wordLookupDFA.o $(OBJ)/writeFile.o\
$(OBJ)/constants.o $(OBJ)/writedb.o $(OBJ)/readdb.o $(OBJ)/wildcards.o\
$(OBJ)/dust.o $(OBJ)/seg.o $(OBJ)/ncbi_gapalign.o $(OBJ)/ncbi_stat.o\
$(OBJ)/ncbi_tree.o $(OBJ)/ncbi_traceback.o $(OBJ)/ncbi_matrix.o\
$(OBJ)/ncbi_composition_adjustment.o $(OBJ)/ncbi_filter.o\
$(OBJ)/alignments_ncbi.o $(OBJ)/dbLookup_compress.o $(OBJ)/blast_kappa.o 


HEADERFILES=$(INC)/utils.h $(INC)/merge.h $(INC)/queryLookup.h $(INC)/dbLookup.h \
$(INC)/search_queryIdx.h $(INC)/search_dbIdx.h $(INC)/alignments_multi.h\
$(INC)/alignments_dbIdx.h $(INC)/alignments_queryIdx.h  $(INC)/alignments.h\
$(INC)/blast.h $(INC)/bytepackGappedScoring.h $(INC)/descriptions.h\
$(INC)/encoding.h $(INC)/fasterGappedExtension_multi.h\
$(INC)/fasterGappedExtension.h $(INC)/gappedExtension_multi.h\
$(INC)/gappedExtension.h $(INC)/nuGappedScoring.h $(INC)/gappedScoring_multi.h\
$(INC)/gappedScoring.h $(INC)/global.h $(INC)/hitMatrix_multi.h\
$(INC)/hitMatrix.h $(INC)/karlin.h $(INC)/memBlocks.h $(INC)/memSingleBlock.h\
$(INC)/nucleotideLookup.h $(INC)/oldGappedScoring.h\
$(INC)/oldSemiGappedScoring.h $(INC)/parameters.h $(INC)/print.h\
$(INC)/PSSMatrix.h $(INC)/qPosList.h $(INC)/readFasta.h $(INC)/readFile.h\
$(INC)/scoreMatrix.h $(INC)/semiGappedScoring_multi.h $(INC)/semiGappedScoring.h\
$(INC)/statistics.h $(INC)/ungappedExtension_multi.h $(INC)/ungappedExtension.h\
$(INC)/wordLookupDFA_multi.h $(INC)/wordLookupDFA.h $(INC)/writeFile.h\
$(INC)/constants.h $(INC)/smithWatermanTraceback.h $(INC)/smithWatermanScoring.h\
$(INC)/tableGappedScoring.h $(INC)/vbyte.h $(INC)/unpack_multi.h $(INC)/unpack.h\
$(INC)/index.h $(INC)/postings.h $(INC)/hashcounter.h $(INC)/writedb.h\
$(INC)/readdb.h $(INC)/search.h $(INC)/wildcards.h $(INC)/dust.h $(INC)/seg.h\
$(NCBI_CORE)/ncbi_gapalign.h $(NCBI_CORE)/ncbi_stat.h $(NCBI_CORE)/ncbi_tree.h\
$(NCBI_CORE)/ncbi_traceback.h $(NCBI_CORE)/ncbi_matrix.h\
$(NCBI_CORE)/ncbi_composition_adjustment.h $(NCBI_CORE)/ncbi_filter.h\
$(INC)/alignments_ncbi.h $(INC)/dbLookup_compress.h\
$(NCBI_CORE)/blast_kappa.h

SRCFILES=$(SRC)/utils.c $(SRC)/merge.c $(SRC)/search_queryIdx.c\
$(SRC)/search_dbIdx.c $(SRC)/queryLookup.c  $(SRC)/dbLookup.c\
$(SRC)/alignments.c $(SRC)/alignments_multi.c $(SRC)/alignments_dbIdx.c\
$(SRC)/alignments_queryIdx.c $(SRC)/indexdb.c $(SRC)/mublastp.c $(SRC)/blast.c\
$(SRC)/bytepackGappedScoring.c $(SRC)/descriptions.c $(SRC)/encoding.c\
$(SRC)/fasterGappedExtension.c $(SRC)/fasterGappedExtension_multi.c\
$(SRC)/gappedExtension_multi.c $(SRC)/gappedExtension.c\
$(SRC)/gappedScoring_multi.c $(SRC)/gappedScoring.c $(SRC)/nuGappedScoring.c\
$(SRC)/global.c $(SRC)/hitMatrix_multi.c $(SRC)/hitMatrix.c $(SRC)/karlin.c\
$(SRC)/memBlocks.c $(SRC)/memSingleBlock.c $(SRC)/nucleotideLookup.c\
$(SRC)/oldGappedScoring.c $(SRC)/oldSemiGappedScoring.c $(SRC)/parameters.c\
$(SRC)/print.c $(SRC)/PSSMatrix.c $(SRC)/qPosList.c $(SRC)/readFasta.c\
$(SRC)/readFile.c $(SRC)/scoreMatrix.c $(SRC)/semiGappedScoring_multi.c\
$(SRC)/semiGappedScoring.c $(SRC)/statistics.c $(SRC)/ungappedExtension_multi.c\
$(SRC)/ungappedExtension.c $(SRC)/wordLookupDFA_multi.c $(SRC)/wordLookupDFA.c\
$(SRC)/writeFile.c $(SRC)/constants.c $(SRC)/smithWatermanTraceback.c\
$(SRC)/smithWatermanScoring.c $(SRC)/tableGappedScoring.c $(SRC)/vbyte.c\
$(SRC)/unpack_multi.c $(SRC)/unpack.c $(SRC)/index.c $(SRC)/postings.c\
$(SRC)/hashcounter.c $(SRC)/writedb.c $(SRC)/readdb.c $(SRC)/search.c\
$(SRC)/wildcards.c $(SRC)/dust.c $(SRC)/seg.c $(NCBI_CORE)/ncbi_gapalign.c\
$(NCBI_CORE)/ncbi_gapalign.c $(NCBI_CORE)/ncbi_tree.c\
$(NCBI_CORE)/ncbi_traceback.c $(NCBI_CORE)/ncbi_matrix.c\
$(NCBI_CORE)/ncbi_composition_adjustment.c $(NCBI_CORE)/ncbi_filter.c\
$(SRC)/alignments_ncbi.c $(SRC)/dbLookup_compress.c $(SRC)/sortdb.c\
$(NCBI_CORE)/blast_kappa.c

all: $(OBJ) mublastp formatdb indexdb sortdb #sampledb querySelector dbinfo blast

$(OBJ):
	mkdir -p $(OBJ)

#blast: $(OBJ)/blast.o $(OBJFILES)
	#$(CC) $(CFLAGS) -o blast $(OBJ)/blast.o $(OBJFILES) $(LDFLAGS) 
mublastp: $(OBJ)/mublastp.o $(OBJFILES)
	$(CC) $(CFLAGS) -o mublastp $(OBJ)/mublastp.o $(OBJFILES) $(LDFLAGS) 
ssearch: $(OBJ)/ssearch.o $(OBJFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -o ssearch $(OBJ)/ssearch.o $(OBJFILES)
blast-debug: $(SRC)/blast.c $(SRCFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -g -o blast-debug $(SRCFILES)
blast1: $(SRC)/blast.c $(SRCFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -DNO_STAGE2 -o blast1 $(SRCFILES)
blast12: $(SRC)/blast.c $(SRCFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -DNO_STAGE3 -o blast12 $(SRCFILES)
blast123: $(SRC)/blast.c $(SRCFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -DNO_STAGE4 -o blast123 $(SRCFILES)
blast-bitlookup: $(SRC)/blast.c $(SRCFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -DBITLOOKUP -o blast-bitlookup $(SRCFILES)
verboseBlast: $(SRC)/blast.c $(SRCFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -DVERBOSE -o verboseBlast $(SRCFILES)
formatdb: $(OBJ)/formatdb.o $(OBJFILES)
	$(CC) $(CFLAGS) -o formatdb $(OBJ)/formatdb.o $(OBJFILES) $(LDFLAGS) 
createindex: $(OBJ)/createindex.o $(OBJFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -o createindex $(OBJ)/createindex.o $(OBJFILES)
readdb: $(OBJ)/readdbApp.o $(OBJFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -o readdb $(OBJ)/readdbApp.o $(OBJFILES)
chooseWilds: $(OBJ)/chooseWilds.o $(OBJFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -o chooseWilds $(OBJ)/chooseWilds.o $(OBJFILES)
dust: $(OBJ)/dustApp.o $(OBJFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -o dust $(OBJ)/dustApp.o $(OBJFILES)
printDescription: $(OBJ)/printDescription.o $(OBJFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -o printDescription $(OBJ)/printDescription.o $(OBJFILES)
cluster: $(OBJ)/cluster.o $(OBJFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -o cluster $(OBJ)/cluster.o $(OBJFILES)
rsdb: $(OBJ)/rsdb.o $(OBJFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -o rsdb $(OBJ)/rsdb.o $(OBJFILES)

$(OBJ)/chooseWilds.o: $(SRC)/chooseWilds.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/chooseWilds.o $(SRC)/chooseWilds.c
$(OBJ)/readdbApp.o: $(SRC)/readdbApp.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/readdbApp.o $(SRC)/readdbApp.c
$(OBJ)/dustApp.o: $(SRC)/dustApp.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/dustApp.o $(SRC)/dustApp.c
$(OBJ)/printDescription.o: $(SRC)/printDescription.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/printDescription.o $(SRC)/printDescription.c
$(OBJ)/formatdb.o: $(SRC)/formatdb.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/formatdb.o $(SRC)/formatdb.c
$(OBJ)/createindex.o: $(SRC)/createindex.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/createindex.o $(SRC)/createindex.c
$(OBJ)/blast.o: $(SRC)/blast.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/blast.o $(SRC)/blast.c
$(OBJ)/cluster.o: $(SRC)/cluster.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/cluster.o $(SRC)/cluster.c
$(OBJ)/rsdb.o: $(SRC)/rsdb.c $(SRC)/identityAlign.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/rsdb.o $(SRC)/rsdb.c

$(OBJ)/alignments.o: $(SRC)/alignments.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/alignments.o $(SRC)/alignments.c
$(OBJ)/bytepackGappedScoring.o: $(SRC)/bytepackGappedScoring.c $(SRC)/fasterBytepackGappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/bytepackGappedScoring.o $(SRC)/bytepackGappedScoring.c
$(OBJ)/constants.o: $(SRC)/constants.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/constants.o $(SRC)/constants.c
$(OBJ)/descriptions.o: $(SRC)/descriptions.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/descriptions.o $(SRC)/descriptions.c
$(OBJ)/encoding.o: $(SRC)/encoding.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/encoding.o $(SRC)/encoding.c
$(OBJ)/gappedExtension.o: $(SRC)/gappedExtension.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/gappedExtension.o $(SRC)/gappedExtension.c
$(OBJ)/fasterGappedExtension.o: $(SRC)/fasterGappedExtension.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/fasterGappedExtension.o $(SRC)/fasterGappedExtension.c
$(OBJ)/gappedScoring.o: $(SRC)/gappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/gappedScoring.o $(SRC)/gappedScoring.c
$(OBJ)/gappedScoring_multi.o: $(SRC)/gappedScoring_multi.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/gappedScoring_multi.o $(SRC)/gappedScoring_multi.c
$(OBJ)/nuGappedScoring.o: $(SRC)/nuGappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/nuGappedScoring.o $(SRC)/nuGappedScoring.c
$(OBJ)/global.o: $(SRC)/global.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/global.o $(SRC)/global.c
$(OBJ)/hitMatrix.o: $(SRC)/hitMatrix.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/hitMatrix.o $(SRC)/hitMatrix.c
$(OBJ)/karlin.o: $(SRC)/karlin.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/karlin.o $(SRC)/karlin.c
$(OBJ)/memBlocks.o: $(SRC)/memBlocks.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/memBlocks.o $(SRC)/memBlocks.c
$(OBJ)/memSingleBlock.o: $(SRC)/memSingleBlock.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/memSingleBlock.o $(SRC)/memSingleBlock.c
$(OBJ)/nucleotideLookup.o: $(SRC)/nucleotideLookup.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/nucleotideLookup.o $(SRC)/nucleotideLookup.c
$(OBJ)/oldGappedScoring.o: $(SRC)/oldGappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/oldGappedScoring.o $(SRC)/oldGappedScoring.c
$(OBJ)/oldSemiGappedScoring.o: $(SRC)/oldSemiGappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/oldSemiGappedScoring.o $(SRC)/oldSemiGappedScoring.c
$(OBJ)/parameters.o: $(SRC)/parameters.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/parameters.o $(SRC)/parameters.c
$(OBJ)/print.o: $(SRC)/print.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/print.o $(SRC)/print.c
$(OBJ)/PSSMatrix.o: $(SRC)/PSSMatrix.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/PSSMatrix.o $(SRC)/PSSMatrix.c
$(OBJ)/qPosList.o: $(SRC)/qPosList.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/qPosList.o $(SRC)/qPosList.c
$(OBJ)/readFasta.o: $(SRC)/readFasta.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/readFasta.o $(SRC)/readFasta.c
$(OBJ)/readFile.o: $(SRC)/readFile.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/readFile.o $(SRC)/readFile.c
$(OBJ)/scoreMatrix.o: $(SRC)/scoreMatrix.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/scoreMatrix.o $(SRC)/scoreMatrix.c
$(OBJ)/semiGappedScoring.o: $(SRC)/semiGappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/semiGappedScoring.o $(SRC)/semiGappedScoring.c
$(OBJ)/statistics.o: $(SRC)/statistics.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/statistics.o $(SRC)/statistics.c
$(OBJ)/ungappedExtension.o: $(SRC)/ungappedExtension.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/ungappedExtension.o $(SRC)/ungappedExtension.c
$(OBJ)/wordLookupDFA.o: $(SRC)/wordLookupDFA.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/wordLookupDFA.o $(SRC)/wordLookupDFA.c
$(OBJ)/writeFile.o: $(SRC)/writeFile.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/writeFile.o $(SRC)/writeFile.c
$(OBJ)/smithWatermanTraceback.o: $(SRC)/smithWatermanTraceback.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/smithWatermanTraceback.o $(SRC)/smithWatermanTraceback.c
$(OBJ)/smithWatermanScoring.o: $(SRC)/smithWatermanScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/smithWatermanScoring.o $(SRC)/smithWatermanScoring.c
$(OBJ)/tableGappedScoring.o: $(SRC)/tableGappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/tableGappedScoring.o $(SRC)/tableGappedScoring.c
$(OBJ)/vbyte.o: $(SRC)/vbyte.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/vbyte.o $(SRC)/vbyte.c
$(OBJ)/unpack.o: $(SRC)/unpack.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/unpack.o $(SRC)/unpack.c
$(OBJ)/index.o: $(SRC)/index.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/index.o $(SRC)/index.c
$(OBJ)/hashcounter.o: $(SRC)/hashcounter.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/hashcounter.o $(SRC)/hashcounter.c
$(OBJ)/postings.o: $(SRC)/postings.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/postings.o $(SRC)/postings.c
$(OBJ)/writedb.o: $(SRC)/writedb.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/writedb.o $(SRC)/writedb.c
$(OBJ)/readdb.o: $(SRC)/readdb.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/readdb.o $(SRC)/readdb.c
$(OBJ)/search.o: $(SRC)/search.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/search.o $(SRC)/search.c
$(OBJ)/wildcards.o: $(SRC)/wildcards.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/wildcards.o $(SRC)/wildcards.c
$(OBJ)/ssearch.o: $(SRC)/blast.c $(HEADERFILES)
	$(CC) $(CFLAGS) -DSSEARCH -c -o $(OBJ)/ssearch.o $(SRC)/blast.c
$(OBJ)/dust.o: $(SRC)/dust.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/dust.o $(SRC)/dust.c
$(OBJ)/seg.o: $(SRC)/seg.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/seg.o $(SRC)/seg.c

$(OBJ)/querySelector.o: $(SRC)/querySelector.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/querySelector.o $(SRC)/querySelector.c
querySelector: $(OBJ)/querySelector.o $(OBJFILES)
	$(CC) $(CFLAGS) -o querySelector $(OBJ)/querySelector.o $(OBJFILES) $(LDFLAGS) 



$(OBJ)/search_dbIdx.o: $(SRC)/search_dbIdx.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/search_dbIdx.o $(SRC)/search_dbIdx.c
#$(OBJ)/search_queryIdx.o: $(SRC)/search_queryIdx.c $(HEADERFILES)
	#$(CC) $(CFLAGS) -c -o $(OBJ)/search_queryIdx.o $(SRC)/search_queryIdx.c
$(OBJ)/dbLookup.o: $(SRC)/dbLookup.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/dbLookup.o $(SRC)/dbLookup.c
$(OBJ)/queryLookup.o: $(SRC)/queryLookup.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/queryLookup.o $(SRC)/queryLookup.c
$(OBJ)/merge.o: $(SRC)/merge.c $(HEADERFILES)
	$(CCPLUS) $(CFLAGS) -c -o $(OBJ)/merge.o $(SRC)/merge.c

$(OBJ)/ungappedExtension_multi.o: $(SRC)/ungappedExtension_multi.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/ungappedExtension_multi.o $(SRC)/ungappedExtension_multi.c
$(OBJ)/alignments_multi.o: $(SRC)/alignments_multi.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/alignments_multi.o $(SRC)/alignments_multi.c
$(OBJ)/alignments_dbIdx.o: $(SRC)/alignments_dbIdx.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/alignments_dbIdx.o $(SRC)/alignments_dbIdx.c
#$(OBJ)/alignments_queryIdx.o: $(SRC)/alignments_queryIdx.c $(HEADERFILES)
	#$(CC) $(CFLAGS) -c -o $(OBJ)/alignments_queryIdx.o $(SRC)/alignments_queryIdx.c
$(OBJ)/unpack_multi.o: $(SRC)/unpack_multi.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/unpack_multi.o $(SRC)/unpack_multi.c

$(OBJ)/wordLookupDFA_multi.o: $(SRC)/wordLookupDFA_multi.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/wordLookupDFA_multi.o $(SRC)/wordLookupDFA_multi.c

$(OBJ)/fasterGappedExtension_multi.o: $(SRC)/fasterGappedExtension_multi.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/fasterGappedExtension_multi.o $(SRC)/fasterGappedExtension_multi.c

$(OBJ)/ncbi_gapalign.o: $(NCBI_CORE)/ncbi_gapalign.c $(HEADERFILES) 
	$(CC) $(CFLAGS) -c -o $(OBJ)/ncbi_gapalign.o $(NCBI_CORE)/ncbi_gapalign.c

$(OBJ)/blast_kappa.o: $(NCBI_CORE)/blast_kappa.c $(HEADERFILES) 
	$(CC) $(CFLAGS) -c -o $(OBJ)/blast_kappa.o $(NCBI_CORE)/blast_kappa.c

$(OBJ)/ncbi_stat.o: $(NCBI_CORE)/ncbi_stat.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/ncbi_stat.o $(NCBI_CORE)/ncbi_stat.c

$(OBJ)/ncbi_filter.o: $(NCBI_CORE)/ncbi_filter.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/ncbi_filter.o $(NCBI_CORE)/ncbi_filter.c

$(OBJ)/ncbi_traceback.o: $(NCBI_CORE)/ncbi_traceback.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/ncbi_traceback.o $(NCBI_CORE)/ncbi_traceback.c

$(OBJ)/ncbi_composition_adjustment.o: $(NCBI_CORE)/ncbi_composition_adjustment.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/ncbi_composition_adjustment.o $(NCBI_CORE)/ncbi_composition_adjustment.c

$(OBJ)/ncbi_matrix.o: $(NCBI_CORE)/ncbi_matrix.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/ncbi_matrix.o $(NCBI_CORE)/ncbi_matrix.c

$(OBJ)/ncbi_tree.o: $(NCBI_CORE)/ncbi_tree.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/ncbi_tree.o $(NCBI_CORE)/ncbi_tree.c

$(OBJ)/alignments_ncbi.o: $(SRC)/alignments_ncbi.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/alignments_ncbi.o $(SRC)/alignments_ncbi.c

$(OBJ)/dbLookup_compress.o: $(SRC)/dbLookup_compress.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/dbLookup_compress.o $(SRC)/dbLookup_compress.c

$(OBJ)/hitMatrix_multi.o: $(SRC)/hitMatrix_multi.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/hitMatrix_multi.o $(SRC)/hitMatrix_multi.c
$(OBJ)/utils.o: $(SRC)/utils.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/utils.o $(SRC)/utils.c
$(OBJ)/indexdb.o: $(SRC)/indexdb.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/indexdb.o $(SRC)/indexdb.c
indexdb: $(OBJ)/indexdb.o $(OBJFILES)
	$(CC) $(CFLAGS) -o indexdb $(OBJ)/indexdb.o $(OBJFILES) $(LDFLAGS)

$(OBJ)/sortdb.o: $(SRC)/sortdb.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/sortdb.o $(SRC)/sortdb.c

sortdb: $(OBJ)/sortdb.o $(OBJFILES)
	$(CC) $(CFLAGS) -o sortdb $(OBJ)/sortdb.o $(OBJFILES) $(LDFLAGS) 


$(OBJ)/dbinfo.o: $(SRC)/dbinfo.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/dbinfo.o $(SRC)/dbinfo.c

dbinfo: $(OBJ)/dbinfo.o $(OBJFILES)
	$(CC) $(CFLAGS) -o dbinfo $(OBJ)/dbinfo.o $(OBJFILES) $(LDFLAGS) 



$(OBJ)/sampledb.o: $(SRC)/sampledb.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/sampledb.o $(SRC)/sampledb.c

sampledb: $(OBJ)/sampledb.o $(OBJFILES)
	$(CC) $(CFLAGS) -o sampledb $(OBJ)/sampledb.o $(OBJFILES) $(LDFLAGS) 

$(OBJ)/gappedExtension_multi.o: $(SRC)/gappedExtension_multi.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/gappedExtension_multi.o $(SRC)/gappedExtension_multi.c

$(OBJ)/semiGappedScoring_multi.o: $(SRC)/semiGappedScoring_multi.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/semiGappedScoring_multi.o $(SRC)/semiGappedScoring_multi.c

$(OBJ)/mublastp.o: $(SRC)/mublastp.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/mublastp.o $(SRC)/mublastp.c

clean:
	rm -f mublastp sortdb blast indexdb formatdb $(OBJ)/*

