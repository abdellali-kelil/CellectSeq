//---------------------------------------------------------------------------
//
#include "Header.hpp"
#include "Pvalue.hpp"

//---------------------------------------------------------------------------
// STEP 0: Parse input arguments
void ParseAgrs(char *arg)
{
	// Verbose
	printf_s("\nSTEP00 : Parse input arguments");
	printf_s("\n==============================\n\n");

	// set nbr of opened files
	_setmaxstdio(8192);

	// Init random values genarator
	srand((unsigned int)time(NULL));

	if (_stricmp(arg, "01 CD151") == 0)
	{
		// antigen id
		PRJID = "01 CD151";

		// results destination
		PATH1 = ".\\RESULT\\01 CD151";

		// input NGS pools information
		DEMUX = ".\\INPUT\\demux - CD151.dat";

		// input positive/negative information
		TESTS = ".\\INPUT\\tests - CD151.dat";

		// input parameters
		NBRTHR =   10;			// Number of simultaneous threads
		MINOCC =    5;			// Min occurrences of sequence in positives
		MINNEG =    0.25f;		// Min frequency of sequences in negatives over positives
		MAXNWC =    0.35f;		// Max frequency of wildcards
		MINFRQ =    0.10f;		// Min frequency of amino acids at each position
		COHIDX =    1.20f;		// Min size effect for p-value significance (Cohen's index)
		CUTOFF =   1E-10f;		// pvalue cutoff

		// reads available
		READ1 = 1;
		READ2 = 1;
		READ3 = 1;
	}

	if (_stricmp(arg, "02 CA9") == 0)
	{
		// antigen id
		PRJID = "02 CA9";

		// results destination
		PATH1 = ".\\RESULT\\02 CA9";

		// input NGS pools information
		DEMUX = ".\\INPUT\\demux - CA9.dat";

		// input positive/negative information
		TESTS = ".\\INPUT\\tests - CA9.dat";

		// input parameters
		NBRTHR =   10;			// Number of simultaneous threads
		MINOCC =    5;			// Min occurrences of sequence in positives
		MINNEG =    0.25f;		// Min frequency of sequences in negatives over positives
		MAXNWC =    0.35f;		// Max frequency of wildcards
		MINFRQ =    0.10f;		// Min frequency of amino acids at each position
		COHIDX =    1.20f;		// Min size effect for p-value significance (Cohen's index)
		CUTOFF =   1E-10f;		// pvalue cutoff

		// reads available
		READ1 = 1;
		READ2 = 1;
		READ3 = 1;
	}

	if (_stricmp(arg, "03 ITGA11") == 0)
	{
		// antigen id
		PRJID = "03 ITGA11";

		// results destination
		PATH1 = ".\\RESULT\\03 ITGA11";

		// input NGS pools information
		DEMUX = ".\\INPUT\\demux - ITGA11.dat";

		// input positive/negative information
		TESTS = ".\\INPUT\\tests - ITGA11.dat";

		// input parameters
		NBRTHR =   10;			// Number of simultaneous threads
		MINOCC =    5;			// Min occurrences of sequence in positives
		MINNEG =    0.25f;		// Min frequency of sequences in negatives over positives
		MAXNWC =    0.35f;		// Max frequency of wildcards
		MINFRQ =    0.10f;		// Min frequency of amino acids at each position
		COHIDX =    1.20f;		// Min size effect for p-value significance (Cohen's index)
		CUTOFF =   1E-10f;		// pvalue cutoff

		// reads available
		READ1 = 1;
		READ2 = 1;
		READ3 = 1;
	}
}

//---------------------------------------------------------------------------
// STEP 1: Sequence frequencies
void Frequence()
{
	// Verbose
	printf_s("\nSTEP 1: Sequence frequencies over pools");
	printf_s("\n=======================================\n\n");

	// Nbr pools
	int NB1 = 0;

	// Buffer
	char BUF[MAXBUF] = "\0";
	char SEL[MAXBUF] = "\0";

	// Headers
	char *for1 = new char[MAXBUF]; *for1 = '\0';
	char *for2 = new char[MAXBUF]; *for2 = '\0';
	char *rev1 = new char[MAXBUF]; *rev1 = '\0';
	char *rev2 = new char[MAXBUF]; *rev2 = '\0';
	char *samp = new char[MAXBUF]; *samp = '\0';
	char *libr = new char[MAXBUF]; *libr = '\0';
	char *cont = new char[MAXBUF]; *cont = '\0';

	// Columns
	char **FOR1 = new char*[1];
	char **FOR2 = new char*[1];
	char **REV1 = new char*[1];
	char **REV2 = new char*[1];
	char **SAMP = new char*[1];
	char **LIBR = new char*[1];
	char **CONT = new char*[1];

	// Streams
	FILE *file0 = NULL;
	FILE *file1 = NULL;
	FILE *file2 = NULL;
	FILE *file3 = NULL;

	////////////////////////////////////////////////////////////////
	///////////////////////////// step 0 ///////////////////////////
	////////////////////////////////////////////////////////////////

	// Open demux file
	file0 = openfile(DEMUX, "rt");

	// Scroll over
	while (!feof(file0))
	{
		// Read current line										// Skeep comments
		if (fgets(BUF, MAXBUF, file0) == NULL) break;				if (strchr(BUF, '#') != NULL) continue;

		// Parse fileds
		if (sscanf_s(BUF, "%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\n]\n", for1, MAXBUF, for2, MAXBUF, rev1, MAXBUF, rev2, MAXBUF, samp, MAXBUF, libr, MAXBUF, cont, MAXBUF) != 7) printf_s("FATAL ERROR : wrong file format !!!\n"), exit(EXIT_FAILURE);

		// Inc samples
		NB1++;

		// Update samples
		FOR1 = Realloc(FOR1, NB1); FOR1[NB1 - 1] = new char[strlen(for1) + 1]; if (strcpy_s(FOR1[NB1 - 1], strlen(for1) + 1, for1) != 0) exit(EXIT_FAILURE);
		FOR2 = Realloc(FOR2, NB1); FOR2[NB1 - 1] = new char[strlen(for2) + 1]; if (strcpy_s(FOR2[NB1 - 1], strlen(for2) + 1, for2) != 0) exit(EXIT_FAILURE);
		REV1 = Realloc(REV1, NB1); REV1[NB1 - 1] = new char[strlen(rev1) + 1]; if (strcpy_s(REV1[NB1 - 1], strlen(rev1) + 1, rev1) != 0) exit(EXIT_FAILURE);
		REV2 = Realloc(REV2, NB1); REV2[NB1 - 1] = new char[strlen(rev2) + 1]; if (strcpy_s(REV2[NB1 - 1], strlen(rev2) + 1, rev2) != 0) exit(EXIT_FAILURE);
		SAMP = Realloc(SAMP, NB1); SAMP[NB1 - 1] = new char[strlen(samp) + 1]; if (strcpy_s(SAMP[NB1 - 1], strlen(samp) + 1, samp) != 0) exit(EXIT_FAILURE);
		LIBR = Realloc(LIBR, NB1); LIBR[NB1 - 1] = new char[strlen(libr) + 1]; if (strcpy_s(LIBR[NB1 - 1], strlen(libr) + 1, libr) != 0) exit(EXIT_FAILURE);
		CONT = Realloc(CONT, NB1); CONT[NB1 - 1] = new char[strlen(cont) + 1]; if (strcpy_s(CONT[NB1 - 1], strlen(cont) + 1, cont) != 0) exit(EXIT_FAILURE);
	}

	// Close stream
	fclose(file0);

	////////////////////////////////////////////////////////////////
	///////////////////////////// step 1 ///////////////////////////
	////////////////////////////////////////////////////////////////

	// Check doubled samples
	for (int i=0; i<NB1; i++)
	{
		for (int j=0; j<NB1; j++)
		{
			// Skeep if same sample
			if (i == j) continue;

			// Check sample labels
			if (_stricmp(FOR1[i], FOR1[j]) == 0 && _stricmp(REV1[i], REV1[j]) == 0) printf("FATAL ERROR : wrong pool labels !!!\n"), exit(EXIT_FAILURE);

			// Check sample labels
			if (_stricmp(FOR2[i], FOR2[j]) == 0 && _stricmp(REV2[i], REV2[j]) == 0) printf("FATAL ERROR : wrong pool barcodes !!!\n"), exit(EXIT_FAILURE);
		}
	}

	////////////////////////////////////////////////////////////////
	///////////////////////////// step 2 ///////////////////////////
	////////////////////////////////////////////////////////////////

	// Scroll over
	for (int i=0; i<NB1; i++)
	{
		// Parsers
		int occ; bool sep; char *buf;

		// sequences
		char sel[MAXBUF] = "\0"; char seq[NBREAD][MAXBUF] = {"\0", "\0", "\0", "\0"};

		// suffix tree
		CTree *Tree = NULL;

		// Load library
		Library(LIBR[i]);

		// Get selection pool
		Selection(SEL, FOR1[i], REV1[i], LIBR[i]);

		//////////////////////////
		///////// step 1 /////////
		//////////////////////////

		// Prepare Suffix Tree
		Tree        = new CTree;
		Tree->IDX   = -1;
		Tree->RES   = '#';
		Tree->OCC   = NULL;
		Tree->BRO   = NULL;
		Tree->SON   = NULL;

		// Get sequences input filename
		file1 = openfile(PATH1, "STEP01", SEL, "dat", "rt");

		// Parsers
		int nbr1 = 0;

		// Scroll over
		while (!feof(file1))
		{
			// Verbose
			if (++nbr1 % 1000 == 0) printf_s("%-15s %-20s step 1 %10d\r", PRJID, SEL, nbr1);

			// Read current line
			if (fgets(BUF, MAXBUF, file1) == NULL) break;

			// skeep header
			if (strchr(BUF, '#') != NULL) continue;

			// Parser
			buf = BUF;

			// Read aa sequences
			if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%s", seq[0], MAXBUF); buf = strchr(buf, '\t') + 1; }
			if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%s", seq[1], MAXBUF); buf = strchr(buf, '\t') + 1; }
			if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%s", seq[2], MAXBUF); buf = strchr(buf, '\t') + 1; }
			if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%s", seq[3], MAXBUF); buf = strchr(buf, '\t') + 1; }

			// Read occurrences
			sscanf_s(buf, "%d", &occ); buf = strchr(buf, '\t') + 1;

			// Init sequence
			BUF[0] = '\0';

			// Separator
			sep = false;

			// Print to corresponding output file
			if (READ3 == 1 && CDRL3 == 1) { if (sep) strcat_s(BUF, MAXBUF, "\t"); strcat_s(BUF, MAXBUF, seq[0]); sep = true; }
			if (READ2 == 1 && CDRH1 == 1) { if (sep) strcat_s(BUF, MAXBUF, "\t"); strcat_s(BUF, MAXBUF, seq[1]); sep = true; }
			if (READ2 == 1 && CDRH2 == 1) { if (sep) strcat_s(BUF, MAXBUF, "\t"); strcat_s(BUF, MAXBUF, seq[2]); sep = true; }
			if (READ1 == 1 && CDRH3 == 1) { if (sep) strcat_s(BUF, MAXBUF, "\t"); strcat_s(BUF, MAXBUF, seq[3]); sep = true; }

			// Add sequence to suffix tree
			SetTree0(Tree, BUF, length(BUF), 0, NB1, 0, 0);
		}

		// Close stream
		fclose(file1);

		// Verbose
		printf_s("%-15s %-20s step 1 %10d\n", PRJID, SEL, nbr1);

		//////////////////////////
		///////// step 2 /////////
		//////////////////////////

		// Scroll over
		for (int j=0; j<NB1; j++)
		{
			// Load library
			if (_stricmp(LIBR[i], LIBR[j]) != 0) continue;

			// Vernbose
			int nbr2 = 0;

			// Get selection pool
			Selection(sel, FOR1[j], REV1[j], LIBR[j]);

			// Get filename
			file2 = openfile(PATH1, "STEP01", sel, "dat", "rt");

			// Scroll over
			while (!feof(file2))
			{
				// Verbose
				if (++nbr2 % 1000 == 0) printf_s("%-15s %-20s step 2 %10d %s\r", PRJID, SEL, nbr2, sel);

				// Read current line
				if (fgets(BUF, MAXBUF, file2) == NULL) break;

				// skeep header
				if (strchr(BUF, '#') != NULL) continue;

				// Parser
				buf = BUF;

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%s", seq[0], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%s", seq[1], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%s", seq[2], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%s", seq[3], MAXBUF); buf = strchr(buf, '\t') + 1; }

				// Read occurrences
				sscanf_s(buf, "%d", &occ); buf = strchr(buf, '\t') + 1;

				// Init sequence
				BUF[0] = '\0';

				// Separator
				sep = false;

				// Print to corresponding output file
				if (READ3 == 1 && CDRL3 == 1) { if (sep) strcat_s(BUF, MAXBUF, "\t"); strcat_s(BUF, MAXBUF, seq[0]); sep = true; }
				if (READ2 == 1 && CDRH1 == 1) { if (sep) strcat_s(BUF, MAXBUF, "\t"); strcat_s(BUF, MAXBUF, seq[1]); sep = true; }
				if (READ2 == 1 && CDRH2 == 1) { if (sep) strcat_s(BUF, MAXBUF, "\t"); strcat_s(BUF, MAXBUF, seq[2]); sep = true; }
				if (READ1 == 1 && CDRH3 == 1) { if (sep) strcat_s(BUF, MAXBUF, "\t"); strcat_s(BUF, MAXBUF, seq[3]); sep = true; }

				// Write sequences
				GetTree1(Tree, BUF, length(BUF), 0, j, occ);
			}

			// Close stream
			fclose(file2);

			// Verbose
			printf_s("%-15s %-20s step 2 %10d %s\n", PRJID, SEL, nbr2, sel);
		}

		//////////////////////////
		///////// step 3 /////////
		//////////////////////////

		// Get output filename
		file3 = openfile(".\\TEMP\\TMP1.dat", "wt");

		// Write sequences with occurrences
		GetTree0(Tree, BUF, 0, NB1, NULL, file3);

		// Close stream
		fclose(file3);

		// free Tree
		EmptyTree(Tree);

		// Sort sequences
		sort(".\\TEMP\\TMP1.dat", ".\\TEMP\\TMP2.dat", "\t", 1, 1 + i + READ1*CDRH3 + READ2*CDRH1 + READ2*CDRH2 + READ3*CDRL3, "gr");

		// Sort and copy sequences
		sprintf_s(BUF, MAXBUF, "copy /b .\\TEMP\\TMP2.dat \"%s\\STEP02\\%s.dat\" > nul", PATH1, SEL); system(BUF);

		// Delete temporary files
		system("del /F /Q .\\TEMP\\*.* > nul");

		// Verbose
		printf_s("\n");
	}

	// free memory
	delete[] for1;
	delete[] for2;
	delete[] rev1;
	delete[] rev2;
	delete[] samp;
	delete[] libr;
	delete[] cont;

	// free memory
	for (int i=0; i<NB1; i++)
	{
		delete[] FOR1[i];
		delete[] FOR2[i];
		delete[] REV1[i];
		delete[] REV2[i];
		delete[] SAMP[i];
		delete[] LIBR[i];
		delete[] CONT[i];
	}

	// free memory
	delete[] FOR1;
	delete[] FOR2;
	delete[] REV1;
	delete[] REV2;
	delete[] SAMP;
	delete[] LIBR;
	delete[] CONT;

	// Verbose
	printf_s("\n");
}

//---------------------------------------------------------------------------
// STEP 2: Discover motifs
void Positives()
{
	// Verbose
	printf_s("\nSTEP 2: Discover motifs in positive pool");
	printf_s("\n========================================\n\n");

	// Nbr pools
	int NB1 = 0;

	// Buffer
	char BUF[MAXBUF] = "\0";
	char COM[MAXBUF] = "\0";

	// Headers
	char *for1 = new char[MAXBUF]; *for1 = '\0';
	char *for2 = new char[MAXBUF]; *for2 = '\0';
	char *rev1 = new char[MAXBUF]; *rev1 = '\0';
	char *rev2 = new char[MAXBUF]; *rev2 = '\0';
	char *samp = new char[MAXBUF]; *samp = '\0';
	char *libr = new char[MAXBUF]; *libr = '\0';
	char *cont = new char[MAXBUF]; *cont = '\0';

	// Columns
	char **FOR1 = new char*[1];
	char **FOR2 = new char*[1];
	char **REV1 = new char*[1];
	char **REV2 = new char*[1];
	char **SAMP = new char*[1];
	char **LIBR = new char*[1];
	char **CONT = new char*[1];

	// Streams
	FILE *file0;
	FILE *file1;
	FILE *file2;
	FILE *file3;
	FILE *file4;
	FILE *file5;
	FILE *file6;
	FILE *file7;
	FILE *file8;
	FILE *file9;
	FILE *file10;
	FILE *file11;
	FILE *file12;
	FILE *file13;

	////////////////////////////////////////////////////////////////
	///////////////////////////// step 0 ///////////////////////////
	////////////////////////////////////////////////////////////////

	// Open demux file
	file0 = openfile(DEMUX, "rt");

	// Scroll over
	while (!feof(file0))
	{
		// Read current line										// Skeep comments
		if (fgets(BUF, MAXBUF, file0) == NULL) break;				if (strchr(BUF, '#') != NULL) continue;

		// Parse fileds
		if (sscanf_s(BUF, "%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\n]\n", for1, MAXBUF, for2, MAXBUF, rev1, MAXBUF, rev2, MAXBUF, samp, MAXBUF, libr, MAXBUF, cont, MAXBUF) != 7) printf_s("FATAL ERROR : wrong file format !!!\n"), exit(EXIT_FAILURE);

		// Inc samples
		NB1++;

		// Update samples
		FOR1 = Realloc(FOR1, NB1); FOR1[NB1 - 1] = new char[strlen(for1) + 1]; if (strcpy_s(FOR1[NB1 - 1], strlen(for1) + 1, for1) != 0) exit(EXIT_FAILURE);
		FOR2 = Realloc(FOR2, NB1); FOR2[NB1 - 1] = new char[strlen(for2) + 1]; if (strcpy_s(FOR2[NB1 - 1], strlen(for2) + 1, for2) != 0) exit(EXIT_FAILURE);
		REV1 = Realloc(REV1, NB1); REV1[NB1 - 1] = new char[strlen(rev1) + 1]; if (strcpy_s(REV1[NB1 - 1], strlen(rev1) + 1, rev1) != 0) exit(EXIT_FAILURE);
		REV2 = Realloc(REV2, NB1); REV2[NB1 - 1] = new char[strlen(rev2) + 1]; if (strcpy_s(REV2[NB1 - 1], strlen(rev2) + 1, rev2) != 0) exit(EXIT_FAILURE);
		SAMP = Realloc(SAMP, NB1); SAMP[NB1 - 1] = new char[strlen(samp) + 1]; if (strcpy_s(SAMP[NB1 - 1], strlen(samp) + 1, samp) != 0) exit(EXIT_FAILURE);
		LIBR = Realloc(LIBR, NB1); LIBR[NB1 - 1] = new char[strlen(libr) + 1]; if (strcpy_s(LIBR[NB1 - 1], strlen(libr) + 1, libr) != 0) exit(EXIT_FAILURE);
		CONT = Realloc(CONT, NB1); CONT[NB1 - 1] = new char[strlen(cont) + 1]; if (strcpy_s(CONT[NB1 - 1], strlen(cont) + 1, cont) != 0) exit(EXIT_FAILURE);
	}

	// Close stream
	fclose(file0);

	////////////////////////////////////////////////////////////////
	///////////////////////////// step 1 ///////////////////////////
	////////////////////////////////////////////////////////////////

	// Check doubled samples
	for (int i=0; i<NB1; i++)
	{
		for (int j=0; j<NB1; j++)
		{
			// Skeep if same sample
			if (i == j) continue;

			// Check sample labels
			if (_stricmp(FOR1[i], FOR1[j]) == 0 && _stricmp(REV1[i], REV1[j]) == 0) printf("FATAL ERROR : wrong pool labels !!!\n"), exit(EXIT_FAILURE);

			// Check sample labels
			if (_stricmp(FOR2[i], FOR2[j]) == 0 && _stricmp(REV2[i], REV2[j]) == 0) printf("FATAL ERROR : wrong pool barcodes !!!\n"), exit(EXIT_FAILURE);
		}
	}

	////////////////////////////////////////////////////////////////
	///////////////////////////// step 2 ///////////////////////////
	////////////////////////////////////////////////////////////////

	// Nbr tests
	int NB2 = 0;

	// Parser, positives and negatives
	char *pos = new char[MAXBUF]; char **POS = new char*[1];
	char *neg = new char[MAXBUF]; char **NEG = new char*[1];
	char *lib = new char[MAXBUF]; char **LIB = new char*[1];
	char *lab = new char[MAXBUF]; char **LAB = new char*[1];

	// Open tests file
	file1 = openfile(TESTS, "rt");

	// Scroll over
	while (!feof(file1))
	{
		// Read current line										// Skeep comments
		if (fgets(BUF, MAXBUF, file1) == NULL) break;				if (strchr(BUF, '#') != NULL) continue;

		// Parse fileds
		if (sscanf_s(BUF, "%[^\t]\t%[^\t]\t%[^\t]\t%[^\n]\n", pos, MAXBUF, neg, MAXBUF, lib, MAXBUF, lab, MAXBUF) != 4) printf_s("FATAL ERROR : wrong file format !!!\n"), exit(EXIT_FAILURE);

		// inc
		NB2++;

		POS = Realloc(POS, NB2); POS[NB2 - 1] = new char[strlen(pos) + 1]; if (strcpy_s(POS[NB2 - 1], strlen(pos) + 1, pos) != 0) exit(EXIT_FAILURE);
		NEG = Realloc(NEG, NB2); NEG[NB2 - 1] = new char[strlen(neg) + 1]; if (strcpy_s(NEG[NB2 - 1], strlen(neg) + 1, neg) != 0) exit(EXIT_FAILURE);
		LIB = Realloc(LIB, NB2); LIB[NB2 - 1] = new char[strlen(lib) + 1]; if (strcpy_s(LIB[NB2 - 1], strlen(lib) + 1, lib) != 0) exit(EXIT_FAILURE);
		LAB = Realloc(LAB, NB2); LAB[NB2 - 1] = new char[strlen(lab) + 1]; if (strcpy_s(LAB[NB2 - 1], strlen(lab) + 1, lab) != 0) exit(EXIT_FAILURE);
	}

	// Close stream
	fclose(file1);

	////////////////////////////////////////////////////////////////
	///////////////////////////// step 3 ///////////////////////////
	////////////////////////////////////////////////////////////////

	for (int i=0; i<NB2; i++)
	{
		// Parsers
		int occ; int idx; bool sep; char *buf;

		// sequences
		char sel[MAXBUF] = "\0"; char SEQ[MAXBUF] = "\0"; char **seq = Alloc2D(NBREAD, MAXBUF, '\0');

		// lengths
		int *len = Alloc1D(NBREAD, 0);

		// Suffix tree
		CTree *Tree = NULL;

		// Load library
		Library(LIB[i]);

		// samples combination
		Combination(COM, POS[i], NEG[i]);

		//////////////////////////
		///////// step 1 /////////
		//////////////////////////

		// Type of each sample, positive/negative
		int *TYP = new int[NB1];

		// Set each sample
		for (int j=0; j<NB1; j++)
		{
			// Get sequences input filename
			sprintf_s(BUF, MAXBUF, "%s_%s", FOR1[j], REV1[j]);

			// Check if positive/negative/none
			if (strstr(POS[i], BUF) != NULL) TYP[j] = +1; else if (strstr(NEG[i], BUF) != NULL) TYP[j] = -1; else TYP[j] = 0;
		}

		//////////////////////////
		///////// step 2 /////////
		//////////////////////////

		// Total number of all sequences
		int *ALL = Alloc1D(NB1, 0);

		// Find total number of sequences, in positives and negatives
		for (int j=0; j<NB1; j++)
		{
			// Skeep if not positive or negative
			if (TYP[j] == 0) continue;

			// Check library
			if (_stricmp(LIB[i], LIBR[j]) != 0) printf("FATAL ERROR : wrong library !!!\n"), exit(EXIT_FAILURE);

			// Verbose
			int nbr0 = 0;

			// Get selection pool
			Selection(sel, FOR1[j], REV1[j], LIBR[j]);

			// Get filename
			file2 = openfile(PATH1, "STEP01", sel, "dat", "rt");

			// Scroll over
			while (!feof(file2))
			{
				// Verbose
				if (++nbr0 % 1000 == 0) printf_s("%-15s P(%s) step 0 %10d %s\r", PRJID, POS[i], nbr0, sel);

				// Read the current line from the input file
				if (fgets(BUF, MAXBUF, file2) == NULL) break;

				// skeep header
				if (strchr(BUF, '#') != NULL) continue;

				// Parse buffer
				buf = BUF;

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) buf = strchr(buf, '\t') + 1;
				if (READ2 == 1 && CDRH1 == 1) buf = strchr(buf, '\t') + 1;
				if (READ2 == 1 && CDRH2 == 1) buf = strchr(buf, '\t') + 1;
				if (READ1 == 1 && CDRH3 == 1) buf = strchr(buf, '\t') + 1;

				// Read occurrences
				sscanf_s(buf, "%d", &occ);

				// Update occurrences
				ALL[j] += occ;
			}

			// Close stream
			fclose(file2);

			// Verbose
			printf_s("%-15s P(%s) step 0 %10d %s\n", PRJID, POS[i], nbr0, sel);
		}

		//////////////////////////
		///////// step 3 /////////
		//////////////////////////

		// Get output file for current variable sequences
		file3 = openfile(".\\TEMP\\TMP1.dat", "wt");
		file4 = openfile(".\\TEMP\\TMP2.dat", "wt");
		file5 = openfile(".\\TEMP\\TMP3.dat", "wt");

		// Scroll over
		for (int j=0; j<NB1; j++)
		{
			// Skeep if not positive
			if (TYP[j] == 0) continue;

			// Check library
			if (_stricmp(LIB[i], LIBR[j]) != 0) printf("FATAL ERROR : wrong library !!!\n"), exit(EXIT_FAILURE);

			// Verbose
			int nbr1 = 0;

			// Get selection pool
			Selection(sel, FOR1[j], REV1[j], LIBR[j]);

			// Get filename
			file6 = openfile(PATH1, "STEP02", sel, "dat", "rt");

			// Scroll over
			while (!feof(file6))
			{
				// Verbose
				if (++nbr1 % 1000 == 0) printf_s("%-15s P(%s) step 1 %10d %s\r", PRJID, POS[i], nbr1, sel);

				// Read current line
				if (fgets(BUF, MAXBUF, file6) == NULL) break;

				// trim buffer
				sscanf_s(BUF, "%[^\n]\n", BUF, MAXBUF);

				// Parser
				buf = BUF;

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%s", seq[0], MAXBUF); buf = strchr(buf, '\t') + 1; len[0] = length(seq[0]); }
				if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%s", seq[1], MAXBUF); buf = strchr(buf, '\t') + 1; len[1] = length(seq[1]); }
				if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%s", seq[2], MAXBUF); buf = strchr(buf, '\t') + 1; len[2] = length(seq[2]); }
				if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%s", seq[3], MAXBUF); buf = strchr(buf, '\t') + 1; len[3] = length(seq[3]); }

				// Init occurrences in positives and negatives
				int occ1 = +INT_MAX; int all1 = +INT_MAX; double frq1 = +DBL_MAX;
				int occ2 = -INT_MAX; int all2 = -INT_MAX; double frq2 = -DBL_MAX;

				// Read occurrences
				for (int k=0; k<NB1; k++)
				{
					// Read occurrences
					sscanf_s(buf, "%d", &occ); buf = strchr(buf, '\t') + 1;

					switch (TYP[k])
					{
						case +1: if (frq1 > (double) occ / ALL[k]) {occ1 = occ; all1 = ALL[k]; frq1 = (double) occ / ALL[k]; } break;
						case -1: if (frq2 < (double) occ / ALL[k]) {occ2 = occ; all2 = ALL[k]; frq2 = (double) occ / ALL[k]; } break;
					}
				}

				// Check type and min occ
				if (TYP[j] == 1 && occ1 > 0)
				{
					// Check occurrences
					while (occ1 >= MINOCC)
					{
						// Check stop codons
						if (strchr(BUF, '*') != NULL) break;

						// Check wildcards
						if (strchr(BUF, WCSYMB) != NULL) break;

						// Get sequences lengths
						if (READ3 == 1 && CDRL3 == 1) fprintf_s(file3, "%d\t", len[0]);
						if (READ2 == 1 && CDRH1 == 1) fprintf_s(file3, "%d\t", len[1]);
						if (READ2 == 1 && CDRH2 == 1) fprintf_s(file3, "%d\t", len[2]);
						if (READ1 == 1 && CDRH3 == 1) fprintf_s(file3, "%d\t", len[3]);

						// write sequence
						fprintf_s(file3, "%s\n", BUF); break;
					}

					// Get sequences lengths
					if (READ3 == 1 && CDRL3 == 1) fprintf_s(file4, "%d\t", len[0]);
					if (READ2 == 1 && CDRH1 == 1) fprintf_s(file4, "%d\t", len[1]);
					if (READ2 == 1 && CDRH2 == 1) fprintf_s(file4, "%d\t", len[2]);
					if (READ1 == 1 && CDRH3 == 1) fprintf_s(file4, "%d\t", len[3]);

					// Write positive sequences
					fprintf_s(file4, "%s\n", BUF);
				}

				// Get sequences lengths
				if (READ3 == 1 && CDRL3 == 1) fprintf_s(file5, "%d\t", len[0]);
				if (READ2 == 1 && CDRH1 == 1) fprintf_s(file5, "%d\t", len[1]);
				if (READ2 == 1 && CDRH2 == 1) fprintf_s(file5, "%d\t", len[2]);
				if (READ1 == 1 && CDRH3 == 1) fprintf_s(file5, "%d\t", len[3]);

				// Write positive and negative sequences
				fprintf_s(file5, "%s\n", BUF);
			}

			// Close stream
			fclose(file6);

			// Verbose
			printf_s("%-15s P(%s) step 1 %10d %s\n", PRJID, POS[i], nbr1, sel);
		}

		// Close stream
		fclose(file3);
		fclose(file4);
		fclose(file5);

		// sort sequences by length
		switch ((READ3 == 1 && CDRL3 == 1) + (READ2 == 1 && CDRH1 == 1) + (READ2 == 1 && CDRH2 == 1) + (READ1 == 1 && CDRH3 == 1))
		{
			case 4 :
			{
				system("sort -T .\\TEMP\\ --parallel=10 -t\"\t\" -k1g -k2g -k3g -k4g .\\TEMP\\TMP1.dat | uniq > .\\TEMP\\TMP4.dat");
				system("sort -T .\\TEMP\\ --parallel=10 -t\"\t\" -k1g -k2g -k3g -k4g .\\TEMP\\TMP2.dat | uniq > .\\TEMP\\TMP5.dat");
				system("sort -T .\\TEMP\\ --parallel=10 -t\"\t\" -k1g -k2g -k3g -k4g .\\TEMP\\TMP3.dat | uniq > .\\TEMP\\TMP6.dat"); break;
			}

			case 3 :
			{
				system("sort -T .\\TEMP\\ --parallel=10 -t\"\t\" -k1g -k2g -k3g .\\TEMP\\TMP1.dat | uniq > .\\TEMP\\TMP4.dat");
				system("sort -T .\\TEMP\\ --parallel=10 -t\"\t\" -k1g -k2g -k3g .\\TEMP\\TMP2.dat | uniq > .\\TEMP\\TMP5.dat");
				system("sort -T .\\TEMP\\ --parallel=10 -t\"\t\" -k1g -k2g -k3g .\\TEMP\\TMP3.dat | uniq > .\\TEMP\\TMP6.dat"); break;
			}

			case 2 :
			{
				system("sort -T .\\TEMP\\ --parallel=10 -t\"\t\" -k1g -k2g .\\TEMP\\TMP1.dat | uniq > .\\TEMP\\TMP4.dat");
				system("sort -T .\\TEMP\\ --parallel=10 -t\"\t\" -k1g -k2g .\\TEMP\\TMP2.dat | uniq > .\\TEMP\\TMP5.dat");
				system("sort -T .\\TEMP\\ --parallel=10 -t\"\t\" -k1g -k2g .\\TEMP\\TMP3.dat | uniq > .\\TEMP\\TMP6.dat"); break;
			}

			case 1 :
			{
				system("sort -T .\\TEMP\\ --parallel=10 -t\"\t\" -k1g .\\TEMP\\TMP1.dat | uniq > .\\TEMP\\TMP4.dat");
				system("sort -T .\\TEMP\\ --parallel=10 -t\"\t\" -k1g .\\TEMP\\TMP2.dat | uniq > .\\TEMP\\TMP5.dat");
				system("sort -T .\\TEMP\\ --parallel=10 -t\"\t\" -k1g .\\TEMP\\TMP3.dat | uniq > .\\TEMP\\TMP6.dat"); break;
			}
		}

		// Output file for current sequences						// Close stream
		file12 = openfile(PATH1, "STEP03", COM, "dat", "wt");		fclose(file12);

		//////////////////////////
		///////// step 4 /////////
		//////////////////////////

		// Get lengths of variable sequences
		int   NB3 = 0;
		int  *NB4 = new int[1];
		int  *NB5 = new int[1];
		int  *LEN = new int[1];
		int **CDR = new int*[NBREAD];

		// init position	// Verbose
		__int64 pos0 = 0;	int nbr2 = 0;

		// Read input file
		file7 = openfile(".\\TEMP\\TMP4.dat", "rt");

		// Get positive sequences
		while (!feof(file7))
		{
			// Verbose
			printf_s("%-15s P(%s) step 2 %10d %s\r", PRJID, POS[i], ++nbr2, sel);

			// Get current sequence
			if (fgets(BUF, MAXBUF, file7) == NULL) break;

			// Parser
			buf = BUF; idx = 0;

			// Read aa sequences
			if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%d", &len[0]); buf = strchr(buf, '\t') + 1; }
			if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%d", &len[1]); buf = strchr(buf, '\t') + 1; }
			if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%d", &len[2]); buf = strchr(buf, '\t') + 1; }
			if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%d", &len[3]); buf = strchr(buf, '\t') + 1; }

			// check first
			if (NB3 == 0) goto Add;

			// Current idx
			idx = 0;

			// Check length for each read
			if (READ3 == 1 && CDRL3 == 1) if (CDR[NB3-1][idx++] != len[0]) goto Add;
			if (READ2 == 1 && CDRH1 == 1) if (CDR[NB3-1][idx++] != len[1]) goto Add;
			if (READ2 == 1 && CDRH2 == 1) if (CDR[NB3-1][idx++] != len[2]) goto Add;
			if (READ1 == 1 && CDRH3 == 1) if (CDR[NB3-1][idx++] != len[3]) goto Add;

			// inc nbr
			NB4[NB3-1] += 1;

			// next
			continue;

		Add :	{
				// Inc sequences
				NB3 += 1;

				// Resize sequences lengths
				LEN = Realloc(LEN, NB3);
				CDR = Realloc(CDR, NB3);
				NB4 = Realloc(NB4, NB3);
				NB5 = Realloc(NB5, NB3);

				// Init sequences occurrences and lengths
				NB4[NB3-1] = 1; NB5[NB3-1] = 0; CDR[NB3-1] = new int[NBREAD];

				// Get occurrences
				LEN[NB3-1] = READ1*CDRH3 + READ2*CDRH1 + READ2*CDRH2 + READ3*CDRL3 - 1;

				// Curent idx
				idx = 0;

				// Update length
				if (READ3 == 1 && CDRL3 == 1) { CDR[NB3-1][idx++] = len[0]; LEN[NB3-1] += len[0]; }
				if (READ2 == 1 && CDRH1 == 1) { CDR[NB3-1][idx++] = len[1]; LEN[NB3-1] += len[1]; }
				if (READ2 == 1 && CDRH2 == 1) { CDR[NB3-1][idx++] = len[2]; LEN[NB3-1] += len[2]; }
				if (READ1 == 1 && CDRH3 == 1) { CDR[NB3-1][idx++] = len[3]; LEN[NB3-1] += len[3]; }

				// Read input file
				file8 = openfile(".\\TEMP\\TMP6.dat", "rt");

				// back to previous position
				_fseeki64(file8, pos0, SEEK_SET);

				// Get positive
				while (!feof(file8))
				{
					// get current position
					pos0 = _ftelli64(file8);

					// Get current sequence
					if (fgets(BUF, MAXBUF, file8) == NULL) break;

					// Parser
					buf = BUF; idx = 0;

					// Read aa sequences
					if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%d", &len[0]); buf = strchr(buf, '\t') + 1; if (CDR[NB3-1][idx] > len[0]) continue; if (CDR[NB3-1][idx] < len[0]) break; idx++; }
					if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%d", &len[1]); buf = strchr(buf, '\t') + 1; if (CDR[NB3-1][idx] > len[1]) continue; if (CDR[NB3-1][idx] < len[1]) break; idx++; }
					if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%d", &len[2]); buf = strchr(buf, '\t') + 1; if (CDR[NB3-1][idx] > len[2]) continue; if (CDR[NB3-1][idx] < len[2]) break; idx++; }
					if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%d", &len[3]); buf = strchr(buf, '\t') + 1; if (CDR[NB3-1][idx] > len[3]) continue; if (CDR[NB3-1][idx] < len[3]) break; idx++; }

					// inc number
					NB5[NB3-1] += 1;
				}

				// Close stream
				fclose(file8);
			}
		}

		// Close stream
		fclose(file7);

		// Verbose
		printf_s("\n");

		//////////////////////////
		///////// step 5 /////////
		//////////////////////////

		// init positions
		__int64 pos1 = 0;
		__int64 pos2 = 0;
		__int64 pos3 = 0;

		// Scroll over
		for (int j=0; j<NB3; j++)
		{
			//////////////////
			///// step 1 /////
			//////////////////

			// Init rwas data containers
			char **SEQ1 = Alloc2D(NB4[j], LEN[j]+1, '\0'); int NBR1 = 0;

			// Read input file
			file9 = openfile(".\\TEMP\\TMP4.dat", "rt");

			// back to previous position
			_fseeki64(file9, pos1, SEEK_SET);

			// Get positive sequences, 
			while (!feof(file9))
			{
				// get current position
				pos1 = _ftelli64(file9);

				// Get current sequence
				if (fgets(BUF, MAXBUF, file9) == NULL) break;

				// trim buffer
				sscanf_s(BUF, "%[^\n]\n", BUF, MAXBUF);

				// Parser
				buf = BUF; idx = 0;

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%d", &len[0]); buf = strchr(buf, '\t') + 1; if (CDR[j][idx] > len[0]) continue; if (CDR[j][idx] < len[0]) break; idx++; }
				if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%d", &len[1]); buf = strchr(buf, '\t') + 1; if (CDR[j][idx] > len[1]) continue; if (CDR[j][idx] < len[1]) break; idx++; }
				if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%d", &len[2]); buf = strchr(buf, '\t') + 1; if (CDR[j][idx] > len[2]) continue; if (CDR[j][idx] < len[2]) break; idx++; }
				if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%d", &len[3]); buf = strchr(buf, '\t') + 1; if (CDR[j][idx] > len[3]) continue; if (CDR[j][idx] < len[3]) break; idx++; }

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%s", seq[0], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%s", seq[1], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%s", seq[2], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%s", seq[3], MAXBUF); buf = strchr(buf, '\t') + 1; }

				// Init sequence
				SEQ[0] = '\0';

				// Separator
				sep = false;

				// Print to corresponding output file
				if (READ3 == 1 && CDRL3 == 1) { if (sep) strcat_s(SEQ, MAXBUF, "_"); strcat_s(SEQ, MAXBUF, seq[0]); sep = true; }
				if (READ2 == 1 && CDRH1 == 1) { if (sep) strcat_s(SEQ, MAXBUF, "_"); strcat_s(SEQ, MAXBUF, seq[1]); sep = true; }
				if (READ2 == 1 && CDRH2 == 1) { if (sep) strcat_s(SEQ, MAXBUF, "_"); strcat_s(SEQ, MAXBUF, seq[2]); sep = true; }
				if (READ1 == 1 && CDRH3 == 1) { if (sep) strcat_s(SEQ, MAXBUF, "_"); strcat_s(SEQ, MAXBUF, seq[3]); sep = true; }

				// Get sequence
				sprintf_s(SEQ1[NBR1++], strlen(SEQ)+1, "%s", SEQ);
			}

			// Close stream
			fclose(file9);

			// check error
			if (NB4[j] != NBR1) exit(EXIT_FAILURE);

			//////////////////
			///// step 2 /////
			//////////////////

			// Prepare Suffix Tree
			Tree        = new CTree;
			Tree->IDX   = -1;
			Tree->RES   = '#';
			Tree->OCC   = NULL;
			Tree->BRO   = NULL;
			Tree->SON   = NULL;

			// Read input file
			file10 = openfile(".\\TEMP\\TMP5.dat", "rt");

			// back to previous position
			_fseeki64(file10, pos2, SEEK_SET);

			// Get positive sequences, scroll over all sequences
			while (!feof(file10))
			{
				// get current position
				pos2 = _ftelli64(file10);

				// Read the current line from the input file
				if (fgets(BUF, MAXBUF, file10) == NULL) break;

				// Parser
				buf = BUF; idx = 0;

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%d", &len[0]); buf = strchr(buf, '\t') + 1; if (CDR[j][idx] > len[0]) continue; if (CDR[j][idx] < len[0]) break; idx++; }
				if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%d", &len[1]); buf = strchr(buf, '\t') + 1; if (CDR[j][idx] > len[1]) continue; if (CDR[j][idx] < len[1]) break; idx++; }
				if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%d", &len[2]); buf = strchr(buf, '\t') + 1; if (CDR[j][idx] > len[2]) continue; if (CDR[j][idx] < len[2]) break; idx++; }
				if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%d", &len[3]); buf = strchr(buf, '\t') + 1; if (CDR[j][idx] > len[3]) continue; if (CDR[j][idx] < len[3]) break; idx++; }

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%s", seq[0], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%s", seq[1], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%s", seq[2], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%s", seq[3], MAXBUF); buf = strchr(buf, '\t') + 1; }

				// Init sequence
				SEQ[0] = '\0';

				// Separator
				sep = false;

				// Print to corresponding output file
				if (READ3 == 1 && CDRL3 == 1) { if (sep) strcat_s(SEQ, MAXBUF, "_"); strcat_s(SEQ, MAXBUF, seq[0]); sep = true; }
				if (READ2 == 1 && CDRH1 == 1) { if (sep) strcat_s(SEQ, MAXBUF, "_"); strcat_s(SEQ, MAXBUF, seq[1]); sep = true; }
				if (READ2 == 1 && CDRH2 == 1) { if (sep) strcat_s(SEQ, MAXBUF, "_"); strcat_s(SEQ, MAXBUF, seq[2]); sep = true; }
				if (READ1 == 1 && CDRH3 == 1) { if (sep) strcat_s(SEQ, MAXBUF, "_"); strcat_s(SEQ, MAXBUF, seq[3]); sep = true; }

				// If negative
				for (int l=0; l<NB1; l++)
				{
					// Read current occurrence
					sscanf_s(buf, "%d", &occ); buf = strchr(buf, '\t') + 1;

					// Update suffix tree
					if (TYP[l] == 1) SetTree2(Tree, SEQ, length(SEQ), 0, NB1, l, occ);
				}
			}

			// Close stream
			fclose(file10);

			// sort suffix tree
			SortTree1(Tree);

			//////////////////
			///// step 3 /////
			//////////////////

			// Init rwas data containers
			char **SEQ2 = Alloc2D(NB5[j], LEN[j]+1, '\0'); int **OCC2 = Alloc2D(NB5[j], NB1, 0); int NBR2 = 0;

			// Read input file
			file11 = openfile(".\\TEMP\\TMP6.dat", "rt");

			// back to previous position
			_fseeki64(file11, pos3, SEEK_SET);

			// Get positive sequences, 
			while (!feof(file11))
			{
				// get current position
				pos3 = _ftelli64(file11);

				// Get current sequence
				if (fgets(BUF, MAXBUF, file11) == NULL) break;

				// trim buffer
				sscanf_s(BUF, "%[^\n]\n", BUF, MAXBUF);

				// Parser
				buf = BUF; idx = 0;

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%d", &len[0]); buf = strchr(buf, '\t') + 1; if (CDR[j][idx] > len[0]) continue; if (CDR[j][idx] < len[0]) break; idx++; }
				if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%d", &len[1]); buf = strchr(buf, '\t') + 1; if (CDR[j][idx] > len[1]) continue; if (CDR[j][idx] < len[1]) break; idx++; }
				if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%d", &len[2]); buf = strchr(buf, '\t') + 1; if (CDR[j][idx] > len[2]) continue; if (CDR[j][idx] < len[2]) break; idx++; }
				if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%d", &len[3]); buf = strchr(buf, '\t') + 1; if (CDR[j][idx] > len[3]) continue; if (CDR[j][idx] < len[3]) break; idx++; }

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%s", seq[0], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%s", seq[1], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%s", seq[2], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%s", seq[3], MAXBUF); buf = strchr(buf, '\t') + 1; }

				// Init sequence
				SEQ[0] = '\0';

				// Separator
				sep = false;

				// Print to corresponding output file
				if (READ3 == 1 && CDRL3 == 1) { if (sep) strcat_s(SEQ, MAXBUF, "_"); strcat_s(SEQ, MAXBUF, seq[0]); sep = true; }
				if (READ2 == 1 && CDRH1 == 1) { if (sep) strcat_s(SEQ, MAXBUF, "_"); strcat_s(SEQ, MAXBUF, seq[1]); sep = true; }
				if (READ2 == 1 && CDRH2 == 1) { if (sep) strcat_s(SEQ, MAXBUF, "_"); strcat_s(SEQ, MAXBUF, seq[2]); sep = true; }
				if (READ1 == 1 && CDRH3 == 1) { if (sep) strcat_s(SEQ, MAXBUF, "_"); strcat_s(SEQ, MAXBUF, seq[3]); sep = true; }

				// Inc nbr of sequences
				NBR2 += 1;

				// Get sequence
				sprintf_s(SEQ2[NBR2-1], strlen(SEQ)+1, "%s", SEQ);

				// Update occurrences
				OCC2[NBR2-1] = Alloc1D(NB1, 0);

				// If negative
				for (int l=0; l<NB1; l++)
				{
					// Read current occurrence
					sscanf_s(buf, "%d", &occ); buf = strchr(buf, '\t') + 1;

					// Check for seqs errs
					OCC2[NBR2-1][l] = occ;
				}
			}

			// Close stream
			fclose(file11);

			// check error
			if (NB5[j] != NBR2) exit(EXIT_FAILURE);

			//////////////////
			///// step 4 /////
			//////////////////

			// Init threads ids
			unsigned  iThread;
			HANDLE   *hThread = new HANDLE[NBRTHR];
			CThread **cThread = new CThread*[NBRTHR];

			// Scroll to max length
			for (int k=1; k <= LEN[j]; k++)
			{
				// Verbose
				printf_s("%-15s P(%s) step 3  itr %2d/%-2d len %2d/%d\r", PRJID, POS[i], j + 1, NB3, k, LEN[j]);

				// Init threads
				for (int l=0; l<NBRTHR; l++)
				{
					// Current thread
					cThread[l]       = new CThread[1];	// Current thread
					cThread[l]->THR  = l;				// Current thread id
					cThread[l]->CUR  = k;				// Current motifs length
					cThread[l]->DIM  = NB1;				// Dimension, number of experience
					cThread[l]->CDR  = CDR[j];			// Variable sequences lengths
					cThread[l]->TYP  = TYP;				// Type (positive, negative or none)
					cThread[l]->ALL  = ALL;				// Total nbr of sequences (occurrences)
					cThread[l]->NB1  = NBR1;			// Total number of sequences in positive pool
					cThread[l]->NB2  = NBR2;			// Total number of sequences in negative pool
					cThread[l]->OCC2 = OCC2;			// Total Number of sequences in pool1 by Pool
					cThread[l]->SEQ1 = SEQ1;			// Sequences in positive pool
					cThread[l]->SEQ2 = SEQ2;			// sequences in negative pool
					cThread[l]->TREE = Tree;			// Suffix tree
				}

				// Launch threads
				for (int l=0; l<NBRTHR; l++)
				{
					// Create handle
					hThread[l] = (HANDLE)_beginthreadex(	NULL,				// Thread security
															0,					// Thread stack size
															&DiscoverMotifs,	// Thread starting address
															(void *)cThread[l],	// Thread start argument
															DETACHED_PROCESS,	// Thread suspended
															&iThread);			// Thread ID

					// Check if thread successfuly created
					if (hThread[l] == 0 || hThread[l] == INVALID_HANDLE_VALUE) printf("Thread Creation Failed\n"), exit(EXIT_FAILURE);
				}

				// Wait for threads to run
				WaitForMultipleObjectsEx(NBRTHR, hThread, TRUE, INFINITE, FALSE);

				// Close all of thread handles
				for (int l=0; l<NBRTHR; l++) CloseHandle(hThread[l]);
			}

			// free memory
			for (int k=0; k<NBRTHR; k++)
			delete[] cThread[k];
			delete[] hThread;

			// free Tree
			EmptyTree(Tree);

			// Verbose
			printf_s("\n");

			//////////////////
			///// step 5 /////
			//////////////////

			// Get output file for current variable sequences
			file12 = openfile(PATH1, "STEP03", COM, "dat", "at");

			// Scroll over
			for (int k=0; k<NBRTHR; k++)
			{
				// Open current output file											// Open file for reading
				sprintf_s(BUF, MAXBUF, ".\\TEMP\\TMP9_%d_%d.dat", LEN[j], k);		file13 = openfile(BUF, "rt");

				while (!feof(file13))
				{
					// Read current line from input file
					if (fgets(BUF, MAXBUF, file13) == NULL) break;

					// trim buffer
					sscanf_s(BUF, "%[^\n]\n", BUF, MAXBUF);

					// Separator
					sep = false;

					// Print to corresponding output file
					if (READ3 == 1 && CDRL3 == 1) { if (sep) if ((buf = strchr(BUF, '_')) != NULL) *buf = '\t'; sep = true; }
					if (READ2 == 1 && CDRH1 == 1) { if (sep) if ((buf = strchr(BUF, '_')) != NULL) *buf = '\t'; sep = true; }
					if (READ2 == 1 && CDRH2 == 1) { if (sep) if ((buf = strchr(BUF, '_')) != NULL) *buf = '\t'; sep = true; }
					if (READ1 == 1 && CDRH3 == 1) { if (sep) if ((buf = strchr(BUF, '_')) != NULL) *buf = '\t'; sep = true; }

					// Write current line to output file
					fprintf_s(file12, "%s\n", BUF);
				}

				// Close stream
				fclose(file13);
			}

			// Close stream
			fclose(file12);

			// free memory
			for (int k=0; k<NBR1; k++)
			{
				delete[] SEQ1[k];
			}

			// free memory
			for (int k=0; k<NBR2; k++)
			{
				delete[] SEQ2[k];
				delete[] OCC2[k];
			}

			// free memory
			delete[] SEQ1;
			delete[] SEQ2;
			delete[] OCC2;

			// Delete temporary files
			system("del /F /Q .\\TEMP\\TMP7*.dat > nul");
			system("del /F /Q .\\TEMP\\TMP8*.dat > nul");
			system("del /F /Q .\\TEMP\\TMP9*.dat > nul");
		}

		// free memory
		for (int j=0; j<NBREAD; j++)
		delete[] seq[j];
		delete[] seq;

		// free memory
		for (int j=0; j<NB3; j++)
		delete[] CDR[j];
		delete[] CDR;
		delete[] LEN;
		delete[] NB4;
		delete[] NB5;

		// free memory
		delete[] TYP;
		delete[] ALL;

		// Delete temporary files
		system("del /F /Q .\\TEMP\\TMP*.dat > nul");

		// Verbose
		printf_s("---------------------------------------------------------------------------------------------------\n\n");
	}

	// free memory
	delete[] for1;
	delete[] for2;
	delete[] rev1;
	delete[] rev2;
	delete[] samp;
	delete[] libr;
	delete[] cont;

	// free memory
	for (int i=0; i<NB1; i++)
	{
		delete[] FOR1[i];
		delete[] FOR2[i];
		delete[] REV1[i];
		delete[] REV2[i];
		delete[] SAMP[i];
		delete[] LIBR[i];
		delete[] CONT[i];
	}

	// free memory
	delete[] FOR1;
	delete[] FOR2;
	delete[] REV1;
	delete[] REV2;
	delete[] SAMP;
	delete[] LIBR;
	delete[] CONT;

	// free memory
	for (int i=0; i<NB2; i++)
	{
		delete[] POS[i];
		delete[] NEG[i];
		delete[] LIB[i];
		delete[] LAB[i];
	}

	// free memory
	delete[] POS; delete[] pos;
	delete[] NEG; delete[] neg;
	delete[] LIB; delete[] lib;
	delete[] LAB; delete[] lab;
}

//---------------------------------------------------------------------------
// STEP 3: Predict sequences
void Predicted()
{
	// Verbose
	printf_s("\nSTEP 3: Predict specific sequences");
	printf_s("\n==================================\n\n");

	// Nbr pools
	int NB1 = 0;

	// Buffer
	char BUF[MAXBUF] = "\0";
	char COM[MAXBUF] = "\0";

	// Headers
	char *for1 = new char[MAXBUF]; *for1 = '\0';
	char *for2 = new char[MAXBUF]; *for2 = '\0';
	char *rev1 = new char[MAXBUF]; *rev1 = '\0';
	char *rev2 = new char[MAXBUF]; *rev2 = '\0';
	char *samp = new char[MAXBUF]; *samp = '\0';
	char *libr = new char[MAXBUF]; *libr = '\0';
	char *cont = new char[MAXBUF]; *cont = '\0';

	// Columns
	char **FOR1 = new char*[1];
	char **FOR2 = new char*[1];
	char **REV1 = new char*[1];
	char **REV2 = new char*[1];
	char **SAMP = new char*[1];
	char **LIBR = new char*[1];
	char **CONT = new char*[1];

	// Streams
	FILE *file0;
	FILE *file1;
	FILE *file2;
	FILE *file3;
	FILE *file4;
	FILE *file5;
	FILE *file6;
	FILE *file7;
	FILE *file8;
	FILE *file9;
	FILE *file10;
	FILE *file11;
	FILE *file12;
	FILE *file13;
	FILE *file14;
	FILE *file15;
	FILE *file16;
	FILE *file17;
	FILE *file18;

	////////////////////////////////////////////////////////////////
	///////////////////////////// step 0 ///////////////////////////
	////////////////////////////////////////////////////////////////

	// Open demux file
	file0 = openfile(DEMUX, "rt");

	// Scroll over
	while (!feof(file0))
	{
		// Read current line										// Skeep comments
		if (fgets(BUF, MAXBUF, file0) == NULL) break;				if (strchr(BUF, '#') != NULL) continue;

		// Parse fileds
		if (sscanf_s(BUF, "%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\n]\n", for1, MAXBUF, for2, MAXBUF, rev1, MAXBUF, rev2, MAXBUF, samp, MAXBUF, libr, MAXBUF, cont, MAXBUF) != 7) printf_s("FATAL ERROR : wrong file format !!!\n"), exit(EXIT_FAILURE);

		// Inc samples
		NB1++;

		// Update samples
		FOR1 = Realloc(FOR1, NB1); FOR1[NB1 - 1] = new char[strlen(for1) + 1]; if (strcpy_s(FOR1[NB1 - 1], strlen(for1) + 1, for1) != 0) exit(EXIT_FAILURE);
		FOR2 = Realloc(FOR2, NB1); FOR2[NB1 - 1] = new char[strlen(for2) + 1]; if (strcpy_s(FOR2[NB1 - 1], strlen(for2) + 1, for2) != 0) exit(EXIT_FAILURE);
		REV1 = Realloc(REV1, NB1); REV1[NB1 - 1] = new char[strlen(rev1) + 1]; if (strcpy_s(REV1[NB1 - 1], strlen(rev1) + 1, rev1) != 0) exit(EXIT_FAILURE);
		REV2 = Realloc(REV2, NB1); REV2[NB1 - 1] = new char[strlen(rev2) + 1]; if (strcpy_s(REV2[NB1 - 1], strlen(rev2) + 1, rev2) != 0) exit(EXIT_FAILURE);
		SAMP = Realloc(SAMP, NB1); SAMP[NB1 - 1] = new char[strlen(samp) + 1]; if (strcpy_s(SAMP[NB1 - 1], strlen(samp) + 1, samp) != 0) exit(EXIT_FAILURE);
		LIBR = Realloc(LIBR, NB1); LIBR[NB1 - 1] = new char[strlen(libr) + 1]; if (strcpy_s(LIBR[NB1 - 1], strlen(libr) + 1, libr) != 0) exit(EXIT_FAILURE);
		CONT = Realloc(CONT, NB1); CONT[NB1 - 1] = new char[strlen(cont) + 1]; if (strcpy_s(CONT[NB1 - 1], strlen(cont) + 1, cont) != 0) exit(EXIT_FAILURE);
	}

	// Close stream
	fclose(file0);

	////////////////////////////////////////////////////////////////
	///////////////////////////// step 1 ///////////////////////////
	////////////////////////////////////////////////////////////////

	// Check doubled samples
	for (int i=0; i<NB1; i++)
	{
		for (int j=0; j<NB1; j++)
		{
			// Skeep if same sample
			if (i == j) continue;

			// Check sample labels
			if (_stricmp(FOR1[i], FOR1[j]) == 0 && _stricmp(REV1[i], REV1[j]) == 0) printf("FATAL ERROR : wrong pool labels !!!\n"), exit(EXIT_FAILURE);

			// Check sample labels
			if (_stricmp(FOR2[i], FOR2[j]) == 0 && _stricmp(REV2[i], REV2[j]) == 0) printf("FATAL ERROR : wrong pool barcodes !!!\n"), exit(EXIT_FAILURE);
		}
	}

	////////////////////////////////////////////////////////////////
	///////////////////////////// step 2 ///////////////////////////
	////////////////////////////////////////////////////////////////

	// Nbr tests
	int NB2 = 0;

	// Parser, positives and negatives
	char *pos = new char[MAXBUF]; char **POS = new char*[1];
	char *neg = new char[MAXBUF]; char **NEG = new char*[1];
	char *lib = new char[MAXBUF]; char **LIB = new char*[1];
	char *lab = new char[MAXBUF]; char **LAB = new char*[1];

	// Open tests file
	file1 = openfile(TESTS, "rt");

	// Scroll over
	while (!feof(file1))
	{
		// Read current line										// Skeep comments
		if (fgets(BUF, MAXBUF, file1) == NULL) break;				if (strchr(BUF, '#') != NULL) continue;

		// Parse fileds
		if (sscanf_s(BUF, "%[^\t]\t%[^\t]\t%[^\t]\t%[^\n]\n", pos, MAXBUF, neg, MAXBUF, lib, MAXBUF, lab, MAXBUF) != 4) printf_s("FATAL ERROR : wrong file format !!!\n"), exit(EXIT_FAILURE);

		// inc
		NB2++;

		POS = Realloc(POS, NB2); POS[NB2 - 1] = new char[strlen(pos) + 1]; if (strcpy_s(POS[NB2 - 1], strlen(pos) + 1, pos) != 0) exit(EXIT_FAILURE);
		NEG = Realloc(NEG, NB2); NEG[NB2 - 1] = new char[strlen(neg) + 1]; if (strcpy_s(NEG[NB2 - 1], strlen(neg) + 1, neg) != 0) exit(EXIT_FAILURE);
		LIB = Realloc(LIB, NB2); LIB[NB2 - 1] = new char[strlen(lib) + 1]; if (strcpy_s(LIB[NB2 - 1], strlen(lib) + 1, lib) != 0) exit(EXIT_FAILURE);
		LAB = Realloc(LAB, NB2); LAB[NB2 - 1] = new char[strlen(lab) + 1]; if (strcpy_s(LAB[NB2 - 1], strlen(lab) + 1, lab) != 0) exit(EXIT_FAILURE);
	}

	// Close stream
	fclose(file1);

	////////////////////////////////////////////////////////////////
	///////////////////////////// step 3 ///////////////////////////
	////////////////////////////////////////////////////////////////

	for (int i=0; i<NB2; i++)
	{
		// Parsers
		int occ; bool sep; char *buf;

		// sequences
		char sel[MAXBUF] = "\0"; char SEQ0[MAXBUF] = "\0"; char SEQ1[MAXBUF] = "\0"; char SEQ2[MAXBUF] = "\0"; char SEQ3[MAXBUF] = "\0"; char **seq = Alloc2D(NBREAD, MAXBUF, '\0');

		// suffix trees
		CTree *Tree1 = NULL;
		CTree *Tree2 = NULL;

		// Load library
		Library(LIB[i]);

		// samples combination
		Combination(COM, POS[i], NEG[i]);

		//////////////////////////
		///////// step 1 /////////
		//////////////////////////

		// Type of each sample, positive/negative
		int *TYP = new int[NB1];

		// Set each sample
		for (int j=0; j<NB1; j++)
		{
			// Get sequences input filename
			sprintf_s(BUF, MAXBUF, "%s_%s", FOR1[j], REV1[j]);

			// Check if positive/negative/none
			if (strstr(POS[i], BUF) != NULL) TYP[j] = +1; else if (strstr(NEG[i], BUF) != NULL) TYP[j] = -1; else TYP[j] = 0;
		}

		//////////////////////////
		///////// step 2 /////////
		//////////////////////////

		// Total number of all sequences
		int *ALL = Alloc1D(NB1, 0);

		// Find total number of sequences, in positives and negatives
		for (int j=0; j<NB1; j++)
		{
			// Skeep if not positive or negative
			if (TYP[j] == 0) continue;

			// Check library
			if (_stricmp(LIB[i], LIBR[j]) != 0) printf("FATAL ERROR : wrong library !!!\n"), exit(EXIT_FAILURE);

			// Verbose
			int nbr0 = 0;

			// Get selection pool
			Selection(sel, FOR1[j], REV1[j], LIBR[j]);

			// Get filename
			file2 = openfile(PATH1, "STEP01", sel, "dat", "rt");

			// Scroll over
			while (!feof(file2))
			{
				// Verbose
				if (++nbr0 % 1000 == 0) printf_s("%-15s P(%s) step 0 %10d %s\r", PRJID, POS[i], nbr0, sel);

				// Read the current line from the input file
				if (fgets(BUF, MAXBUF, file2) == NULL) break;

				// skeep header
				if (strchr(BUF, '#') != NULL) continue;

				// Parse buffer
				buf = BUF;

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) buf = strchr(buf, '\t') + 1;
				if (READ2 == 1 && CDRH1 == 1) buf = strchr(buf, '\t') + 1;
				if (READ2 == 1 && CDRH2 == 1) buf = strchr(buf, '\t') + 1;
				if (READ1 == 1 && CDRH3 == 1) buf = strchr(buf, '\t') + 1;

				// Read occurrences
				sscanf_s(buf, "%d", &occ);

				// Update occurrences
				ALL[j] += occ;
			}

			// Close stream
			fclose(file2);

			// Verbose
			printf_s("%-15s P(%s) step 0 %10d %s\n", PRJID, POS[i], nbr0, sel);
		}

		//////////////////////////
		///////// step 3 /////////
		//////////////////////////

		// Get output file for current variable sequences
		file3 = openfile(".\\TEMP\\TMP1.dat", "wt");

		// Scroll over
		for (int j=0; j<NB1; j++)
		{
			// Skeep if not positive
			if (TYP[j] != 1) continue;

			// Check library
			if (_stricmp(LIB[i], LIBR[j]) != 0) printf("FATAL ERROR : wrong library !!!\n"), exit(EXIT_FAILURE);

			// Verbose
			int nbr1 = 0;

			// Get selection pool
			Selection(sel, FOR1[j], REV1[j], LIBR[j]);

			// Get filename
			file4 = openfile(PATH1, "STEP02", sel, "dat", "rt");

			// Scroll over
			while (!feof(file4))
			{
				// Verbose
				if (++nbr1 % 1000 == 0) printf_s("%-15s P(%s) step 1 %10d %s\r", PRJID, POS[i], nbr1, sel);

				// Read current line
				if (fgets(BUF, MAXBUF, file4) == NULL) break;

				// trim buffer
				sscanf_s(BUF, "%[^\n]\n", BUF, MAXBUF);

				// Check wildcards and stop codons
				if (strchr(BUF, WCSYMB) != NULL || strchr(BUF, '*') != NULL) continue;

				// Parser
				buf = BUF;

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%s", seq[0], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%s", seq[1], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%s", seq[2], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%s", seq[3], MAXBUF); buf = strchr(buf, '\t') + 1; }

				// Init occurrences in positives and negatives
				int occ1 = +INT_MAX; int all1 = +INT_MAX; double frq1 = +DBL_MAX;
				int occ2 = -INT_MAX; int all2 = -INT_MAX; double frq2 = -DBL_MAX;

				// Read occurrences
				for (int k=0; k<NB1; k++)
				{
					// Read occurrences
					sscanf_s(buf, "%d", &occ); buf = strchr(buf, '\t') + 1;

					switch (TYP[k])
					{
						case +1: if (frq1 > (double) occ / ALL[k]) {occ1 = occ; all1 = ALL[k]; frq1 = (double) occ / ALL[k]; } break;
						case -1: if (frq2 < (double) occ / ALL[k]) {occ2 = occ; all2 = ALL[k]; frq2 = (double) occ / ALL[k]; } break;
					}
				}

				// Separator
				sep = false;

				// write sequences
				if (READ3 == 1 && CDRL3 == 1) { if (sep) fprintf_s(file3, "\t"); fprintf_s(file3, "%s", seq[0]); sep = true; }
				if (READ2 == 1 && CDRH1 == 1) { if (sep) fprintf_s(file3, "\t"); fprintf_s(file3, "%s", seq[1]); sep = true; }
				if (READ2 == 1 && CDRH2 == 1) { if (sep) fprintf_s(file3, "\t"); fprintf_s(file3, "%s", seq[2]); sep = true; }
				if (READ1 == 1 && CDRH3 == 1) { if (sep) fprintf_s(file3, "\t"); fprintf_s(file3, "%s", seq[3]); sep = true; }

				// write candidate
				fprintf_s(file3, "\t%d\t%d\t%d\t%d\t%.4e\t%.4e\n", occ1, occ2, all1, all2, frq1, frq2);
			}

			// Close stream;
			fclose(file4);

			// Verbose
			printf_s("%-15s P(%s) step 1 %10d %s\n", PRJID, POS[i], nbr1, sel);
		}

		// Close stream
		fclose(file3);

		// Remove duplicat sequences
		system("sort -T .\\TEMP\\ --parallel=10 .\\TEMP\\TMP1.dat | uniq > .\\TEMP\\TMP2.dat");

		//////////////////////////
		///////// step 4 /////////
		//////////////////////////

		// Init suffix tree
		Tree1        = new CTree;
		Tree1->IDX   = -1;
		Tree1->RES   = '#';
		Tree1->OCC   = NULL;
		Tree1->BRO   = NULL;
		Tree1->SON   = NULL;

		// Verbose
		int nbr2 = 0;

		// Get output file for current variable sequences
		file5 = openfile(".\\TEMP\\TMP2.dat", "rt");

		while (!feof(file5))
		{
			// Verbose
			if (++nbr2 % 1000 == 0) printf_s("%-15s P(%s) step 2 %10d\r", PRJID, POS[i], nbr2);

			// Read current line
			if (fgets(BUF, MAXBUF, file5) == NULL) break;

			// Parser
			buf = BUF;

			// Read aa sequences
			if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%s", seq[0], MAXBUF); buf = strchr(buf, '\t') + 1; }
			if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%s", seq[1], MAXBUF); buf = strchr(buf, '\t') + 1; }
			if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%s", seq[2], MAXBUF); buf = strchr(buf, '\t') + 1; }
			if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%s", seq[3], MAXBUF); buf = strchr(buf, '\t') + 1; }

			// Init sequence
			SEQ0[0] = '\0';

			// Separator
			sep = false;

			// Print to corresponding output file
			if (READ3 == 1 && CDRL3 == 1) { if (sep) strcat_s(SEQ0, MAXBUF, "\t"); strcat_s(SEQ0, MAXBUF, seq[0]); sep = true; }
			if (READ2 == 1 && CDRH1 == 1) { if (sep) strcat_s(SEQ0, MAXBUF, "\t"); strcat_s(SEQ0, MAXBUF, seq[1]); sep = true; }
			if (READ2 == 1 && CDRH2 == 1) { if (sep) strcat_s(SEQ0, MAXBUF, "\t"); strcat_s(SEQ0, MAXBUF, seq[2]); sep = true; }
			if (READ1 == 1 && CDRH3 == 1) { if (sep) strcat_s(SEQ0, MAXBUF, "\t"); strcat_s(SEQ0, MAXBUF, seq[3]); sep = true; }

			// Read	occurrence
			sscanf_s(buf, "%d", &occ); buf = strchr(buf, '\t') + 1;

			// Update suffix tree
			SetTree0(Tree1, SEQ0, length(SEQ0), 0, 2, 0, occ);

			// Read occurrence
			sscanf_s(buf, "%d", &occ); buf = strchr(buf, '\t') + 1;

			// Update suffix tree
			SetTree0(Tree1, SEQ0, length(SEQ0), 0, 2, 1, occ);
		}

		// Close stream;
		fclose(file5);

		// Verbose
		printf_s("%-15s P(%s) step 2 %10d\n", PRJID, POS[i], nbr2);

		//////////////////////////
		///////// step 5 /////////
		//////////////////////////

		// Get motifs filename
		sprintf_s(BUF, MAXBUF, "%s\\STEP03\\%s.dat", PATH1, COM);

		// Check motifs file
		if(_access_s(BUF, 0) == 0)
		{
			// Init suffix tree
			Tree2        = new CTree;
			Tree2->IDX   = -1;
			Tree2->RES   = '#';
			Tree2->OCC   = NULL;
			Tree2->BRO   = NULL;
			Tree2->SON   = NULL;

			// Suffix Tree and parsers
			int nbr3 = 0;

			// Get input file
			file6 = openfile(PATH1, "STEP03", COM, "dat", "rt");

			while (!feof(file6))
			{
				// Verbose
				if (++nbr3 % 1000 == 0) printf_s("%-15s P(%s) step 3 %10d\r", PRJID, POS[i], nbr3);

				// Read current line
				if (fgets(BUF, MAXBUF, file6) == NULL) break;

				// Parser
				buf = BUF;

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%s", seq[0], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%s", seq[1], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%s", seq[2], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%s", seq[3], MAXBUF); buf = strchr(buf, '\t') + 1; }

				// Init sequence
				SEQ0[0] = '\0';

				// Separator
				sep = false;

				// Print to corresponding output file
				if (READ3 == 1 && CDRL3 == 1) { if (sep) strcat_s(SEQ0, MAXBUF, "\t"); strcat_s(SEQ0, MAXBUF, seq[0]); sep = true; }
				if (READ2 == 1 && CDRH1 == 1) { if (sep) strcat_s(SEQ0, MAXBUF, "\t"); strcat_s(SEQ0, MAXBUF, seq[1]); sep = true; }
				if (READ2 == 1 && CDRH2 == 1) { if (sep) strcat_s(SEQ0, MAXBUF, "\t"); strcat_s(SEQ0, MAXBUF, seq[2]); sep = true; }
				if (READ1 == 1 && CDRH3 == 1) { if (sep) strcat_s(SEQ0, MAXBUF, "\t"); strcat_s(SEQ0, MAXBUF, seq[3]); sep = true; }

				// Check motif's frequency in positives and negatives
				for (int k=0; k<NB1; k++)
				{
					// Read current occurrence
					sscanf_s(buf, "%d", &occ); buf = strchr(buf, '\t') + 1;

					// Update suffix tree
					SetTree0(Tree2, SEQ0, length(SEQ0), 0, NB1, k, occ);
				}
			}

			// Close stream;
			fclose(file6);

			// Verbose
			printf_s("%-15s P(%s) step 3 %10d\n", PRJID, POS[i], nbr3);
		}

		//////////////////////////
		///////// step 6 /////////
		//////////////////////////

		// Get output file for current variable sequences
		file7 = openfile(".\\TEMP\\TMP2.dat", "rt");

		// Get output file for current variable sequences
		file8 = openfile(".\\TEMP\\TMP3.dat", "wt");

		// Verbose
		int nbr4 = 0;

		// Get positive sequences, scroll over all sequences
		while (!feof(file7))
		{
			// Init Cohen's coefficient
			double EFS = 0.0f;

			// init frequency
			double FRQ = 0.0f;

			// Init ranking pvalue
			double PVR = 1.0f;

			// Init frequency and ranking averages
			double AVF1 = 0.0f; double AVR1 = 0.0f;
			double AVF2 = 0.0f; double AVR2 = 0.0f;

			// Init occurrences in positives and negatives
			int OCC1 = 0; int ALL1 = 0; double FRQ1 = 0.0f;
			int OCC2 = 0; int ALL2 = 0; double FRQ2 = 0.0f;

			// Init frequency and ranking standard deviations
			double SDF1 = 0.0f; double SDR1 = 0.0f; double VAR1 = 0.0f;
			double SDF2 = 0.0f; double SDR2 = 0.0f; double VAR2 = 0.0f;

			// counters
			int NB3=0; int NB4=0; int NB5=0; int NB6=0; int NB7=0; int NB8=0;

			///////////////////
			///// step 01 /////
			///////////////////

			// Verbose
			printf_s("%-15s P(%s) step 4 %10d\r", PRJID, POS[i], ++nbr4);

			// Read current line from input file
			if (fgets(BUF, MAXBUF, file7) == NULL) break;

			// Parse buffer
			buf = BUF;

			// Read aa sequences
			if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%s", seq[0], MAXBUF); buf = strchr(buf, '\t') + 1; }
			if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%s", seq[1], MAXBUF); buf = strchr(buf, '\t') + 1; }
			if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%s", seq[2], MAXBUF); buf = strchr(buf, '\t') + 1; }
			if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%s", seq[3], MAXBUF); buf = strchr(buf, '\t') + 1; }

			// Read current occurrence
			sscanf_s(buf, "%d", &OCC1); buf = strchr(buf, '\t') + 1;

			// Read current occurrence
			sscanf_s(buf, "%d", &OCC2); buf = strchr(buf, '\t') + 1;

			// Read current occurrence
			sscanf_s(buf, "%d", &ALL1); buf = strchr(buf, '\t') + 1;

			// Read current occurrence
			sscanf_s(buf, "%d", &ALL2); buf = strchr(buf, '\t') + 1;

			// Read current occurrence
			sscanf_s(buf, "%lf", &FRQ1); buf = strchr(buf, '\t') + 1;

			// Read current occurrence
			sscanf_s(buf, "%lf", &FRQ2); buf = strchr(buf, '\t') + 1;

			// Check fold change frequencies for specificity
			if (FRQ2 > FRQ1 * MINNEG) continue;

			///////////////////
			///// step 02 /////
			///////////////////

			// Init sequence
			SEQ1[0] = '\0';

			// Separator
			sep = false;

			// Print to corresponding output file
			if (READ3 == 1 && CDRL3 == 1) { if (sep) strcat_s(SEQ1, MAXBUF, "\t"); strcat_s(SEQ1, MAXBUF, seq[0]); sep = true; }
			if (READ2 == 1 && CDRH1 == 1) { if (sep) strcat_s(SEQ1, MAXBUF, "\t"); strcat_s(SEQ1, MAXBUF, seq[1]); sep = true; }
			if (READ2 == 1 && CDRH2 == 1) { if (sep) strcat_s(SEQ1, MAXBUF, "\t"); strcat_s(SEQ1, MAXBUF, seq[2]); sep = true; }
			if (READ1 == 1 && CDRH3 == 1) { if (sep) strcat_s(SEQ1, MAXBUF, "\t"); strcat_s(SEQ1, MAXBUF, seq[3]); sep = true; }

			// Get output file
			file9 = openfile(".\\TEMP\\TMP4.dat", "wt");

			// Gett all motifs
			GetTree5(Tree2, SEQ1, length(SEQ1), 0, NB1, SEQ0, &NB3, file9);

			// Close stream
			fclose(file9);

			// Check min occ
			if (NB3 < MINOCC) goto Exit1;

			///////////////////
			///// step 03 /////
			///////////////////

			// Get sequences input filename
			file10 = openfile(".\\TEMP\\TMP4.dat", "rt");

			// Scroll over
			while (!feof(file10))
			{
				// Get current sequence
				if (fgets(BUF, MAXBUF, file10) == NULL) break;

				// Parse buffer
				buf = BUF;

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%s", seq[0], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%s", seq[1], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%s", seq[2], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%s", seq[3], MAXBUF); buf = strchr(buf, '\t') + 1; }

				// Init sequence
				SEQ0[0] = '\0';

				// Separator
				sep = false;

				// Print to corresponding output file
				if (READ3 == 1 && CDRL3 == 1) { if (sep) strcat_s(SEQ0, MAXBUF, "\t"); strcat_s(SEQ0, MAXBUF, seq[0]); sep = true; }
				if (READ2 == 1 && CDRH1 == 1) { if (sep) strcat_s(SEQ0, MAXBUF, "\t"); strcat_s(SEQ0, MAXBUF, seq[1]); sep = true; }
				if (READ2 == 1 && CDRH2 == 1) { if (sep) strcat_s(SEQ0, MAXBUF, "\t"); strcat_s(SEQ0, MAXBUF, seq[2]); sep = true; }
				if (READ1 == 1 && CDRH3 == 1) { if (sep) strcat_s(SEQ0, MAXBUF, "\t"); strcat_s(SEQ0, MAXBUF, seq[3]); sep = true; }

				// Gett all motifs
				GetTree6(Tree1, SEQ0, length(SEQ0), 0);
			}

			// Close stream
			fclose(file10);

			// Get output file
			file11 = openfile(".\\TEMP\\TMP5.dat", "wt");

			// Gett all motifs
			GetTree7(Tree1, length(SEQ1), 0, 2, SEQ0, &NB4, TYP, file11);

			// Close stream
			fclose(file11);

			// Check min occ
			if (NB4 < MINOCC) goto Exit2;

			///////////////////
			///// step 04 /////
			///////////////////

			// Get sequences input filename
			file12 = openfile(".\\TEMP\\TMP4.dat", "rt");

			// Get sequences input filename
			file13 = openfile(".\\TEMP\\TMP6.dat", "wt");

			while (!feof(file12))
			{
				// Read current line
				if (fgets(BUF, MAXBUF, file12) == NULL) break;

				// Parser
				buf = BUF;

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) buf = strchr(buf, '\t') + 1;
				if (READ2 == 1 && CDRH1 == 1) buf = strchr(buf, '\t') + 1;
				if (READ2 == 1 && CDRH2 == 1) buf = strchr(buf, '\t') + 1;
				if (READ1 == 1 && CDRH3 == 1) buf = strchr(buf, '\t') + 1;

				// Init averages in positives and negatives
				double frq1 = +DBL_MAX;
				double frq2 = -DBL_MAX;

				// If negative
				for (int k=0; k < NB1; k++)
				{
					// Read occurrences
					sscanf_s(buf, "%d", &occ); buf = strchr(buf, '\t') + 1;

					switch (TYP[k])
					{
						case +1: frq1 = __MIN(frq1, (double) occ / ALL[k]); break;
						case -1: frq2 = __MAX(frq2, (double) occ / ALL[k]); break;
					}
				}

				// Write frequencies in positives
				fprintf_s(file13, "%e\t+1\n", frq1);

				// Write frequencies in negatives
				fprintf_s(file13, "%e\t-1\n", frq2);

				// Inc nbr of motifs
				NB5 += 1;
			}

			// Close stream
			fclose(file12);
			fclose(file13);

			// Sort sequences by pvalues
			sort(".\\TEMP\\TMP6.dat", ".\\TEMP\\TMP7.dat", "\t", 2, 1, "g", 2, "g");

			// Check nbr of motifs
			if (NB3 != NB5) printf("FATAL ERROR : wrong nbr of motifs !!!\n"), exit(EXIT_FAILURE);

			///////////////////
			///// step 05 /////
			///////////////////

			// Get sequences input filename
			file14 = openfile(".\\TEMP\\TMP7.dat", "rt");

			while (!feof(file14))
			{
				// Parse fileds
				if (fscanf_s(file14, "%lf\t%d", &FRQ, &occ) != 2) break;

				// Inc nbr of motifs
				NB6 += 1;

				// Update average
				switch (occ)
				{
					case +1: AVF1 += FRQ; AVR1 += NB6; break;
					case -1: AVF2 += FRQ; AVR2 += NB6; break;
				}
			}

			// Close stream
			fclose(file14);

			// Get frequency average
			AVF1 = NB5 > 0 ? (AVF1 / NB5) : 0.0f;
			AVF2 = NB5 > 0 ? (AVF2 / NB5) : 0.0f;

			// Get rank average
			AVR1 = NB5 > 0 ? (AVR1 / NB5) : 0.0f;
			AVR2 = NB5 > 0 ? (AVR2 / NB5) : 0.0f;

			// Check nbr of motifs
			if (NB5 != NB6 / 2) printf("FATAL ERROR : wrong nbr of motifs !!!\n"), exit(EXIT_FAILURE);

			///////////////////
			///// step 06 /////
			///////////////////

			// Get sequences input filename
			file15 = openfile(".\\TEMP\\TMP7.dat", "rt");

			while (!feof(file15))
			{
				// Parse fileds
				if (fscanf_s(file15, "%lf\t%d", &FRQ, &occ) != 2) break;

				// Inc nbr of motifs
				NB7 += 1;

				// Update average
				switch (occ)
				{
					case +1: SDF1 += pow(AVF1 - FRQ, 2) / NB5; SDR1 += pow(AVR1 - (double)NB7, 2) / NB5; VAR1 += pow(AVF1 - FRQ, 2) / ((double)NB5 - 1); break;
					case -1: SDF2 += pow(AVF2 - FRQ, 2) / NB5; SDR2 += pow(AVR2 - (double)NB7, 2) / NB5; VAR2 += pow(AVF2 - FRQ, 2) / ((double)NB5 - 1); break;
				}
			}

			// Close stream
			fclose(file15);

			// Get frequency standard deviation
			SDF1 = NB5 > 0 ? sqrt(SDF1) : 0.0f;
			SDF2 = NB5 > 0 ? sqrt(SDF2) : 0.0f;

			// Get rank standard deviation
			SDR1 = NB5 > 0 ? sqrt(SDR1) : 0.0f;
			SDR2 = NB5 > 0 ? sqrt(SDR2) : 0.0f;

			// Check nbr of motifs
			if (NB6 != NB7) printf("FATAL ERROR : wrong nbr of motifs !!!\n"), exit(EXIT_FAILURE);

			///////////////////
			///// step 07 /////
			///////////////////

			// Calculate pvalue
			PVR = pvalue(AVR1, AVR2, SDR1, SDR2, NB5, NB5);

			// Caulculate Cohen's coefficient
			EFS = Effect(AVF1, AVF2, VAR1, VAR2, NB5, NB5);

			// Check pvalue significance for specificity
			if (PVR > CUTOFF || EFS < COHIDX) goto Exit3;

			///////////////////
			///// step 08 /////
			///////////////////

			// Init positional frequencies
			int   NBR = 0;
			int** OCC = Alloc2D(AASYMB + 1, length(SEQ1), 0);

			// Get sequences input filename
			file16 = openfile(".\\TEMP\\TMP5.dat", "rt");

			// Look for residues in current sequence
			while (!feof(file16))
			{
				// Get current sequence
				if (fgets(BUF, MAXBUF, file16) == NULL) break;

				// Parse buffer
				buf = BUF;

				// Read aa sequences
				if (READ3 == 1 && CDRL3 == 1) { sscanf_s(buf, "%s", seq[0], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH1 == 1) { sscanf_s(buf, "%s", seq[1], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ2 == 1 && CDRH2 == 1) { sscanf_s(buf, "%s", seq[2], MAXBUF); buf = strchr(buf, '\t') + 1; }
				if (READ1 == 1 && CDRH3 == 1) { sscanf_s(buf, "%s", seq[3], MAXBUF); buf = strchr(buf, '\t') + 1; }

				// Init sequence
				SEQ2[0] = '\0';

				// Separator
				sep = false;

				// Print to corresponding output file
				if (READ3 == 1 && CDRL3 == 1) { if (sep) strcat_s(SEQ2, MAXBUF, "\t"); strcat_s(SEQ2, MAXBUF, seq[0]); sep = true; }
				if (READ2 == 1 && CDRH1 == 1) { if (sep) strcat_s(SEQ2, MAXBUF, "\t"); strcat_s(SEQ2, MAXBUF, seq[1]); sep = true; }
				if (READ2 == 1 && CDRH2 == 1) { if (sep) strcat_s(SEQ2, MAXBUF, "\t"); strcat_s(SEQ2, MAXBUF, seq[2]); sep = true; }
				if (READ1 == 1 && CDRH3 == 1) { if (sep) strcat_s(SEQ2, MAXBUF, "\t"); strcat_s(SEQ2, MAXBUF, seq[3]); sep = true; }

				// Read current occurrence
				sscanf_s(buf, "%d", &occ); buf = strchr(buf, '\t') + 1;

				// Update aa frequencies
				for (int l=0; l < length(SEQ2); l++)
				{
					// Update aa frequency
					if (SEQ2[l] == '\t') continue;

					// Update aa frequency
					OCC[aa2id(SEQ2[l])][l] += occ;
				}

				// Update total occurrences
				NBR += occ;

				// Inc nbr of sequences
				NB8 += 1;
			}

			// Close stream
			fclose(file16);

			// Check nbr of motifs
			if (NB4 != NB8) printf("FATAL ERROR : wrong nbr of motifs !!!\n"), exit(EXIT_FAILURE);

			///////////////////
			///// step 09 /////
			///////////////////

			// Write family sequence
			for (int k=0, l=0; k < length(SEQ1); k++)
			{
				// Check if another sequence
				if (SEQ1[k] == '\t')
				{
					// Separator
					SEQ3[l++] = '\t', SEQ3[l] = '\0'; continue;
				}

				// Init current frequencies
				int    a1 = 0;
				int    a2 = 0;
				float* AA = Alloc1D(AASYMB, 0.0f);

				// Get residues frequencies
				for (int m=0; m < AASYMB; m++)
				{
					AA[m] += (float) OCC[m][k] / NBR;
				}

				for (int m=0; m < AASYMB; m++)
				{
					a1 += AA[m] >= 1.0f - MINFRQ;
					a2 += AA[m] >= MINFRQ - 0.0f;
				}

				if (a1 == 1)
				{
					for (int m=0; m < AASYMB; m++)
					{
						if (AA[m] >= 1.0f - MINFRQ) SEQ3[l++] = id2aa(m), SEQ3[l] = '\0';
					}
				}
				else
				if (a2 == 0)
				{
					SEQ3[l++] = '.', SEQ3[l] = '\0';
				}
				else
				{
					// Open bracket
					SEQ3[l++] = '[', SEQ3[l] = '\0';

					for (int m=0; m < AASYMB; m++)
					{
						if (AA[m] >= MINFRQ) SEQ3[l++] = id2aa(m), SEQ3[l] = '\0';
					}

					// Close bracket
					SEQ3[l++] = ']', SEQ3[l] = '\0';
				}

				// free memory
				delete[] AA;
			}

			// free memory
			for (int k=0; k < AASYMB + 1; k++)
			delete[] OCC[k];
			delete[] OCC;

			///////////////////
			///// step 10 /////
			///////////////////

			// current candidate
			fprintf_s(file8, "%s\t%d\t%d\t%d\t%s\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e", SEQ1, OCC1, OCC2, NB4, SEQ3, AVF1, SDF1, AVF2, SDF2, AVF2 > 0.0f ? AVF1/AVF2 : 0.0f);

			// write specificity pvalue
			fprintf_s(file8, "\t%.4e\n", PVR);

			// delete temp files
			system("del /F /Q .\\TEMP\\TMP4.dat > nul");
			system("del /F /Q .\\TEMP\\TMP5.dat > nul");
			system("del /F /Q .\\TEMP\\TMP6.dat > nul");
			system("del /F /Q .\\TEMP\\TMP7.dat > nul");

			// next
			continue;

			///////////////////
			///// step 11 /////
			///////////////////

		Exit1:	{
				// delete temp files
				system("del /F /Q .\\TEMP\\TMP4.dat > nul");

				// next
				continue;
			}

			///////////////////
			///// step 12 /////
			///////////////////

		Exit2:	{
				// delete temp files
				system("del /F /Q .\\TEMP\\TMP4.dat > nul");
				system("del /F /Q .\\TEMP\\TMP5.dat > nul");

				// next
				continue;
			}

			///////////////////
			///// step 13 /////
			///////////////////

		Exit3:	{
				// delete temp files
				system("del /F /Q .\\TEMP\\TMP4.dat > nul");
				system("del /F /Q .\\TEMP\\TMP5.dat > nul");
				system("del /F /Q .\\TEMP\\TMP6.dat > nul");
				system("del /F /Q .\\TEMP\\TMP7.dat > nul");

				// next
				continue;
			}
		}

		// Close stream
		fclose(file7);
		fclose(file8);

		// free Tree
		EmptyTree(Tree1);
		EmptyTree(Tree2);

		// free memory
		for (int j=0; j<NBREAD; j++)
		delete[] seq[j];
		delete[] seq;

		// free memory
		delete[] TYP;
		delete[] ALL;

		//////////////////////////
		///////// step 7 /////////
		//////////////////////////

		// Sort sequences by pvalues
		sort(".\\TEMP\\TMP3.dat", ".\\TEMP\\TMP8.dat", "\t", 4, 9 + 2 * (READ1 * CDRH3 + READ2 * CDRH1 + READ2 * CDRH2 + READ3 * CDRL3), "g",
															   10 + 2 * (READ1 * CDRH3 + READ2 * CDRH1 + READ2 * CDRH2 + READ3 * CDRL3), "g",
															    1 + 1 * (READ1 * CDRH3 + READ2 * CDRH1 + READ2 * CDRH2 + READ3 * CDRL3), "gr",
															    3 + 1 * (READ1 * CDRH3 + READ2 * CDRH1 + READ2 * CDRH2 + READ3 * CDRL3), "gr"); break;

		// Get output file for current variable sequences
		file17 = openfile(".\\TEMP\\TMP8.dat", "rt");

		// Read input file
		file18 = openfile(PATH1, "STEP04", COM, "csv", "wt");

		// Print to corresponding output file
		if (READ3 == 1 && CDRL3 == 1) fprintf_s(file18, "#L3 AA SEQ,");
		if (READ2 == 1 && CDRH1 == 1) fprintf_s(file18, "#H1 AA SEQ,");
		if (READ2 == 1 && CDRH2 == 1) fprintf_s(file18, "#H2 AA SEQ,");
		if (READ1 == 1 && CDRH3 == 1) fprintf_s(file18, "#H3 AA SEQ,");

		// Write headers
		fprintf_s(file18, "#OBS IN POS (MIN),#OBS IN NEG (MAX),#HOM IN POS (MIN),");

		// Print to corresponding output file
		if (READ3 == 1 && CDRL3 == 1) fprintf_s(file18, "#L3 AA MOT,");
		if (READ2 == 1 && CDRH1 == 1) fprintf_s(file18, "#H1 AA MOT,");
		if (READ2 == 1 && CDRH2 == 1) fprintf_s(file18, "#H2 AA MOT,");
		if (READ1 == 1 && CDRH3 == 1) fprintf_s(file18, "#H3 AA MOT,");

		// Write headers
		fprintf_s(file18, "#FREQ IN POS (AVR),#FREQ IN POS (STD),#FREQ IN NEG (AVR),#FREQ IN NEG (STD),#FREQ (POS/NEG),#PVAL\n");

		// Scroll over
		while (!feof(file17))
		{
			// Read current line
			if (fgets(BUF, MAXBUF, file17) == NULL) break;

			// trim buffer
			sscanf_s(BUF, "%[^\n]\n", BUF, MAXBUF);

			// Convert tab to csv
			for (int k=0; k < length(BUF); k++) { if (BUF[k] == '\t') BUF[k] = ','; }

			// wtite sequence
			fprintf(file18, "%s\n", BUF);
		}

		// Close streams
		fclose(file17);
		fclose(file18);

		//////////////////////////
		///////// step 9 /////////
		//////////////////////////

		// Delete temporary files
		system("del /F /Q .\\TEMP\\TMP*.* > nul");

		// Verbose
		printf_s("%-15s P(%s) step 4 %10d\n", PRJID, POS[i], nbr4);

		// Verbose
		printf_s("\n");
	}

	////////////////////////////////////////////////////////////////
	///////////////////////////// step 4 ///////////////////////////
	////////////////////////////////////////////////////////////////

	// free memory
	delete[] for1;
	delete[] for2;
	delete[] rev1;
	delete[] rev2;
	delete[] samp;
	delete[] libr;
	delete[] cont;

	// free memory
	for (int i=0; i<NB1; i++)
	{
		delete[] FOR1[i];
		delete[] FOR2[i];
		delete[] REV1[i];
		delete[] REV2[i];
		delete[] SAMP[i];
		delete[] LIBR[i];
		delete[] CONT[i];
	}

	// free memory
	delete[] FOR1;
	delete[] FOR2;
	delete[] REV1;
	delete[] REV2;
	delete[] SAMP;
	delete[] LIBR;
	delete[] CONT;

	// free memory
	for (int i=0; i<NB2; i++)
	{
		delete[] POS[i];
		delete[] NEG[i];
		delete[] LIB[i];
		delete[] LAB[i];
	}

	// free memory
	delete[] POS; delete[] pos;
	delete[] NEG; delete[] neg;
	delete[] LIB; delete[] lib;
	delete[] LAB; delete[] lab;
}

//---------------------------------------------------------------------------
// CellectSeq pipeline
void pipeline(char *project)
{
	// STEP 0: Parse input arguments
	ParseAgrs(project);

	// STEP 1: Sequence frequencies
	Frequence();

	// STEP 2: Find motifs in positives
	Positives();

	// STEP 3: Predict specific sequences
	Predicted();
}

//---------------------------------------------------------------------------
// main routine
void main()
{
	pipeline("01 CD151");
	pipeline("02 CA9");
	pipeline("03 ITGA11");
}
