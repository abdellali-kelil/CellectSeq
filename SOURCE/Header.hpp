//---------------------------------------------------------------------------
// stack size
#pragma comment (linker, "/STACK:0x500000")

//---------------------------------------------------------------------------
// headers
#include <Windows.h>
#include <process.h>
#include <limits.h>
#include <stdio.h>
#include <io.h>

// constants
constexpr auto MAXBUF = 1024;		// Max buffer length
constexpr auto AASYMB =   20;		// Number of amino acids
constexpr auto NASYMB =    4;		// Number of nucleic acids
constexpr auto NBREAD =    4;		// Max nbr of CDR reads
constexpr auto WCSYMB =   'X';		// Symbol sequencing error

int    NBRTHR = 0;						// Number of simultaneous threads
int    MINOCC = 0;						// Min occurrences of sequence in positives
float  MINNEG = 0.0f;					// Min frequency of motifs in negatives over positives
float  MAXNWC = 0.0f;					// Max wildcards frequency
float  MINFRQ = 0.0f;					// Min frequency of amino acids at each position
double COHIDX = 0.0f;					// Min size effect for p-value significance
double CUTOFF = 0.0f;					// pvalue cutoff

// Reads considered
char READ1 = 0;
char READ2 = 0;
char READ3 = 0;

// CDRs considered
char CDRL3 = 0;
char CDRH1 = 0;
char CDRH2 = 0;
char CDRH3 = 0;

// Working files/folders
char *PRJID = NULL;					// Project Folder (antigen)
char *PATH1 = NULL;					// Results Folder (antigen)
char *DEMUX = NULL;					// Demux file
char *TESTS = NULL;					// Tests file

//---------------------------------------------------------------------------
// Structur : Suffixe Tree
typedef struct CTree
{
	int    IDX;
	char   RES;
	int   *OCC;
	CTree *BRO;
	CTree *SON;
} CTree;

//---------------------------------------------------------------------------
// Add new to tree
CTree* AddBRO(CTree *Tree)
{
	// get first
	CTree *ptr = Tree;

	// get last
	while (ptr->BRO != NULL) { ptr = ptr->BRO; }

	// add new
	ptr->BRO = new CTree;

	// get new
	ptr = ptr->BRO;
	
	// set no next
	ptr->BRO = NULL;
	ptr->SON = NULL;

	// send new
	return ptr;
}

//---------------------------------------------------------------------------
// Add new to tree
CTree* AddSON(CTree *Tree)
{
	// get first
	CTree *ptr = Tree;

	// get last
	while (ptr->SON != NULL) { ptr = ptr->SON; }

	// add new
	ptr->SON = new CTree;

	// get new
	ptr = ptr->SON;
	
	// set no next
	ptr->BRO = NULL;
	ptr->SON = NULL;

	// send new
	return ptr;
}

//---------------------------------------------------------------------------
// Structur : Suffixe Tree
typedef struct CThread
{
	int    THR;		// Current thread
	int    CUR;		// Current length of motifs
	int    DIM;		// Dimension, number of pools
	int    NB1;		// Total number of sequences in positive pool
	int    NB2;		// Total number of sequences in negative pool
	int   *ALL;		// Total nbr of sequences (occurrences) in each pool
	int   *CDR;		// Lengths of CDR sequences
	int   *TYP;		// Type of each dimention, positive/negative
	int  **OCC2;	// Total Number of sequences in each pool
	char **SEQ1;	// Sequences in positive pool
	char **SEQ2;	// sequences in negative pool
	CTree *TREE;	// Suffix tree
} CThread;

//---------------------------------------------------------------------------
// Swap two variables
template<class typ> inline void SWAP(typ &a, typ &b) { typ t = a; a = b; b = t; }

//---------------------------------------------------------------------------
// MAX 2 values
template<class typ> inline typ __MAX(typ a, typ b) { return a > b ? a : b; }

//---------------------------------------------------------------------------
// MAX 3 values
template<class typ> inline typ __MAX(typ a, typ b, typ c) { return __MAX(a,b) > c ? __MAX(a,b) : c; }

//---------------------------------------------------------------------------
// MIN 2 values
template<class typ> inline typ __MIN(typ a, typ b) { return a < b ? a : b; }

//---------------------------------------------------------------------------
// MIN 3 values
template<class typ> inline typ __MIN(typ a, typ b, typ c) { return __MIN(a,b) < c ? __MIN(a,b) : c; }

//---------------------------------------------------------------------------
// Resize memory allocation
template<class typ> inline typ* Realloc(typ *src, int size)
{
	// alloc memory
	typ *dst = new typ[size];

	// check if null
	if (src == NULL) return dst;

	// copy source
	errno_t err = memcpy_s(dst, size * sizeof(typ), src, size * sizeof(typ));

	// check error
	if (err != 0) { printf_s("MEMORY REALLOCATION EEROR !!!\n"); exit(EXIT_FAILURE); }

	// free mem
	delete[] src;

	// return mem
	return dst;
}

//---------------------------------------------------------------------------
// 1 dimension(s) memory allocation
template<class typ> inline typ* Alloc1D(int size, typ value)
{
	// Allocate memory
	typ *ptr = new typ[size];

	// Init values
	for (int i = 0; i<size; i++) ptr[i] = value;

	// Return new memory allocation
	return ptr;
}

//---------------------------------------------------------------------------
// 2 dimension(s) allocation
template<class typ> inline typ** Alloc2D(int size1, int size2, typ value)
{
	// Allocate memory
	typ **ptr = new typ*[size1];

	// Init values
	for (int i = 0; i<size1; i++)
	{
		// Allocate memory
		ptr[i] = new typ[size2];

		// init values
		for (int j = 0; j<size2; j++) ptr[i][j] = value;
	}

	// Return new memory allocation
	return ptr;
}

//---------------------------------------------------------------------------
// 3 dimension(s) memory allocation
template<class typ> inline typ*** Alloc3D(int size1, int size2, int size3, typ value)
{
	// Allocate memory
	typ ***ptr = new typ**[size1];

	// Init values
	for (int i = 0; i<size1; i++)
	{
		// Allocate memory
		ptr[i] = new typ*[size2];

		for (int j = 0; j<size2; j++)
		{
			// Allocate memory
			ptr[i][j] = new typ[size3];

			// set values
			for (int k = 0; k<size3; k++) ptr[i][j][k] = value;
		}
	}

	// Return new memory allocation
	return ptr;
}

//---------------------------------------------------------------------------
// Put stop char at token positions
inline void strtoks(char *str, const char *tok)
{
	char *cur = NULL;
	char *nex = NULL;

	// get first token:
	cur = strtok_s(str, tok, &nex);

	// get next tokens
	while (cur != NULL) { cur = strtok_s(NULL, tok, &nex); }
}

//---------------------------------------------------------------------------
// Calculate length of char*
inline int length(char *str)
{
	// Init length
	int len = 0; char *ptr = str;

	// Scroll over
	while (*ptr != '\0' && *ptr != '\n' && *ptr != '\r') { ptr++; len++; }

	// Return length
	return len;
}

// ---------------------------------------------------------------------------
// Open file stream
FILE *openfile(char *path1, char *path2, char *file, char *ext, const char *p)
{
	// Buffer
	char BUF[MAXBUF] = "\0";

	if (length(file) > 200)
	{
		// Set to long filename format
		sprintf_s(BUF, MAXBUF, "\\\\?\\%s\\%s\\%s.%s", path1, path2, file, ext);
	}
	else
	{
		// Set to short filename format
		sprintf_s(BUF, MAXBUF, "%s\\%s\\%s.%s", path1, path2, file, ext);
	}

	// Open for read
	FILE *stream = _fsopen(BUF, p, _SH_DENYNO);

	// Check stream error
	if (stream == NULL) { printf_s("\nFile \"%s\" not found !!!\n", file); exit(EXIT_FAILURE); }

	// Return opened file
	return stream;
}

// ---------------------------------------------------------------------------
// Open file stream
FILE *openfile(char *file, const char *p)
{
	// Buffer
	char BUF[MAXBUF] = "\0";

	if (length(file) > 200)
	{
		// Set to long filename format
		sprintf_s(BUF, MAXBUF, "\\\\?\\%s", file);
	}
	else
	{
		// Set to short filename format
		sprintf_s(BUF, MAXBUF, "%s", file);
	}

	// Open for read
	FILE *stream = _fsopen(BUF, p, _SH_DENYNO);

	// Check stream error
	if (stream == NULL) { printf_s("\nFile \"%s\" not found !!!\n", file); exit(EXIT_FAILURE); }

	// Return opened file
	return stream;
}

//---------------------------------------------------------------------------
// Mapt amino acid one letter symbol to id
inline char aa2id(char aa)
{
	// wildcard symbol
	if (aa == WCSYMB) return 20;

	// other symbols
	switch (aa)
	{
		// standard aa.
		case 'A': return  0;
		case 'C': return  1;
		case 'D': return  2;
		case 'E': return  3;
		case 'F': return  4;
		case 'G': return  5;
		case 'H': return  6;
		case 'I': return  7;
		case 'K': return  8;
		case 'L': return  9;
		case 'M': return 10;
		case 'N': return 11;
		case 'P': return 12;
		case 'Q': return 13;
		case 'R': return 14;
		case 'S': return 15;
		case 'T': return 16;
		case 'V': return 17;
		case 'W': return 18;
		case 'Y': return 19;
		case 'X': return 20;
		case '*': return 20;
		case '.': return 20;
		case '_': return 20;

		// Invalid residue
		default : printf_s("FATAL ERROR : wrong AA !!!\n"); exit(EXIT_FAILURE);
	}
}

//---------------------------------------------------------------------------
// Map amino acid id to one letter symbol
inline char id2aa(char id)
{
	switch (id)
	{
		// standard aa.
		case  0: return 'A';
		case  1: return 'C';
		case  2: return 'D';
		case  3: return 'E';
		case  4: return 'F';
		case  5: return 'G';
		case  6: return 'H';
		case  7: return 'I';
		case  8: return 'K';
		case  9: return 'L';
		case 10: return 'M';
		case 11: return 'N';
		case 12: return 'P';
		case 13: return 'Q';
		case 14: return 'R';
		case 15: return 'S';
		case 16: return 'T';
		case 17: return 'V';
		case 18: return 'W';
		case 19: return 'Y';
		case 20: return '.';

		// Invalid residue
		default : printf_s("FATAL ERROR : wrong AA !!!\n"); exit(EXIT_FAILURE);
	}
}

// ---------------------------------------------------------------------------
// Get selection pool
void Selection(char *sel, char *str1, char *str2, char *str3)
{
	// format id
	if (_stricmp(str3, "NA") == 0)
	{
		sprintf_s(sel, MAXBUF, "%s_%s", str1, str2);
	}
	else
	{
		sprintf_s(sel, MAXBUF, "%s_%s(%s)", str1, str2, str3);
	}
}

// ---------------------------------------------------------------------------
// Get samples combination id
void Combination(char *com,char *pos, char *neg)
{
	// format id
	sprintf_s(com, MAXBUF, "P(%s)N(%s)", pos, neg);
}

// ---------------------------------------------------------------------------
// Discover motifs in positives
void Library(char *lib)
{
	if (_stricmp(lib, "LibF") == 0)
	{
		// CDR sequences
		CDRL3 = 1;
		CDRH1 = 1;
		CDRH2 = 1;
		CDRH3 = 1;
	}
}

//---------------------------------------------------------------------------
// file multiple keys sorting
void sort(char *inp, char*out, char *arg, int nbr, ...)
{
	// buffer
	char BUF[MAXBUF] = "\0";

	// Init buffer
	sprintf_s(BUF, MAXBUF, "sort -T .\\TEMP\\ --parallel=40 -t\"%s\"", arg);

	// init args
	va_list ap; va_start(ap, 2*nbr);

	// read args
	for(int i=0; i<nbr; i++) 
	{
		char key[MAXBUF] = "\0";
		char typ[MAXBUF] = "\0";

		// sort key
		sprintf_s(key, "%d", va_arg(ap, int));

		// sort jey type
		sprintf_s(typ, "%s", va_arg(ap, char*));

		strcat_s(BUF, MAXBUF, " ");
		strcat_s(BUF, MAXBUF, "-k");
		strcat_s(BUF, MAXBUF, key);
		strcat_s(BUF, MAXBUF, ",");
		strcat_s(BUF, MAXBUF, key);
		strcat_s(BUF, MAXBUF, typ);
	}

	// end args
	va_end(ap);

	// Input file
	strcat_s(BUF, MAXBUF, " ");
	strcat_s(BUF, MAXBUF, inp);

	// separator
	strcat_s(BUF, MAXBUF, " ");
	strcat_s(BUF, MAXBUF, ">");

	// Input file
	strcat_s(BUF, MAXBUF, " ");
	strcat_s(BUF, MAXBUF, out);

	// run sort
	system(BUF);
}

// ---------------------------------------------------------------------------
// Add sequence to suffix tree
void SetTree0(CTree *Tree, char *SEQ, int LEN, int POS, int DIM, int CUR, int OCC)
{
	// init match
	bool found = false;

	// get first
	CTree *SON = Tree->SON;

	// Scroll over all
	while (SON != NULL)
	{
		// check if match
		if (SON->RES == SEQ[POS])
		{
			// Check max length
			if (POS == LEN - 1)
			{
				// Check suffix
				if (SON->IDX < 0)
				{
					// Set index
					SON->IDX = +1;

					// Check dimension
					if (DIM > 0)
					{
						// Init occurrences
						SON->OCC = Alloc1D(DIM, 0);
					}
				}

				// Check dimension
				if (DIM > 0)
				{
					// Update occurrences
					SON->OCC[CUR] += OCC;
				}
			}
			else
			{
				// Next
				SetTree0(SON, SEQ, LEN, POS + 1, DIM, CUR, OCC);
			}

			// set match and exit
			found = true; break;
		}

		// next
		SON = SON->BRO;
	}

	// check if match
	if (found == false)
	{
		// init new
		CTree *ADD;

		// check if first new
		if (Tree->SON == NULL)
		{
			// add new
			ADD = AddSON(Tree);
		}
		else
		{
			// add next
			ADD = AddBRO(Tree->SON);
		}

		// New
		ADD->IDX = -1;
		ADD->RES = SEQ[POS];
		ADD->OCC = NULL;
		ADD->BRO = NULL;
		ADD->SON = NULL;

		// Check max length
		if (POS == LEN - 1)
		{
			// Set index
			ADD->IDX = +1;

			// Check dimention
			if (DIM > 0)
			{
				ADD->OCC      = Alloc1D(DIM, 0);
				ADD->OCC[CUR] = OCC;
			}
		}
		else
		{
			// Next
			SetTree0(ADD, SEQ, LEN, POS + 1, DIM, CUR, OCC);
		}
	}
}

// ---------------------------------------------------------------------------
// Add motif to suffix tree
void SetTree1(CTree *Tree, char *MOT, int LEN, int POS, int DIM, int NWC, int CUR, int *CDR)
{
	// get first
	CTree *SON = Tree->SON;

	// Keep track of matche
	char resid;
	bool find1 = false;
	bool find2 = false;

	// Scroll over all
	while (SON != NULL)
	{
		// Check if current residue matches with tree suffix
		if (SON->RES == MOT[POS] || SON->RES == '.')
		{
			// Get current residue
			resid = SON->RES;

			// Match found
			find1 = find1 || resid == '.';
			find2 = find2 || resid == MOT[POS];

			// Check for end and occurrences in positives
			if (POS < LEN - 1) SetTree1(SON, MOT, LEN, POS + 1, DIM, (resid == '_') ? 0 : NWC + (resid == '.'), (resid == '_') ? (CUR + 1) : CUR, CDR);
		}

		// next
		SON = SON->BRO;
	}

	// Check if we need to add wildcard
	if (find1 == false && POS == LEN - 1)
	{
		// Check wildcard conditions
		if (NWC + 1 <= MAXNWC * CDR[CUR] && MOT[POS] != '_')
		{
			// init new
			CTree *ADD;

			// check if first new
			if (Tree->SON == NULL)
			{
				// add new
				ADD = AddSON(Tree);
			}
			else
			{
				// add next
				ADD = AddBRO(Tree->SON);
			}

			// New
			ADD->IDX =  1;
			ADD->RES = '.';
			ADD->OCC = Alloc1D(DIM, 0);
			ADD->BRO = NULL;
			ADD->SON = NULL;
		}
	}

	// Check if we need to add residue
	if (find2 == false && POS == LEN - 1)
	{
		// Check wildcards
		if (NWC <= MAXNWC * CDR[CUR] || MOT[POS] == '_')
		{
			// init new
			CTree *ADD;

			// check if first new
			if (Tree->SON == NULL)
			{
				// add new
				ADD = AddSON(Tree);
			}
			else
			{
				// add next
				ADD = AddBRO(Tree->SON);
			}

			// New
			ADD->IDX =  1;
			ADD->RES = MOT[POS];
			ADD->OCC = Alloc1D(DIM, 0);
			ADD->BRO = NULL;
			ADD->SON = NULL;
		}
	}
}

// ---------------------------------------------------------------------------
// Add motif to suffix tree
void SetTree2(CTree *Tree, char *SEQ, int LEN, int POS, int DIM, int CUR, int OCC)
{
	// get first
	CTree *SON = Tree->SON;

	// Init no match
	bool found = false;

	// Scroll over all
	while (SON != NULL)
	{
		// Check if residues match
		if (SON->RES == SEQ[POS])
		{
			// Check motif's end
			SON->OCC[CUR] += OCC;

			// Check for end
			if (POS < LEN - 1) SetTree2(SON, SEQ, LEN, POS + 1, DIM, CUR, OCC);

			// Match found, exit
			found = true; break;
		}

		// next
		SON = SON->BRO;
	}

	// Check if match found
	if (found == false)
	{
		// init new
		CTree *ADD;

		// check if first new
		if (Tree->SON == NULL)
		{
			// add new
			ADD = AddSON(Tree);
		}
		else
		{
			// add next
			ADD = AddBRO(Tree->SON);
		}

		// New
		ADD->IDX = -1;
		ADD->RES = SEQ[POS];
		ADD->OCC = NULL;
		ADD->BRO = NULL;
		ADD->SON = NULL;

		// Check last aa
		ADD->OCC      = Alloc1D(DIM, 0);
		ADD->OCC[CUR] = OCC;

		// Check for end
		if (POS < LEN - 1) SetTree2(ADD, SEQ, LEN, POS + 1, DIM, CUR, OCC);
	}
}

// ---------------------------------------------------------------------------
// Get a motif from Suffix Tree
void GetTree0(CTree *Tree, char *SEQ, int POS, int DIM, int *NBR, FILE *file)
{
	// get first son
	CTree *SON = Tree->SON;

	// Scroll over sons
	while(SON != NULL)
	{
		// Add current residue to motif
		SEQ[POS + 0] = SON->RES;
		SEQ[POS + 1] = '\0';

		// Write if end of sequence
		if (SON->IDX > 0)
		{
			// Inc occurrences
			if (NBR != NULL) *NBR += 1;

			// Write sequences
			fprintf_s(file, "%s", SEQ);

			// Write sequence occurrences
			for (int j=0; j<DIM; j++) fprintf_s(file, "\t%d", SON->OCC[j]); fprintf_s(file, "\n");
		}

		// Next branch
		GetTree0(SON, SEQ, POS + 1, DIM, NBR, file);

		// Next
		SON = SON->BRO;
	}
}

// ---------------------------------------------------------------------------
// Add motif to suffix tree
void GetTree1(CTree *Tree, char *MOT, int LEN, int POS, int CUR, int OCC)
{
	// get first son
	CTree *SON = Tree->SON;

	// Scroll over sons
	while(SON != NULL)
	{
		// Check if current residue matches with tree suffix
		if ((SON->RES == MOT[POS]) || (SON->RES == '.' && MOT[POS] != '\t' && MOT[POS] != '_'))
		{
			// Search next
			if (POS == LEN - 1)
			{
				if (SON->IDX > 0) SON->OCC[CUR] += OCC;
			}
			else
			{
				GetTree1(SON, MOT, LEN, POS + 1, CUR, OCC);
			}
		}

		// Next
		SON = SON->BRO;
	}
}

// ---------------------------------------------------------------------------
// Add motif to suffix tree
void GetTree2(CTree *Tree, char *MOT, int LEN, int POS, int DIM, int *OCC, int *TYP)
{
	// get first son
	CTree *SON = Tree->SON;

	// Scroll over sons
	while(SON != NULL)
	{
		// Check suffix order 
		if (aa2id(SON->RES) < aa2id(MOT[POS])) return;

		// Check if current residue matches with tree suffix
		if (SON->RES == MOT[POS] || SON->RES == '.')
		{
			// Search next
			if (POS == LEN - 1)
			{
				// Check if leaf
				if (SON->IDX > 0)
				{
					for (int j=0; j<DIM; j++) if (TYP[j] != 0) SON->OCC[j] += OCC[j];
				}
			}
			else
			{
				// Next
				GetTree2(SON, MOT, LEN, POS + 1, DIM, OCC, TYP);
			}
		}

		// Next
		SON = SON->BRO;
	}
}

// ---------------------------------------------------------------------------
// Get a motif from Suffix Tree
void GetTree3(CTree *Tree, char *SEQ, int LEN, int POS, int DIM, int *TYP, int *NBR, FILE *file)
{
	// get first son
	CTree *SON1 = Tree->SON;

	// Scroll over sons
	while (SON1 != NULL)
	{
		// Get current residue
		SEQ[POS + 0] = SON1->RES;
		SEQ[POS + 1] = '\0';

		// Write if end of sequence
		if (POS < LEN - 1) goto Next1;

		// Check wildcard
		if (SON1->RES != '.') goto Next2;

		// Init nbr of aa
		int NAA = 0;

		CTree *SON2 = Tree->SON;

		// Scroll over sons
		while (SON2 != NULL)
		{
			// Skeep if same
			if (aa2id(SON2->RES) < AASYMB)
			{
				// Init occurrences
				int occ1 = +INT_MAX;
				int occ2 = -INT_MAX;

				// Find aa min occurrences in pos
				for (int k=0; k<DIM; k++)
				{
					switch (TYP[k])
					{
						case +1: occ1 = __MIN(occ1, SON2->OCC[k]); break;
						case -1: occ2 = __MAX(occ2, SON2->OCC[k]); break;
					}
				}

				// Update frequent aa
				NAA += occ1 >= MINOCC || occ2 >= MINOCC;
			}

			// Next
			SON2 = SON2->BRO;
		}

		// Skeep if not wildcard
		if (NAA < 2) goto Next0;

	Next2 :	{
			// inc motifs
			*NBR += 1;

			// Write sequence and occurrences
			fprintf_s(file, "%s", SEQ); for (int j=0; j<DIM; j++) fprintf_s(file, "\t%d", (int)SON1->OCC[j]); fprintf_s(file, "\n");
		}

	Next1 :	{
			// Next branch
			GetTree3(SON1, SEQ, LEN, POS + 1, DIM, TYP, NBR, file);
		}

	Next0:	{
			// Next
			SON1 = SON1->BRO;
		}
	}
}

// ---------------------------------------------------------------------------
// Get a motif from Suffix Tree
void GetTree4(CTree *Tree, char *MOT, int LEN, int POS, int DIM, int *TYP, char *SEQ, int ***OCC)
{
	// get first son
	CTree *SON = Tree->SON;

	// Scroll over sons
	while (SON != NULL)
	{
		// Check suffix order
		if (aa2id(SON->RES) > aa2id(MOT[POS])) return;

		// Check if current residue matches with tree suffix
		if (SON->RES == MOT[POS] || MOT[POS] == '.')
		{
			SEQ[POS + 0] = SON->RES;
			SEQ[POS + 1] = '\0';

			// Check end of sequence
			if (POS == LEN - 1)
			{
				for (int j=0; j<LEN; j++)
				{
					// Skeep if not wildcard
					if (MOT[j] == '.')
					{
						// Scroll over dim
						for (int k=0; k<DIM; k++)
						{
							// Skeep if not positive
							if (TYP[k] != 1) continue;

							// Update each amino acid frequency
							OCC[j][k][aa2id(SEQ[j])] += SON->OCC[k];
						}
					}
				}
			}
			else
			{
				// Search next
				GetTree4(SON, MOT, LEN, POS + 1, DIM, TYP, SEQ, OCC);
			}
		}

		// Next
		SON = SON->BRO;
	}
}

// ---------------------------------------------------------------------------
// Get a motif from Suffix Tree
void GetTree5(CTree *Tree, char *MOT, int LEN, int POS, int DIM, char *SEQ, int *NBR, FILE *file)
{
	// get first son
	CTree *SON = Tree->SON;

	// Scroll over sons
	while (SON != NULL)
	{
		// Check if matches suffix tree
		if ((SON->RES == MOT[POS]) || (SON->RES == '.' && MOT[POS] != '\t' && MOT[POS] != '_'))
		{
			SEQ[POS + 0] = SON->RES;
			SEQ[POS + 1] = '\0';

			// Search next
			if (POS == LEN - 1)
			{
				// Check if motif
				if (SON->IDX > 0)
				{
					// Inc nbr of motifs
					*NBR += 1;

					// Write sequences
					fprintf_s(file, "%s", SEQ);

					// Write sequence occurrences
					for (int j=0; j<DIM; j++)
					fprintf_s(file, "\t%d", (int)SON->OCC[j]);

					// New line
					fprintf_s(file, "\n");
				}
			}
			else
			{
				GetTree5(SON, MOT, LEN, POS + 1, DIM, SEQ, NBR, file);
			}
		}

		// Next
		SON = SON->BRO;
	}
}

// ---------------------------------------------------------------------------
// Get a motif from Suffix Tree
void GetTree6(CTree *Tree, char *MOT, int LEN, int POS)
{
	// get first son
	CTree *SON = Tree->SON;

	// Scroll over sons
	while (SON != NULL)
	{
		// Check if matches suffix tree
		if ((SON->RES == MOT[POS]) || (SON->RES != '\t' && MOT[POS] == '.'))
		{
			// Search next
			if (POS == LEN - 1)
			{
				// Set sufix
				if (SON->IDX == 1) SON->IDX = 2;
			}
			else
			{
				GetTree6(SON, MOT, LEN, POS + 1);
			}
		}

		// Next
		SON = SON->BRO;
	}
}

// ---------------------------------------------------------------------------
// Get a motif from Suffix Tree
void GetTree7(CTree *Tree, int LEN, int POS, int DIM, char *SEQ, int *NBR, int *TYP, FILE *file)
{
	// get first son
	CTree *SON = Tree->SON;

	// Scroll over sons
	while (SON != NULL)
	{
		SEQ[POS + 0] = SON->RES;
		SEQ[POS + 1] = '\0';

		// Check if motif
		if (POS == LEN - 1)
		{
			if (SON->IDX == 2)
			{
				// Set back suffix
				SON->IDX = 1;

				// Init occ in pos/neg
				int occ1 = +INT_MAX;
				int occ2 = -INT_MAX;

				for (int j=0; j<DIM; j++)
				{
					switch (TYP[j])
					{
						case +1: occ1 = __MIN(occ1, SON->OCC[j]); break;
						case -1: occ2 = __MAX(occ2, SON->OCC[j]); break;
					}
				}

				// Inc nbr of sequences
				*NBR += 1;

				// Write sequences
				fprintf_s(file, "%s", SEQ);

				// Write sequence occurrences
				for (int j=0; j<DIM; j++) fprintf_s(file, "\t%d", (int)SON->OCC[j]); fprintf_s(file, "\n");
			}
		}
		else
		{
			// Next
			GetTree7(SON, LEN, POS + 1, DIM, SEQ, NBR, TYP, file);
		}

		// Next
		SON = SON->BRO;
	}
}

// ---------------------------------------------------------------------------
// Sort suffix tree (decreasing)
void SortTree1(CTree *Tree)
{
	// get first
	CTree *SON1 = Tree->SON;

	// Scroll over
	while (SON1 != NULL)
	{
		// get first
		CTree *SON2 = SON1->BRO;

		// Scroll over
		while (SON2 != NULL)
		{
			// check and swap residues
			if (aa2id(SON1->RES) > aa2id(SON2->RES))
			{
				SWAP(SON1->IDX, SON2->IDX);
				SWAP(SON1->RES, SON2->RES);
				SWAP(SON1->OCC, SON2->OCC);
				SWAP(SON1->SON, SON2->SON);
			}

			// Next
			SON2 = SON2->BRO;
		}

		// sort next
		SortTree1(SON1);

		// Next
		SON1 = SON1->BRO;
	}
}

// ---------------------------------------------------------------------------
// Sort suffix tree (decreasing)
void SortTree2(CTree *Tree)
{
	// get first
	CTree *SON1 = Tree->SON;

	// Scroll over
	while (SON1 != NULL)
	{
		// get first
		CTree *SON2 = SON1->BRO;

		// Scroll over
		while (SON2 != NULL)
		{
			// check and swap residues
			if (aa2id(SON1->RES) < aa2id(SON2->RES))
			{
				SWAP(SON1->IDX, SON2->IDX);
				SWAP(SON1->RES, SON2->RES);
				SWAP(SON1->OCC, SON2->OCC);
				SWAP(SON1->SON, SON2->SON);
			}

			// Next
			SON2 = SON2->BRO;
		}

		// sort next
		SortTree2(SON1);

		// Next
		SON1 = SON1->BRO;
	}
}

// ---------------------------------------------------------------------------
// Free Suffix Tree
void EmptyTree(CTree *Tree)
{
	// next
	if (Tree != NULL)
	{
		// Scroll horizontally
		EmptyTree(Tree->BRO);

		// Scroll vertically
		EmptyTree(Tree->SON);

		// free memory
		delete[] Tree->OCC; delete Tree->BRO; delete Tree->SON;
	}
}

// ---------------------------------------------------------------------------
// Discover motifs in positives
static unsigned __stdcall DiscoverMotifs(void *arg)
{
	// Get local pointer of arg
	CThread *ptr = (CThread*)arg;

	// Check for invalid threads
	if (ptr->THR > 0 && ptr->CUR == 1) return 0;

	// Counters
	int nbr0 = 0; int nbr1 = 0; int nbr2 = 0; double cut = 0.0f;

	// Parsers
	char BUF[MAXBUF]; char MOT[MAXBUF]; char SEQ[MAXBUF];

	// suffix tree
	CTree *Tree = NULL;

	////////////////////////////
	////////// step 1 //////////
	////////////////////////////

	// Init Suffix Tree
	Tree        = new CTree;
	Tree->IDX   = -1;
	Tree->RES   = '#';
	Tree->OCC   = NULL;
	Tree->BRO   = NULL;
	Tree->SON   = NULL;

	// Build suffix tree with previous motifs
	if (ptr->CUR > 1)
	{
		// Scroll over threads
		for (int i=0; i<NBRTHR; i++)
		{
			// Check for invalid threads
			if (i > 0 && ptr->CUR == 2) break;

			// Open input file																		// Check if file exist
			sprintf_s(BUF, MAXBUF, ".\\TEMP\\TMP9_%d_%d.dat", ptr->CUR - 1, i);					FILE *file1 = openfile(BUF, "rt");

			while (!feof(file1))
			{
				// Read current line
				if (fgets(BUF, MAXBUF, file1) == NULL) break;

				// Check if current thread
				if ((++nbr0) % NBRTHR != ptr->THR) continue;

				// Read current motifs
				sscanf_s(BUF, "%s", &MOT, MAXBUF);

				// Add motif in the suffix tree
				SetTree0(Tree, MOT, length(MOT), 0, 0, 0, 0);
			}

			// Close stream;
			fclose(file1);
		}
	}

	////////////////////////////
	////////// step 2 //////////
	////////////////////////////

	// Build siffix tree, scroll over all positives
	for (int i=0; i<ptr->NB1; i++)
	{
		// Add motif in the suffix tree
		SetTree1(Tree, ptr->SEQ1[i], ptr->CUR, 0, ptr->DIM, 0, 0, ptr->CDR);
	}

	// sort suffix tree
	SortTree2(Tree);

	////////////////////////////
	////////// step 3 //////////
	////////////////////////////

	// Build siffix tree, scroll over all positives
	for (int i=0; i<ptr->NB2; i++)
	{
		// Update occurrences
		GetTree2(Tree, ptr->SEQ2[i], ptr->CUR, 0, ptr->DIM, ptr->OCC2[i], ptr->TYP);
	}

	////////////////////////////
	////////// step 4 //////////
	////////////////////////////

	// Open output file																				// Check if file exist
	sprintf_s(BUF, MAXBUF, ".\\TEMP\\TMP7_%d_%d.dat", ptr->CUR, ptr->THR);							FILE *file2 = openfile(BUF, "wt");

	// Get a motif from Suffix Tree
	GetTree3(Tree, BUF, ptr->CUR, 0, ptr->DIM, ptr->TYP, &nbr1, file2);

	// Close output streams
	fclose(file2);

	// free Tree
	EmptyTree(Tree);

	////////////////////////////
	////////// step 5 //////////
	////////////////////////////

	// Open input file																				// Open file for reading
	sprintf_s(BUF, MAXBUF, ".\\TEMP\\TMP7_%d_%d.dat", ptr->CUR, ptr->THR);							FILE *file3 = openfile(BUF, "rt");

	// Open input file																				// Check if file exist
	sprintf_s(BUF, MAXBUF, ".\\TEMP\\TMP8_%d_%d.dat", ptr->CUR, ptr->THR);							FILE *file4 = openfile(BUF, "wt");

	while (!feof(file3))
	{
		// Read current line
		if (fgets(BUF, MAXBUF, file3) == NULL) break;

		// trim buffer
		sscanf_s(BUF, "%[^\n]\n", BUF, MAXBUF);

		// trim buffer
		sscanf_s(BUF, "%[^\t]\t", MOT, MAXBUF);

		// Check for wildcards
		if (memchr(MOT, '.', length(MOT)) == NULL || MOT[length(MOT) - 1] == '_') { fprintf_s(file4, "%s\n", BUF); continue; }

		// Init amino acids frequencies
		int ***OCC = Alloc3D(ptr->CUR, ptr->DIM, AASYMB + 1, 0);

		// Add motif in the suffix tree
		GetTree4(ptr->TREE, MOT, ptr->CUR, 0, ptr->DIM, ptr->TYP, SEQ, OCC);

		// Init nbr of residues with min frequency
		int NAA = INT_MAX;

		// Check amino acids frequencies at wildcards
		for (int i=0; i<ptr->CUR; i++)
		{
			// Skeep if not wildcard
			if (MOT[i] != '.') continue;

			// Init nbr of aa
			int naa = 0;

			// Check nbr of significant aa
			for (int j=0; j<AASYMB; j++)
			{
				// Init occurrences
				int occ1 = +INT_MAX;
				int occ2 = -INT_MAX;

				// Find aa min occurrences in pos
				for (int k=0; k<ptr->DIM; k++)
				{
					switch (ptr->TYP[k])
					{
						case +1: occ1 = __MIN(occ1, OCC[i][k][j]); break;
						case -1: occ2 = __MAX(occ2, OCC[i][k][j]); break;
					}
				}

				// Update frequent aa
				naa += occ1 >= MINOCC || occ2 >= MINOCC;
			}

			// Update freq aa
			if ((NAA = __MIN(NAA, naa)) < 2) break;
		}

		// Check aa frequencies
		if (NAA >= 2) { nbr2 += 1; fprintf_s(file4, "%s\n", BUF); }

		// Free memory
		for (int j=0; j<ptr->CUR; j++)
		{
			// Free memory
			for (int k=0; k<ptr->DIM; k++)
			delete[] OCC[j][k];
			delete[] OCC[j];
		}

		// Free memory
		delete[] OCC;
	}

	// Close stream
	fclose(file3);
	fclose(file4);

	// sampling cutoff (25%)
	cut = 25 * (1 - exp(-pow(nbr2, 2)/pow(1E+04, 2)));

	////////////////////////////
	////////// step 6 //////////
	////////////////////////////

	// Open input file																				// Open file for reading
	sprintf_s(BUF, MAXBUF, ".\\TEMP\\TMP8_%d_%d.dat", ptr->CUR, ptr->THR);							FILE *file5 = openfile(BUF, "rt");

	// Open input file																				// Check if file exist
	sprintf_s(BUF, MAXBUF, ".\\TEMP\\TMP9_%d_%d.dat", ptr->CUR, ptr->THR);							FILE *file6 = openfile(BUF, "wt");

	while (!feof(file5))
	{
		// Read current line
		if (fgets(BUF, MAXBUF, file5) == NULL) break;

		// trim buffer
		sscanf_s(BUF, "%[^\n]\n", BUF, MAXBUF);

		// trim buffer
		sscanf_s(BUF, "%[^\t]\t", MOT, MAXBUF);

		// Check for wildcards separator and sampling
		if (memchr(MOT, '.', length(MOT)) == NULL || MOT[length(MOT) - 1] == '_' || rand() % 100 >= cut - 1) fprintf_s(file6, "%s\n", BUF);
	}

	// Close stream
	fclose(file5);
	fclose(file6);

	// terminate thread
	_endthreadex(0);

	// Return value
	return 0;
}
