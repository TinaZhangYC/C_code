#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<unistd.h>
#include<ctype.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<search.h>

#include "variant.h"

#define MyBufLen 10240

extern char myBuffer[MyBufLen + 1];

// public functions
char * saveString(const char * str) {
	char * ptr = myMalloc(char, strlen(str) + 1);
	strcpy(ptr, str);
	return ptr;
}

////////////////////////////////////////////////
GeneticCode geneticCodes[21] = {
	{ 'A', 4, "GCA", "GCC", "GCG", "GCT" },
	{ 'C', 2, "TGC", "TGT" },
	{ 'D', 2, "GAC", "GAT" },
	{ 'E', 2, "GAA", "GAG" },
	{ 'F', 2, "TTC", "TTT" },
	{ 'G', 4, "GGA", "GGC", "GGG", "GGT" },
	{ 'H', 2, "CAC", "CAT" },
	{ 'I', 3, "ATA", "ATC", "ATT" },
	{ 'K', 2, "AAA", "AAG" },
	{ 'L', 6, "TTA", "TTG", "CTA", "CTC", "CTG", "CTT" },
	{ 'M', 1, "ATG" },
	{ 'N', 2, "AAC", "AAT" },
	{ 'P', 4, "CCA", "CCC", "CCG", "CCT" },
	{ 'Q', 2, "CAA", "CAG" },
	{ 'R', 6, "AGA", "AGG", "CGA", "CGC", "CGG", "CGT" },
	{ 'S', 6, "AGC", "AGT", "TCA", "TCC", "TCG", "TCT" },
	{ 'T', 4, "ACA", "ACC", "ACG", "ACT" },
	{ 'V', 4, "GTA", "GTC", "GTG", "GTT" },
	{ 'W', 1, "TGG" },
	{ 'Y', 2, "TAC", "TAT" },
	{ '*', 3, "TAA", "TGA", "TAG" }	};

// IUPAC nucleotide code	Base
// A	Adenine
// C	Cytosine
// G	Guanine
// T (or U)	Thymine (or Uracil)
IUPAC_map iupac[14] = {
	{'A', 1, "A"}, {'C', 1, "C"}, {'G', 1, "G"}, {'T', 1, "T"}, 
	{'R', 2, "AG"}, {'Y', 2, "CT"}, {'S', 2, "GC"}, {'W', 2, "AT"}, {'K', 2, "GT"}, {'M', 2, "AC"},
	{'B', 3, "CGT"}, {'D', 3, "AGT"}, {'H', 3, "ACT"}, {'V', 3, "ACG"}};

char getNucleotide(char symbol, char ref) {
	for(int i=0; i<14; ++i)
		if(iupac[i].symbol == symbol) {
			if(iupac[i].num == 1) return iupac[i].map[0];
			for(int j=0; j<iupac[i].num; ++j)
				if(iupac[i].map[j] != ref)
					return iupac[i].map[j];
		}			
	return 'X';
}			
			
char getAminoAcid(const char * codon) {
	for(int i=0; i<21; ++i) 
		for(int j=0; j<geneticCodes[i].numCodons; ++j)
			if(!strncmp(geneticCodes[i].codons[j], codon, 3))
				return geneticCodes[i].aa;
	
	return 'X'; // wrong?
}


char * isSynonymous(const char * codon1, const char * codon2, char *symbol) {
	char x = getAminoAcid(codon1);
	char y = getAminoAcid(codon2);
	
/*	
	if(y == '*') 
	{
		sprintf(symbol,"%s","Premature");
		return "Premature";
	}
*/
	if(x == y && x != 'X') 
	{
		sprintf(symbol,"Synonymous: %c->%c", x, y);
		return symbol;
	}
	else 
	{ 
		sprintf(symbol, "Nonsynonymous: %c->%c", x, y); 
		
		//printf("%s->%s\t%s\n",codon1, codon2, symbol);
		return symbol; 
	}
}

void getCodon(const char * cds, int strand, int offset, char * ptr) {
	for(int i = 0; i < 3; ++i)
		* ptr ++ = cds[offset - offset % 3 + i];
}
				
void printAAs(const char * cds) {
	int len = strlen(cds);
	int l = 0;
	while(l < len) {
		printf("%c",getAminoAcid(cds + l));
		l += 3;
	}
}

char getACGT(char ch) {
	switch(ch) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		default : break;	
	}
	return 'N';
}
			
			
void getComplementDNA(char * pCDS) {
	int len = strlen(pCDS), i;
	for(i=0; i<len; ++i)
		pCDS[i] = getACGT(pCDS[i]);
		
	for(i=0; i<len/2; ++i) {
		char ch = pCDS[len - i - 1];
		pCDS[len - i - 1] = pCDS[i];
		pCDS[i] = ch;
	}
}
		
////////////////////////////////////////////////
void printLine(int num, char symbol) {
    for(int i=0; i<num; ++i)
        printf("%c", symbol);
    printf("\n");
}


void printmRNAVariants(const struct Chromosome * pGenome, char * fn) {
	FILE * fp = fopen(fn, "w+");
	fprintf(fp, "Chr\tPos\tType\tREF\tALT\tFname\tEname\tStrand\tEtype\tVarOffset\tFeatureLen\tNucleChange\tCodonChange\tClass\tDescription\n");
	Chromosome * pChr = pGenome;
	while(pChr != NULL) {
		Feature * pF = pChr->pFeatures;
		while(pF != NULL) {
			Transcript * pT = pF->pTranscripts;
			while(pT != NULL) {
				Exon * pE = pT->pExons;
				while(pE != NULL) {
				    // only dealing with utr5, cds, and utr3
					if(!strcmp(pE->type, "exon:")) { pE = pE->next; continue; }
					int len;
					if(!strcmp(pE->type, "CDS:")) len = pT->cdsLen;
					else if(!strcmp(pE->type, "intron:")) len = pT->intronLen;
					else if(!strcmp(pE->type, "five_prime_utr:"))  len = pT->utr5Len;
					else if(!strcmp(pE->type, "three_prime_utr:")) len = pT->utr3Len;
					
					Variant * pV = pE->pVariants;
					
					//printf("chr:%s\ttranscript:%s\n",pChr->chr, pT->name);
					//printVariants(pV);
					
					while(pV != NULL) {
						fprintf(fp, "%s\t%u\t%s\t%s\t%s\t%s\t%s\t%c\t%s\t", pChr->chr, pV->pos, pV->type, pV->ref, pV->alt, pF->name, pE->name, (pF->strand == 1 ? '+' : '-'), pE->type);
						//int AANo=ceil(pV->offset / 3);
						fprintf(fp, "%d\t%d\t", pV->offset, len);
						if(!strcmp(pE->type, "CDS:")) {
							if(!strcmp(pV->type, "SNP")) { // Snp  
									char codon[4], altc[4];
									for(int i = 0; i < 3; ++i) {
										codon[i] = pT->cds[pV->offset - pV->offset % 3 + i];
										altc[i] = codon[i];
									}	
									codon[3] = altc[3] = '\0';
									
									char s = getNucleotide(pV->alt[0], pV->ref[0]);
									altc[pV->offset % 3] = (pF->strand == 1 ? s : getACGT(s));
									char symbol[128];
									isSynonymous(codon, altc, symbol);
									fprintf(fp, "%c<->%c\t%s<->%s\t%s\t%s\n", pT->cds[pV->offset], altc[pV->offset % 3], codon, altc, symbol, pT->desc);
								}
						   else // Indel
						   	fprintf(fp, "\t\t\t%s\n", pT->desc);						
						} 
						else 
							fprintf(fp, "\t\t\t%s\n", pT->desc);
						pV = pV->next;
					}//while PV
					pE = pE->next;	
				}//while PE		
				pT = pT->next;
			}//while PT
			Other_region * pO = pF->pOther_region;
			while( pO != NULL )
			{
				Variant * pV = pO->pVariants;
                while(pV != NULL) {
                        fprintf(fp, "%s\t%u\t%s\t%s\t%s\t%s\t%s\t%c\t%s\t", pChr->chr, pV->pos, pV->type, pV->ref, pV->alt, pF->name, pO->name, (pF->strand == 1 ? '+' : '-'), pO->type);
						//fprintf(fp, "%s\t%s\t", pO->type, pV->alt); 
                        fprintf(fp, "%d\t%d\t\t\t\t%s\n", pV->offset, pO->regionLen, pO->desc);	
					pV = pV->next;
				}//while PV
				pO = pO->next;
			}//while pO
			pF = pF->next;
		}//while PF
		pChr = pChr->next;
	}//while pChr
	fclose(fp);
}				

////////////////////////////////////////////////////////
// main entry
int main(int argc, char ** argv)
{
	if(argc < 6) { printf("Usage: ./annVar 192.168.3.165  Ecoli_pop_pan  zhangyc	zhangyc123  Ecoli_pop_pan_varAnn_raw.txt\n	host:	the name/IP of host\n	db:	the name of DB\n	user:	the user name\n	passwd:	the password of user to access the DB\n	output:	the genetic variants output file\n"); return 1; }
	 
	printf(">>> Parameters: ================================\n");
	int i;
	for(i=1; i<argc; ++i)
		printf("   %d) %s\n", i, argv[i]);
	
	printf("<<< ============================================\n\n");
	const char * hostName = saveString(argv[1]);
 	const char * dbName = saveString(argv[2]);
	const char * userName = saveString(argv[3]);
	const char * userPasswd = saveString(argv[4]);
	//const char * popName = saveString(argv[5]);
	const char * output = saveString(argv[5]);
	
	// Debug
	//int test=atoi(argv[7]);
	//printf("test atoi:%u\t %u\n",test,argv[7]);

	time_t start, end;
	
	
    /////////////////////////////////////////////////////////	
  	time(&start);
  	connectDB(hostName, dbName, userName, userPasswd);
  	Chromosome * pGenome = createGenome();


 	createFeaturesFromDB(pGenome);  	

	//walkGenome(pGenome);

	createTranscriptsFromDB(pGenome);
	createOtherRegionFromDB(pGenome);
  	createGeneticVariants(pGenome);
  	createSequencesFromDB(pGenome);
  	getAllCDSs(pGenome);
	// printAllCDSs("CDS.fasta", pGenome);
	printmRNAVariants(pGenome, output);
  	// walkGenome(pGenome);
    
	freeGenome(pGenome);
	return 0;
}

//////////////////////////////////////////////////////////
