#include <stdio.h>
#include <stdlib.h>
#include<stdint.h>
#include </pub/include/mysql/mysql.h>

#define myMalloc(type, n) (type *)malloc((n) * sizeof (type))
#define AAMaxCodons 6
extern char * saveString (const char *ptr);
extern char * saveChr(const char * str);

//////////////////////////////////////////////////////////////
typedef struct {
	char aa;
	int numCodons;
	char codons[AAMaxCodons][4];
}GeneticCode;

typedef struct {
	char symbol;
	int num;
	char * map;
}IUPAC_map;	

////////////////////////////////////////////////////////////
// based on the Variants table in mysql
typedef struct { 	// based on the Variants table in mysql
	//char *pop;		// the population ID
	char * chr;  	// the chromosome name
	int pos;		// the position of variant
	char * type; 	// the type of variant
	char * ref;		// the reference bases in this position
	char * alt;		// the alteration bases in this position
	int offset;		// the offset from start to variant position
	struct Variant * next;
}Variant;

typedef struct {
	int id, start, end, phase;
	char * type;
	struct Variant * pVariants;
	struct Exon * next;
	char * name;
}Exon;

typedef struct {
	char * name;
	int id, start, end;
	int regionLen;
	char * type;
	struct Variant * pVariants;
	struct Other_region * next;
	char *desc;
}Other_region;
	
typedef struct {
	char * name;
	int id, start, end;
	int utr5Len, utr3Len, cdsLen, intronLen;
	char * cds;
	char * aas;
	struct  Exon * pExons;
	struct  Transcript * next;
	char *desc;
}Transcript;

// based on the table in GBrowse
typedef struct {
	char * name;
	int seqid, id, start, end, strand, tier, bin, indexed;
	char * type;
	int numTranscripts;
	int numProteins;
	int numCDSs; // numCDSs <= numProteins <= numTranscript
	int numOthers;
	struct Transcript * pTranscripts;
	struct Other_region * pOther_region;
	struct Feature * next;
}Feature;
	
typedef struct {
	char * chr;
	int id, start, end, strand, tier, bin, indexed, seqid;
	long numVariants;
	int numFeatures;
	struct Variant * pVariants;
	struct Feature * pFeatures;
	char * pDNA;
	struct Chromosome * next;
}Chromosome;

////////////////////////////////////////////
// functions for genetic code
void getComplementDNA(char * pCDS);
char getAminoAcid(const char * codon);

void printVariants(Variant * pVariants);
Chromosome * getChromosome(char * chr, Chromosome * pGenome);
void printTranscript(Transcript * pT);
void printFeature( Feature * pFeature, int isDetailed);
void printExons( Exon * pExons, const char * pDNA);
void printTranscripts( Transcript * pTranscripts, const char * pDNA);
void printFeatures( Feature * pFeatures, const char * pDNA);
void printSequences(FILE * fp, const char * seq, int len);

void walkGenome( Chromosome * pGenome);
int getNumCDSs( Feature * pF);
int countFeatures( Feature * pFeatures, char * FType, int * numTs);
int * getNumTranscripts( Feature * pFeatures, int numFeatures, char * type);

 //Variant * createRegionVariants(char *pop, int offset, int len, int strand, char *chr, int start, int end);
 Variant * createRegionVariants(int offset, int len, int strand, char *chr, int start, int end);
//void createGeneticVariants( Chromosome * pGenome, char * pop);
void createGeneticVariants( Chromosome * pGenome);
 Other_region * createOther(int id, int *otherl);
 Exon * createExons(int id, int *utr5l, int *cdsl, int *utr3l, int* intronl);
 Transcript * createTranscripts(int id, int * pNumTranscripts);
 Feature * createFeatures( Feature * pFeatures,  Chromosome * pc, char * ft);
 Chromosome * createGenome();
void createTranscriptsFromDB( Chromosome * pGenome);
void createOtherRegionFromDB( Chromosome * pGenome);
void createFeaturesFromDB( Chromosome * pGenome);
void createSequencesFromDB( Chromosome * pGenome);

void freeVariants( Variant * pVariants);
void freeExons( Exon * pExons);
void freeOther_region( Other_region * pOther_region);
void freeTranscripts( Transcript * pTranscripts);
void freeFeatures( Feature * pFeatures);
void freeGenome( Chromosome * pGenome);

void createGenomeDNASequence(char * dnaFileName,  Chromosome * pGenome);
char * getCDS( Transcript * pT, int strand, const char * pDNA);
char * getAAs(const char * cds);
void getAllCDSs( Chromosome * pGenome);

void printCDS(FILE * cdsOut, const char * cds, int len);
void printAllCDSs(const char * fn,  Chromosome * pGenome, int len);

///////////////////////////////////////////////////////////////////////////////
// MySQL operations: defined in mysql.cpp
extern MYSQL * mysqlConnection, mysqlMysql;
extern MYSQL_RES * mysqlResult, * mysqlResult0;
extern MYSQL_ROW mysqlRow, mysqlRow0;

// int connectDB(const char * dbName);
int connectDB(char * hostName, char * dbName, char * user, char * passwd);

void disconnectDB(const char * dbName);
void pingAndReconnectToDB();
int executeSQL(const char * sql);
///////////////////////////////////////////////////////////////
