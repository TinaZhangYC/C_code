#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdint.h>
#include<unistd.h>
#include<ctype.h>
#include<math.h>
#include<time.h>
#include<search.h>

#include "variant.h"

#define MyBufLen 10240
char myBuffer[MyBufLen + 1];

/////////////////////////////////////////////////////////////////////
////////////////////////////////////////
// functions		
void printVariants(Variant * pVariants) {
	Variant * pV = pVariants;
	while(pV != NULL) { 
		printf("   Chr:%s\tPos:%u\tType:%s\tRef:%s\tAlt:%s\n", pV->chr, pV->pos, pV->type, pV->ref, pV->alt);
		pV = pV->next;
	}		           		
}

Chromosome * getChromosome(char * chr, Chromosome * pGenome) {
	Chromosome * pChr = pGenome;
	while(pChr) {
	if(!strcmp(chr, pChr->chr)) return pChr;
		pChr = pChr->next;
	}
	return NULL;
}
			
void printTranscript(Transcript * pT) {
	printf("       -> %s [%u, %u]\n", pT->name, pT->start, pT->end);
	Exon * e = pT->pExons;
	int n = 0;
	while(e != NULL) {
		printf("		%u: %s, %u, %d\n", ++n, e->type, e->start, e->end);
		e = e->next;
	}	
}
	
void printFeature(Feature * pFeature, int isDetailed) {
	printf("    -> %s(%c): [%lu, %lu], %u transcripts; (type: %s)\n", 
	        pFeature->name, pFeature->strand, pFeature->start, pFeature->end, pFeature->numTranscripts, pFeature->type);
	
	if(!isDetailed) return;
		
	Transcript * pT = pFeature->pTranscripts;
	while(pT != NULL) {
		printTranscript(pT);
		pT = pT->next;
	}	
}

int getNumCDSs(Feature * pF) {
	int num = 0;
	Transcript * pT = pF->pTranscripts;
	while(pT) {
		if(pT->aas != NULL) ++ num;
		pT = pT->next;
	}
	return num;
}		
	
void printExons(Exon * pExons, const char * pDNA) {
	Exon * pE = pExons;
	
	while(pE != NULL) {
		printf("			%s, [%u, %u] ", pE->type, pE->start, pE->end);
		printVariants(pE->pVariants);
		printf("\n");	
		pE = pE->next;
	}
	printf("\n");
}

void printTranscripts(Transcript * pTranscripts, const char * pDNA) {
	Transcript * pT = pTranscripts;
	int n = 0;
	while(pT != NULL) {
		printf("		%3u: %s, [%lu, %lu]:\n", ++n, pT->name, pT->start, pT->end);
		printExons(pT->pExons, pDNA);
		pT = pT->next;
	}
}

void printFeatures(Feature * pFeatures, const char * pDNA) {
	Feature * pF = pFeatures;
	int n = 0;
	while(pF != NULL) {
		printf("	%5u: %s[%c](%s), [%lu, %lu], %u transcripts.\n", ++n, 
		       pF->name, (pF->strand == 1 ? '+' : '-'), pF->type, pF->start, pF->end, pF->numTranscripts);
		printTranscripts(pF->pTranscripts, pDNA);
		pF = pF->next;
	}
}

void printSequences(FILE * fp, const char * seq, int len) {
	int i = 0;
	for(i = 0; i < strlen(seq) / len; ++i) {
		for(int j = 0; j < len; ++ j)
			fprintf(fp, "%c", seq[len * i + j]);
		fprintf(fp, "\n");
	}
	if(strlen(seq) % len) {
		for(int j = 0; seq[j] != '\0'; ++ j)
			fprintf(fp, "%c", seq[len * i + j]);
		fprintf(fp, "\n");
	}
}		
	 					 
void walkGenome(Chromosome * pGenome) {
	Chromosome * pChr = pGenome;
	while(pChr != NULL) {
		printf("Chr: %s, %u features, %u variants.\n", pChr->chr, pChr->numFeatures, pChr->numVariants);
		printFeatures(pChr->pFeatures, pChr->pDNA);
		
		// FILE * fp = fopen("DNA.fasta", "wt");
		// printSequences(fp, pChr->pDNA);
		// fclose(fp);
		
		pChr = pChr->next;
	}
}

int countFeatures(Feature * pFeatures, char * FType, int *numTs) {
	int n = 0;
	Feature * pF = pFeatures;
	while(pF != NULL) {
		if(!strcmp(pF->type, FType)) {
			++ n;
			numTs += pF->numTranscripts;
		}	
		pF = pF->next;
	}
	return n;
}

int * getNumTranscripts(Feature * pFeatures, int numFeatures, char * type) {
	int * pNumTranscripts = myMalloc(int, numFeatures);
	int n = 0;
	Feature * pF = pFeatures;
	while(pF != NULL) {
		if(!strcmp(pF->type, type))
			pNumTranscripts[n++] = pF->numTranscripts;
		pF = pF->next;
	}
	return pNumTranscripts;
}
						
//////////////////////////////////////////////////////////		
// creating Variants from db
 Variant * createRegionVariants(int offset, int len, int strand, char *chr, int start, int end) {
	
	//printf("Enter createRegionVariants\n");
	
	sprintf(myBuffer, "SELECT *  FROM Variants AS A WHERE A.Chr = '%s' AND (A.Pos BETWEEN %u AND %u) ORDER BY A.Pos", chr, start, end);

	// printf(myBuffer);
  	executeSQL(myBuffer);
  	mysqlResult = mysql_store_result(mysqlConnection);
  
  	Variant * pVariants = NULL, * prevVariant;
	while((mysqlRow = mysql_fetch_row(mysqlResult)) != NULL)
  	{
  		Variant * pVariant = myMalloc( Variant, 1);
		pVariant->next = NULL;
  		pVariant->chr  = saveString(mysqlRow[0]); // chromosome name
  		pVariant->pos  = atoi(mysqlRow[1]);		  // position
  		pVariant->type = saveString(mysqlRow[2]); // mutation type
  		pVariant->ref  = saveString(mysqlRow[3]); // the reference bases in this position
  		pVariant->alt = saveString(mysqlRow[4]); // the alteration bases in this position
  		
		pVariant->offset = offset + pVariant->pos - start;  // mutation pos in CDS region, exon1-->exon2-->exon3
		//pVariant->offset = offset + pVariant->pos - start + 1;  // mutation pos in CDS region, exon1-->exon2-->exon3
		if(strand == -1) pVariant->offset = len - pVariant->offset - 1;
		//if(strand == -1) pVariant->offset = len - pVariant->offset + 1;
  		
		// adding to the list
		if(pVariants == NULL) pVariants = pVariant;
	  	else prevVariant->next = pVariant;
	
	 	prevVariant = pVariant;
    }
  
  mysql_free_result(mysqlResult);
    	  		 
  return pVariants;
}

// creating Variants from db
void createGeneticVariants( Chromosome * pGenome) {
  	Chromosome * pChr = pGenome;
	while(pChr)
  	{
		Feature * pF = pChr->pFeatures;
		while(pF != NULL) {
			Transcript * pT = pF->pTranscripts;
			Other_region * pR = pF->pOther_region;
			while(pT != NULL) { 
				int utr5Offset = 0, cdsOffset = 0, utr3Offset = 0, intronOffset = 0;
				 Exon * pE = pT->pExons;
				while(pE != NULL) { 
					if(!strcmp(pE->type, "exon:")) { pE = pE->next; continue; }
					
					int len = pE->end - pE->start + 1;  					
					if(!strcmp(pE->type, "five_prime_utr:")) {
  						pE->pVariants = createRegionVariants(utr5Offset, pT->utr5Len, pF->strand, pChr->chr, pE->start, pE->end);
						utr5Offset += len;
					}	
					else if(!strcmp(pE->type, "CDS:")) {
  						pE->pVariants = createRegionVariants(cdsOffset, pT->cdsLen, pF->strand, pChr->chr, pE->start, pE->end);					
						cdsOffset += len;
					}	
					else if(!strcmp(pE->type, "three_prime_utr:")) {
  						pE->pVariants = createRegionVariants(utr3Offset, pT->utr3Len, pF->strand, pChr->chr, pE->start, pE->end);					
						utr3Offset += len;
					} else { // intron	
						pE->pVariants = createRegionVariants(intronOffset, pT->intronLen, pF->strand, pChr->chr, pE->start, pE->end); 
						intronOffset += len;
					}	
  					
  					pE = pE->next;
  				}	
  				pT = pT->next;
  			}
			int RBSoffset = 0, STSoffset = 0, repRegoffset = 0, repoffset = 0, sRNAoffset = 0, termoffset = 0, ncRNAoffset = 0, MEoffset = 0, tRNAoffset = 0, GISoffset = 0, rRNAoffset = 0, tmRNAoffset = 0;
			while(pR != NULL)
			{
				int len = pR->end - pR->start + 1;
				if(!strcmp(pR->type, "RBS:")){
					pR->pVariants = createRegionVariants(RBSoffset, pR->regionLen, pF->strand, pChr->chr, pR->start, pR->end); 
					//repOrioffset += len;
				}
				if(!strcmp(pR->type, "STS:")){
					pR->pVariants = createRegionVariants(STSoffset, pR->regionLen, pF->strand, pChr->chr, pR->start, pR->end); 
				}
				if(!strcmp(pR->type, "repeat_region:")){
					pR->pVariants = createRegionVariants(repRegoffset, pR->regionLen, pF->strand, pChr->chr, pR->start, pR->end); 
					//repRegoffset += len;
				}
				if(!strcmp(pR->type, "rep_origin:")){
					pR->pVariants = createRegionVariants(repoffset, pR->regionLen, pF->strand, pChr->chr, pR->start, pR->end); 
				}
				if(!strcmp(pR->type, "sRNA:")){
					pR->pVariants = createRegionVariants(sRNAoffset, pR->regionLen, pF->strand, pChr->chr, pR->start, pR->end); 
					
					//LTRoffset += len;
				}
				if(!strcmp(pR->type, "terminator:")){
					pR->pVariants = createRegionVariants(termoffset, pR->regionLen, pF->strand, pChr->chr, pR->start, pR->end);
				}
				if(!strcmp(pR->type, "ncRNA:")){
					pR->pVariants = createRegionVariants(ncRNAoffset, pR->regionLen, pF->strand, pChr->chr, pR->start, pR->end);
				}
				if(!strcmp(pR->type, "mobile_element:")){
					pR->pVariants = createRegionVariants(MEoffset, pR->regionLen, pF->strand, pChr->chr, pR->start, pR->end);
				}
				if(!strcmp(pR->type, "tRNA:")){
					pR->pVariants = createRegionVariants(tRNAoffset, pR->regionLen, pF->strand, pChr->chr, pR->start, pR->end);
				}
				if(!strcmp(pR->type, "GIS:")){
					pR->pVariants = createRegionVariants(GISoffset, pR->regionLen, pF->strand, pChr->chr, pR->start, pR->end);
				}
				if(!strcmp(pR->type, "rRNA:")){
					pR->pVariants = createRegionVariants(rRNAoffset, pR->regionLen, pF->strand, pChr->chr, pR->start, pR->end);
				}
				if(!strcmp(pR->type, "tmRNA:")){
					pR->pVariants = createRegionVariants(tmRNAoffset, pR->regionLen, pF->strand, pChr->chr, pR->start, pR->end);
				}
				pR = pR->next;
			}

  			pF = pF->next;	
		}
		pChr = pChr->next;
	}
}
////////////////////////////////////////////////////////////////////////////////////////

// creating Features from db
 Other_region * createOther(int id, int *otherl)
{	
	sprintf(myBuffer, "SELECT C.name, B.tag, A.id, A.start, A.end, D.attribute_value FROM feature AS A, typelist AS B, name AS C, attribute AS D, attributelist AS E WHERE  A.id = %u AND B.id = A.typeid AND A.id = C.id AND A.id = D.id AND D.attribute_id = E.id AND E.tag = 'Description' ORDER BY A.start", id);
  	executeSQL(myBuffer);
  
  	mysqlResult0 = mysql_store_result(mysqlConnection);
  
  	int otherLen = 0;
  	 Other_region * pOthers = NULL, * prevOther;
	while((mysqlRow0 = mysql_fetch_row(mysqlResult0)) != NULL)
  	{
  		 Other_region * pOther = myMalloc( Other_region, 1);
		pOther->next = NULL;
		pOther->name = saveString(mysqlRow0[0]);
  		pOther->type = saveString(mysqlRow0[1]);
  		pOther->id = atoi(mysqlRow0[2]);
		pOther->start = atoi(mysqlRow0[3]);
  		pOther->end   = atoi(mysqlRow0[4]);
		pOther->desc = saveString(mysqlRow0[5]);
  		pOther->pVariants = NULL;
		pOther->regionLen	= pOther->end - pOther->start +1;
		otherLen += pOther->regionLen;
  	
		
		// adding to the list
		if(pOthers == NULL) pOthers = pOther;
	  	else prevOther->next = pOther;
	
	 	prevOther = pOther;
   }
  
  mysql_free_result(mysqlResult0);
  otherl = otherLen;
  
  return pOthers;
}

 Exon * createExons(int id, int *utr5l, int *cdsl, int *utr3l, int *intronl) {
	
	//printf("Enter createExons\n");
	
	
	sprintf(myBuffer, "SELECT C.tag, A.id, A.start, A.end, D.name FROM feature AS A, parent2child AS B, typelist AS C, name AS D WHERE  B.id = %u AND A.id = B.child AND C.id = A.typeid AND D.id = A.id ORDER BY A.start", id);
  	executeSQL(myBuffer);
  
  	mysqlResult0 = mysql_store_result(mysqlConnection);
  
  	int utr5Len = 0, utr3Len = 0, cdsLen = 0, intronLen = 0;
  	 Exon * pExons = NULL, * prevExon, * prevE = NULL;
	while((mysqlRow0 = mysql_fetch_row(mysqlResult0)) != NULL)
  	{
  		Exon * pExon = myMalloc( Exon, 1);
		pExon->next = NULL;
  		pExon->type = saveString(mysqlRow0[0]);
  		pExon->id = atoi(mysqlRow0[1]);
		pExon->start = atoi(mysqlRow0[2]);
  		pExon->end   = atoi(mysqlRow0[3]);
  		pExon->name = saveString(mysqlRow0[4]);
  		pExon->pVariants = NULL;

		//printf("Exon id:%u\ttype:%s\tstart:%u\tend:%u\n",pExon->id, pExon->type, pExon->start, pExon->end);


		if(!strcmp(pExon->type, "five_prime_utr:"))
			utr5Len += pExon->end - pExon->start + 1;
		else if(!strcmp(pExon->type, "CDS:"))
			cdsLen += pExon->end - pExon->start + 1;
		else if(!strcmp(pExon->type, "three_prime_utr:"))
			utr3Len += pExon->end - pExon->start + 1;
		else { // !strcmp(pExon->type, "exon:")
			if(prevE != NULL) { // inserting an Intron
				 Exon * pE = myMalloc( Exon, 1);
				pE->next = NULL;
  				pE->type = saveString("intron:");
  				pE->id = -1;
				pE->start = prevE->start + 1;
  				pE->end   = pExon->start - 1;
  				pE->pVariants = NULL;
				prevExon->next = pE;
				prevExon = pE;
				intronLen += pE->end - pE->start + 1;
			}	
			prevE = pExon; 	
	 	}
	
		//printf("Exon Type:%s\tStart:%d\tEnd:%d\n",pExon->type, pExon->start, pExon->end);

  		// adding to the list
		if(pExons == NULL) pExons = pExon;
	  	else prevExon->next = pExon;
	
	 	prevExon = pExon;
  }
  
  mysql_free_result(mysqlResult0);
  *utr5l = utr5Len;
  *cdsl = cdsLen;
  *utr3l = utr3Len;
  *intronl = intronLen;
  
  return pExons;
}

 Transcript * createTranscripts(int id, int * pNumTranscripts) {
	sprintf(myBuffer, "SELECT C.name, A.id, A.start, A.end, D.attribute_value FROM feature AS A, parent2child AS B, name AS C, attribute AS D, attributelist AS E WHERE B.id = %u AND A.id = B.child AND C.id = A.id AND A.id = D.id AND D.attribute_id = E.id AND E.tag = 'Description'", id);
  	executeSQL(myBuffer);
  	mysqlResult = mysql_store_result(mysqlConnection);
  	* pNumTranscripts = mysql_num_rows(mysqlResult);
  
  	 Transcript * pTranscripts = NULL, * prevTranscript;
	while((mysqlRow = mysql_fetch_row(mysqlResult)) != NULL)
  	{
  		Transcript * pTranscript = myMalloc( Transcript, 1);
		pTranscript->next = NULL;
  		pTranscript->name = saveString(mysqlRow[0]);
  		pTranscript->id = atoi(mysqlRow[1]);
  		pTranscript->start = atoi(mysqlRow[2]);
  		pTranscript->end   = atoi(mysqlRow[3]);
  		pTranscript->desc = saveString(mysqlRow[4]);
		
  		pTranscript->pExons = createExons(pTranscript->id, &pTranscript->utr5Len, &pTranscript->cdsLen, &pTranscript->utr3Len, &pTranscript->intronLen);
  		pTranscript->cds = NULL;
  		pTranscript->aas = NULL;

		//printf("Transcript name:%s\tstart:%d\tend:%d\tCDS len:%d\n",pTranscript->name, pTranscript->start, pTranscript->end, pTranscript->cdsLen);

  		// adding to the list
		if(pTranscripts == NULL) pTranscripts = pTranscript;
		else prevTranscript->next = pTranscript;
	
		prevTranscript = pTranscript;
   }
  
  mysql_free_result(mysqlResult);
  
  return pTranscripts;
}
	
 Feature * createFeatures( Feature * pFeatures,  Chromosome * pc, char * ft) {
	//int len=strlen(pc->chr);
	//char *Chr=myMalloc(char, len-4 + 1);
	//Chr=saveChr(pc->chr);

	//printf("Enter createFeatures\npc->chr:%s\n",pc->chr);

	sprintf(myBuffer, "SELECT B.name, A.seqid, A.id, A.start, A.end, A.strand, A.tier, A.bin, A.indexed FROM feature AS A, name AS B, typelist AS C, locationlist AS D WHERE A.seqid = D.id AND D.seqname = '%s' AND C.tag = '%s' AND A.typeid = C.id AND B.id = A.id", pc->chr, ft);
  	executeSQL(myBuffer);
  	mysqlResult = mysql_store_result(mysqlConnection);
  	pc->numFeatures += mysql_num_rows(mysqlResult);
  
  	if(mysql_num_rows(mysqlResult)==0)
	{
		printf("The child number is 0!\n");
		return pFeatures;
	}

  	// printf("     Chr %s: %lu features\n", chr, * pNumFeatures);
  	Feature * prevFeature;
	while((mysqlRow = mysql_fetch_row(mysqlResult)) != NULL)
    {
  		Feature * pFeature = myMalloc( Feature, 1);
		pFeature->next = NULL;
  		pFeature->type = saveString(ft);

  		pFeature->name = saveString(mysqlRow[0]);
  		pFeature->seqid = atoi(mysqlRow[1]);
		pFeature->id = atoi(mysqlRow[2]);
  		pFeature->start = atoi(mysqlRow[3]);
  		pFeature->end   = atoi(mysqlRow[4]);
  		pFeature->strand = atoi(mysqlRow[5]);
  		pFeature->pTranscripts = NULL;
  		pFeature->pOther_region = NULL;
  		pFeature->numProteins = 0;
  		pFeature->numCDSs = 0;
  
		//printf("Feature name:%s\ttype:%s\n", pFeature->name, pFeature->type);

  		// adding to the list
		if(pFeatures == NULL) pFeatures = pFeature;
	  	else prevFeature->next = pFeature;
	
	 	prevFeature = pFeature;
  	}
  
  	mysql_free_result(mysqlResult);
  	return pFeatures;
}

 Chromosome * createGenome() {
	executeSQL("SELECT C.name, B.id, B.start, B.end, B.strand, B.tier, B.bin, B.indexed, B.seqid from typelist as A, feature as B, name as C where A.tag = 'chromosome:ADD' AND A.id = B.typeid AND C.id = B.id");
  	mysqlResult0 = mysql_store_result(mysqlConnection);
  
  	printf("   There are %u chromosomes in the DB.\n", mysql_num_rows(mysqlResult0));
  	Chromosome * pGenome = NULL, * prevChr;
	while((mysqlRow0 = mysql_fetch_row(mysqlResult0)) != NULL)
  	{
  		Chromosome * pChr = myMalloc( Chromosome, 1);
		pChr->next = NULL;
		pChr->pFeatures = NULL;
		pChr->pDNA = NULL;
  		pChr->chr = saveString(mysqlRow0[0]);
		pChr->id = atoi(mysqlRow0[1]);
  		pChr->start = atoi(mysqlRow0[2]);
  		pChr->end = atoi(mysqlRow0[3]);
  		pChr->strand = atoi(mysqlRow0[4]);
  		pChr->tier = atoi(mysqlRow0[5]);
  		pChr->bin = atoi(mysqlRow0[6]);
  		pChr->indexed = atoi(mysqlRow0[7]);
  		pChr->seqid = atoi(mysqlRow0[8]);
  		
		pChr->numVariants=0;
		pChr->numFeatures=0;
  		pChr->pVariants = NULL;

	  	
		//printf("chr:%s\tid:%d\n",pChr->chr,pChr->id);
		
		// adding to the list
		if(pGenome == NULL) pGenome = pChr;
		else prevChr->next = pChr;
	
		prevChr = pChr;
		// break; // for testing
  	}
  
  	mysql_free_result(mysqlResult0);
  
  	return pGenome;
}

void createTranscriptsFromDB( Chromosome * pGenome) {
	 Chromosome * pChr = pGenome;
	while(pChr != NULL) {
		time_t start, end;
	  	time(&start);
		int numTranscripts = 0;
		Feature * pF = pChr->pFeatures;
		
		//printf("Enter createTranscriptsFromDB\tChr:%s\n",pChr->chr);

		while(pF != NULL) {
			//printf("Feature type:%s\n",pF->type);
			if(!strcmp(pF->type, "gene:"))
			{
				//printf("Feature type:%s\tID:%d\n", pF->type, pF->id);
				pF->pTranscripts = createTranscripts(pF->id, &(pF->numTranscripts));
				numTranscripts += pF->numTranscripts;
			}
			pF = pF->next;
		}	

	  	time(&end);
  		printf("     Chr %s: %u features, %u transcripts costed %ld seconds.\n", 
  	  	pChr->chr, pChr->numFeatures, numTranscripts, end - start);

		pChr = pChr->next;
	}
}

void createOtherRegionFromDB( Chromosome * pGenome) {
	 Chromosome * pChr = pGenome;
	while(pChr != NULL) {
		time_t start, end;
	  	time(&start);
		int numOthers = 0;
		 Feature * pF = pChr->pFeatures;
		while(pF != NULL) 
		{
		 	if(strcmp(pF->type, "gene:"))
		 	{
				pF->pOther_region = createOther(pF->id, &(pF->numOthers));
				numOthers += pF->numOthers;
		 	}
			pF = pF->next;
		}

	  	time(&end);
  		printf("     Chr %s: %u features, %u other regions costed %ld seconds.\n", 
  	  	pChr->chr, pChr->numFeatures, numOthers, end - start);

		pChr = pChr->next;
	}
}

void createFeaturesFromDB( Chromosome * pGenome) {
	 Chromosome * pChr = pGenome;
	while(pChr != NULL) {
		
		//printf("Enter createFeaturesFromDB\n");

		pChr->pFeatures=createFeatures(pChr->pFeatures, pChr, "gene:");
		pChr->pFeatures=createFeatures(pChr->pFeatures, pChr, "RBS:");
		pChr->pFeatures=createFeatures(pChr->pFeatures, pChr, "STS:");
		pChr->pFeatures=createFeatures(pChr->pFeatures, pChr, "repeat_region:");
		pChr->pFeatures=createFeatures(pChr->pFeatures, pChr, "rep_origin:");
		pChr->pFeatures=createFeatures(pChr->pFeatures, pChr, "sRNA:");
		pChr->pFeatures=createFeatures(pChr->pFeatures, pChr, "terminator:");
		pChr->pFeatures=createFeatures(pChr->pFeatures, pChr, "ncRNA:");
		pChr->pFeatures=createFeatures(pChr->pFeatures, pChr, "mobile_element:");
		pChr->pFeatures=createFeatures(pChr->pFeatures, pChr, "tRNA:");
		pChr->pFeatures=createFeatures(pChr->pFeatures, pChr, "GIS:");
		pChr->pFeatures=createFeatures(pChr->pFeatures, pChr, "rRNA:");
		pChr->pFeatures=createFeatures(pChr->pFeatures, pChr, "tmRNA:");
		
		printf("The feature num in Chr %s:%d\n", pChr->chr, pChr->numFeatures);

		pChr = pChr->next;
	}
}

void createSequencesFromDB( Chromosome * pGenome) {
	 Chromosome * pChr = pGenome;
	while(pChr != NULL) {
		pChr->pDNA = myMalloc(char, pChr->end - pChr->start + 2);
		
		sprintf(myBuffer, "SELECT sequence FROM sequence WHERE id = %u ORDER BY id, offset", pChr->seqid);
  		executeSQL(myBuffer);
  		mysqlResult = mysql_store_result(mysqlConnection);
		int l = 0;
		char  * ptr = pChr->pDNA;
		while((mysqlRow = mysql_fetch_row(mysqlResult)) != NULL)
    	{
    		l = strlen(mysqlRow[0]);
			memcpy(ptr, mysqlRow[0], l);
			ptr += l;
		}
		
		* ptr = '\0';
		mysql_free_result(mysqlResult);

		pChr = pChr->next;
	}
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void freeVariants( Variant * pVariants) {
	 Variant * p = pVariants, * t;
	while(p != NULL) {
		t = p->next;
		free(p->chr);
		free(p->type);
		free(p->ref);
		free(p->alt);
		free(p);
		p = t;
	}
}

void freeExons( Exon * pExons) {
	 Exon * p = pExons, * t;
	while(p != NULL) {
		t = p->next;
		free(p->type);
		freeVariants(p->pVariants);
		free(p);
		p = t;
	}
}

void freeOther_region( Other_region * pOther_region) {
	 Other_region * p = pOther_region, * t;
	while(p != NULL) {
		t = p->next;
		free(p->name);
		free(p->type);
		freeVariants(p->pVariants);
		free(p);
		p = t;
	}
}

void freeTranscripts( Transcript * pTranscripts) {
	 Transcript * p = pTranscripts, * t;
	while(p != NULL) {
		t = p->next;
		freeExons(p->pExons);
		free(p->name);
		if(p->aas) free(p->aas);
		if(p->cds) free(p->cds);
		free(p);
		p = t;
	}
}

void freeFeatures( Feature * pFeatures) {
	 Feature * p = pFeatures, * t;
	while(p != NULL) {
		t = p->next;
		free(p->name);
		free(p->type);
		freeTranscripts(p->pTranscripts);
		freeOther_region(p->pOther_region);
		free(p);
		p = t;
	}
}

void freeGenome( Chromosome * pGenome) {
	printf(">>> Free the genome...\n");
	 Chromosome * p = pGenome, * t;
	while(p != NULL) {
		t = p->next;
		if(p->pDNA) free(p->pDNA);
		freeVariants(p->pVariants);
		freeFeatures(p->pFeatures);
		free(p);
		p = t;
	}
}

/////////////////////////////////////////////////////////////////////////////
// DNA sequences
void createGenomeDNASequence(char * dnaFileName,  Chromosome * pGenome) {
	FILE * fp = fopen(dnaFileName, "rt");
	if(fp == NULL) return;
	
	 Chromosome * pChr;
	printf("   Creating DNA sequences...\n");
	while(1) {
		fgets(myBuffer, MyBufLen, fp);
		if(feof(fp)) break;
		
		if(myBuffer[0] == '>') {
			// cout << myBuffer;
			// if(myBuffer[1] == 'M') { cout << "   " << myBuffer << endl; fclose(fp); return; }

			// >1 dna:chromosome chromosome:GRCh37:1:1:249250621:1
			// >2 dna:chromosome chromosome:GRCh37:2:1:243199373:1
			// >Y dna:chromosome chromosome:GRCh37:Y:2649521:59034049:1
			char * ptr = strtok(myBuffer, ":");
		
			ptr = strtok(NULL, ":"); ptr = strtok(NULL, ":"); ptr = strtok(NULL, ":");
			
			pChr = ( Chromosome *)getChromosome(ptr, pGenome);
			if(pChr == NULL) {
				printf("      >>> Chromosome %s not found!\n"); 
				fclose(fp); 
				return; 
			}
							
			ptr = strtok(NULL, ":");
			long start = 1;
			start = atoi(ptr);
				
			ptr = strtok(NULL, ":");
			long seqLen = atol(ptr);
			if(myBuffer[1] == 'Y') seqLen += start - 1;
			
			pChr->pDNA = myMalloc(char, seqLen+1);
			if(pChr->pDNA == NULL) 
			{
				printf(">>> No enough memory!\n"); 
				exit(1); 
			}
			char * p = pChr->pDNA;
			while(1) {
				char ch = fgetc(fp);
				if(ch == '\n' || ch == '\r' || ch == '\0') break;
				if(p - pChr->pDNA >= seqLen) break;
				* p ++ = ch;
			}
			// cout << "L = " << p - pChr->pDNA << endl; 
			* p = '\0';
			printf("      Chr %s, seqLen=%d, readLen=%d\n", pChr->chr, seqLen, p - pChr->pDNA);
		}
	} // end of while
	fclose(fp);
}			

////////////////////////////////////////////////////////////////////////////
// for CDS
char * getCDS(Transcript * pT, int strand, const char * pDNA) {
	long len = 0;
	Exon * pE = pT->pExons;
	while(pE) {
		if(!strcmp(pE->type, "CDS:"))
			len += pE->end - pE->start + 1;
		pE = pE->next;
	}
	
	char * pCDS = myMalloc(char, len + 4); // try to include stop codon 
	char * p = pCDS;
	pE = pT->pExons;
	while(pE) {
		if(!strcmp(pE->type, "CDS:")) {
			int l = pE->end - pE->start + 1;
			strncpy(p, pDNA + pE->start - 1, l); 
			p += l;
		}
		pE = pE->next;
	}
	
	* p = '\0';
	if(strand == -1) {
		// printExons(pT->cds);
		getComplementDNA(pCDS);
	}
	
	return pCDS;
}

char * getAAs(const char * cds) {
	int len = strlen(cds);
	
	char * aas = myMalloc(char, len/3 + 1);
	int l = 0, i = 0;
	while(l < len) {
		aas[i] = getAminoAcid(cds + l);
		l += 3;
		++ i;
	}
	aas[len/3] = '\0';
	return aas;
}

void getAllCDSs(Chromosome * pGenome) {
	printf("    Geting all CDSs...\n");
	// ofstream os("CDS.err");
	int numPs = 0;
	int numCDSs = 0;
	Chromosome * pChr = pGenome;
	while(pChr) {
		if(pChr->pDNA == NULL) { // mitochondrion
			pChr = pChr->next; continue; }
		
		int nF = 0;
		 Feature * pF = pChr->pFeatures;
		while(pF) {
			pF->numProteins = 0;
			pF->numCDSs = 0;

			if(!strcmp(pF->type, "gene:")) {
				// cout << "       " << ++ nF << ") " << pF->name << ": ";
				 Transcript * pT = pF->pTranscripts;
				while(pT) {
					// cout << " +";
					pT->aas = NULL;
					{ // protein_coding
						++ pF->numProteins; // protein-coding
						// char * pCDS = getCDS(pT, pF->strand, pChr->pDNA);
						pT->cds = getCDS(pT, pF->strand, pChr->pDNA);
						if(pT->cds) {
							pT->aas = getAAs(pT->cds);
							if(pT->aas) {
								++ pF->numCDSs;
								// cout << '*';
							}	
							// free(pCDS);
						}
					}	
					pT = pT->next;
				}
				// cout << endl;
				numPs += pF->numProteins;
				numCDSs += pF->numCDSs;
			} // end of if
			pF = pF->next;
		}	
		pChr = pChr->next;
	}
	printf("<<< There are %d mRNAs, and %d CDSs annotated from protein_coding genes\n",numPs, numCDSs);
}

void printCDS(FILE * cdsOut, const char * cds, int len) {
	for(int l=0; l<strlen(cds); ++l) {
		fprintf(cdsOut, "%s",cds[l]);
		if((l+1)% len == 0) fprintf(cdsOut, "\n");
	}
	if(strlen(cds) % len) fprintf(cdsOut, "\n");
}			
		
void printAllCDSs(const char * fn, Chromosome * pGenome, int len) {
	printf("    Output all CDSs to file: %s", fn);
	FILE *cdsOut;
	if(!(cdsOut=fopen(fn,"w+")))
    {
        printf("Cannot open the CDS out file !!\n");
        exit(1);
    }	
	Chromosome * pChr = pGenome;
	while(pChr) {
		Feature * pF = pChr->pFeatures;
		while(pF) {
			Transcript * pT = pF->pTranscripts;
			while(pT) {
				if(pT->aas) {
					fprintf(cdsOut, ">%d\t[%s]\t%s\t%c\t%d\t|%d\n", pT->id, pF->name, pChr->chr, (pF->strand == 1 ? '+' : '-'), pF->start, strlen(pT->aas));  
					printCDS(cdsOut, pT->aas, len);
				}
				pT = pT->next;
			}
			pF = pF->next;
		}
		pChr = pChr->next;
	}
	fclose(cdsOut);
}					 
		
//////////////////////////////////////////////////////////////
