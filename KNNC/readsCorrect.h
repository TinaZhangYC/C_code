#ifndef CREATEKMER_H_
#define CREATEKMER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h>
#include "hashTable.h"

#define MAX_KMER_LEN 100
#define MAX_REFER_LEN  0x7ffffffff
#define MAX_LEN        0xffffff
#define MAX_SEQ_NUM    0xffffffff
#define NAME_LEN	   500
#define HASH_TABLE_MAX_SIZE  0x7ffffff


typedef struct {
	unsigned int k;  //the length of k-mer
	unsigned int l; //lfreq 
	unsigned int h; //hfreq weight
	unsigned int a; //acceptThre weight
	unsigned int d; //discardThre weight
	unsigned int f; //fix weight 
	unsigned int r;  //define the refiltering level
}genome_net_opt;

/*typedef struct kmerNode{
	long long int Id;
	int counter;
	long long int location[MAX_LOCATION_HITS];
	char seq[MAX_KMER_LEN];
	struct kmerNode*  Next;
}*kmerList;*/

#ifdef __cplusplus
 extern "C"{
 #endif
 	long long int power(int a,int n);
	long long int  getId(char *seq);
	char*  substr(char*  str,long int   start,long int lenght);
 
 	//creat kmerNet
	int createKmerNet(char * ReferStr,unsigned int k,char *revStr);

	//correction process
	  //filter
	unsigned long long int insertArray(unsigned long int *kmerFreq,unsigned long int freq,unsigned long long  int n);
	int cmp ( const void *a , const void *b);
	void countFivenum(unsigned long int *kmerFreq,unsigned long int *fiveNum);
	int prefilter(char * ReferStr,unsigned int k,unsigned long int hfreqThre,unsigned long int lfreqThre,unsigned int acceptThre,unsigned int discardThre);
	int filter(char * ReferStr,unsigned int k,char * revStr,unsigned long int *kmerFreq,unsigned long int *fiveNum);
	 //fix
	unsigned int pcmp(char *p1,char *p2,unsigned int k);
	unsigned int fixPos(char* p1,char * p2,unsigned int k);
	int correct(char *ReferStr,unsigned int k,unsigned int fixlevel,unsigned long int lfreqThre,unsigned long int hfreqThre );
	 //refilter
	int refilter(char * ReferStr,char *oriReferStr,unsigned int k,unsigned int refilterLevel,unsigned long int lfreqThre);
	//output
     void printKmerNet(char * ReferStr,unsigned int k,FILE *kmerNetOut);
	 void printAcceptReads(char * ReferStr,unsigned int k,FILE *accReadsOut);
	 void printsubAcceptReads(char * oriReferStr,unsigned int k,FILE *subaccReadsOut);
	 void printDiscardReads(char * oriReferStr,unsigned int k,FILE *disReadsOut);
	 void printprecorReads(char * ReferStr,unsigned int k,FILE *precorReadsOut);
	 void printCorrectReads(char * ReferStr,unsigned int k,FILE *corReadsOut);

   //Initiation
   	unsigned long long int Initiation(char * oriReferStr,unsigned int k,char *readsLibName);

   //Scheduling
 int schedule(unsigned int k,char *readsLibName,unsigned int lfreq,unsigned int hfreq,unsigned int accThre,unsigned int disThre,unsigned int fixscale,unsigned int refilterscale);

 #ifdef __cplusplus
}
#endif

#endif

