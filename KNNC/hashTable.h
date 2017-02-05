#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define HASH_TABLE_MAX_SIZE  0x7ffffff
#define MAX_TIMES	  1000
#define MAX_KEY_LEN   100
#define MAX_NEIGH_SIZE 75

typedef struct readNode
{
	unsigned long int readId;
	struct readNode *rNext;
}*ReadLine;

typedef struct neighborNode
{
	char * neighbor; //string array
	unsigned long int frequency;
	struct  neighborNode *nNext;  //point to next neighbor
}*Neighbor;

typedef struct locateNode
{
	unsigned long long int location; //the position on the reference
	struct  locateNode *lNext;  //point to next location on the reference 
}*Location;

//define the hashNode type
typedef struct HashNode
{
	//int nHashA;
	//int nHashB;
	char *pKey;  //store the kmer on the genome
	Neighbor neighbors;
	unsigned int neighnum;
	Location  locations;
	unsigned long int frequency;		//value3: record the times of this kmer appearing in the reference
	struct HashNode  * pNext; //handle hash confiliction: point to next node 
}*Hashtable;

#ifdef __cplusplus
	extern "C"{
	#endif
	void hash_table_init();//Initiate hash table
	void hash_table_insert(char *referStr,char* pkey,char * revStr,unsigned int k,unsigned long long int position); //insert skey-nvalue
	void hash_table_remove(char *referStr,char* pkey,char *revStr,unsigned int k); //delete key-value
	long long int hash_table_hash_str(char* pkey,unsigned int k);
	unsigned long int RSHash(char * str,unsigned int k);
	void prepareCryptTable();
	unsigned long long int HashString(char *lpszFileName,unsigned int k, unsigned long dwHashType );
	int pstrcmp(char*pkey1,char *pkey2,unsigned int k);
	void revSeq(char * str,char *revStr,unsigned int k);
	Hashtable hash_table_lookup(char* pkey,unsigned int k) ;//look up the node of skey and return the hashnode pointer or NULL
	void hash_table_print(char *ReferStr);
	void hash_table_release();//free the memory of hash table

	#ifdef __cplusplus
}
#endif

#endif


