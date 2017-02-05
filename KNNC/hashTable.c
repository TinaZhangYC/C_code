#include "hashTable.h"


Hashtable hashTable[HASH_TABLE_MAX_SIZE]; //array of hash table 
Hashtable orihashTable[HASH_TABLE_MAX_SIZE];
unsigned long long int hash_table_size; //the number of elements in the hash table
//unsigned long int cryptTable[0x500];

//initialize hash table
void hash_table_init()
{
	hash_table_size = 0;
	//memset(hashTable, 0, sizeof(Hashtable) * HASH_TABLE_MAX_SIZE);
}


//string hash function
long long int hash_table_hash_str(char* pkey,unsigned int k)
{
	char *p = pkey;
	long long int h = *p;
	unsigned i;
	if(h)
	{
		for(i=0;i<k;i++,p++)
		 //h = (h << 2) - h + *p;
			h+=*p;
	}
	return h;
}

unsigned long int RSHash(char * str,unsigned int k)
{
	char *pstr;
	pstr=str;
	unsigned long int b = 378551 ;
	unsigned long int a = 63689 ;
	unsigned int i;
	unsigned long long  int hash = 0 ;
	for(i=0;i<k;i++)
	{
		hash = hash * a + *pstr;
		a *= b;
		pstr++;
	}
	return (hash & 0x7FFFFFFF);
}

//create cryptTable
/*void prepareCryptTable()
{
	unsigned long int seed = 0x00100001, index1 = 0, index2 = 0, i;
	for( index1 = 0; index1 < 0x100; index1++ )
	{
		for( index2 = index1, i = 0; i < 5; i++, index2 += 0x100 )
		{
			unsigned long int temp1, temp2;
			seed = (seed * 125 + 3) % 0x2AAAAB;
			temp1 = (seed & 0xFFFF) << 0x10;													   
			seed = (seed * 125 + 3) % 0x2AAAAB;
			temp2 = (seed & 0xFFFF);
			cryptTable[index2] = ( temp1 | temp2 );
		}
	}
} 

//hash function
unsigned long long int HashString(char *lpszFileName, unsigned int k,unsigned long dwHashType )
{
    const  char *key=lpszFileName;
	unsigned long long int seed1 = 0x7FED7FED;
	unsigned long long int seed2 = 0xEEEEEEEE;
	int ch;
	unsigned int i;
	//unsigned int k;
 
	for(i=0;i<k;i++)
	{
		ch = toupper(*key++);
		seed1 = cryptTable[(dwHashType << 8) + ch] ^ (seed1 + seed2);
		seed2 = ch + seed1 + seed2 + (seed2 << 5) + 3;
	}
	return seed1;
}*/

int pstrcmp(char*pkey1,char*pkey2,unsigned int k)
{
	 char *p1,*p2;
	 p1=pkey1;
	 p2=pkey2;
	 unsigned int i;
	 for(i=0;i<k;i++)
	 {
	 	if(*p1!=*p2)
		{
		  return 1;
		}
		else
		{
			p1++;
			p2++;
		}
	 }
	  return 0;
}

void revSeq(char * str,char *revStr,unsigned int k)
 {
     unsigned int i=0;
     char base;
     char *pstr=str+k-1;
     char *prStr=revStr;
	 //memset(revStr,'\0',(k+1));
     for(i=0;i<k;i++)
     {
         base=*pstr;
		 pstr--;
         switch(base)
         {
             case 'A':
                     (*prStr)='T';
                     prStr++;
                     break;
             case 'T':
                     (*prStr)='A';
					 prStr++;
                     break;
             case 'C':
                     (*prStr)='G';
                     prStr++;
                     break;
             case 'G':
                     (*prStr)='C';
                     prStr++;
					 break;
         };
     }//for :i
}

//insert key-value into hash table
void hash_table_insert(char * ReferStr,char* pkey,char *revStr,unsigned int k,unsigned long long int position)
{
	//unsigned int i;
	/*if(hash_table_size >= HASH_TABLE_MAX_SIZE)
	{
		printf("out of hash table memory!\n");
		return;
	}*/
	unsigned long  int  nHash =RSHash(pkey,k);
    //unsigned long long int  nHashA = HashString(skey, HASH_A );
    //unsigned long long int  nHashB = HashString(skey, HASH_B );
    unsigned long  int  nHashStart = nHash % HASH_TABLE_MAX_SIZE;
    unsigned long  int  pos = nHashStart;
	

    memset(revStr,'\0',(k+1));
	Hashtable pHead = hashTable[pos];
	Hashtable tmpHead;
	tmpHead=pHead;
	while(pHead)
	{
		/*if (pHead->nHashA == nHashA && pHead->nHashB == nHashB )
			{
				pHead->location[pHead->times++]=position;
				return ;
			 }*/
			 //posKey=ReferStr+pHead->posValue;
	
			//printf("posValue:%d",pHead->pKey-ReferStr);
			
		if(pstrcmp(pHead->pKey,pkey,k)== 0)
		{
			//pHead->location[pHead->times++]=position;
		   Location lNode;
		   while(!(lNode=calloc(1,sizeof(struct locateNode))))
            {
		         printf("malloc memory error for location node! try again!!\n");
	        }
			//memset(lNode,0,sizeof(struct locateNode));
			lNode->location=position;
			lNode->lNext=pHead->locations->lNext;
			pHead->locations->lNext=lNode;
			pHead->frequency++;
			  return ;
		}//if:pstrcmp
		tmpHead=pHead;
		pHead = pHead->pNext;
	}
    /////////////////////////////////////////////////
    //reverse complement sequence
      revSeq(pkey,revStr,k);
      unsigned long int prePos=pos;
      nHash=RSHash(revStr,k);
      nHashStart=nHash % HASH_TABLE_MAX_SIZE;
      pos=nHashStart;
      Hashtable pHeadR=hashTable[pos];
	  while(pHeadR)
      {
	       if(pstrcmp(pHeadR->pKey,revStr,k)== 0)
            {
				//pHead->location[pHead->times++]=position;
             	//printf("again!!\n");
				Location lRNode;
				while(!(lRNode=calloc(1,sizeof(struct locateNode))))
		        {
	                 printf("malloc memory error for location node! try again!!\n");
	            }
		        //memset(lRNode,0,sizeof(struct locateNode));
		        lRNode->location=position;
	            lRNode->lNext=pHeadR->locations->lNext;
			    pHeadR->locations->lNext=lRNode;
             	pHeadR->frequency++;
               	return ;
         	}
         pHeadR = pHeadR->pNext;
      }
	//////////////////////////////////////////////////
	pos=prePos;
	Hashtable pNewNode;
	while(!(pNewNode=calloc(1,sizeof(struct HashNode))))
    {
          printf("malloc memory error for hash node! try again!!\n");                 
	}
	//memset(pNewNode,0,sizeof(struct HashNode));
	Location lNode;
    //printf("Insert a node !\n");
    while(!(lNode=calloc(1,sizeof(struct locateNode))))
    {
	     printf("malloc memory error for location node! try again!!\n");
    }
	//memset(lNode,0,sizeof(struct locateNode));
    lNode->location=position;
	lNode->lNext=NULL;
    //locateNode->lNext=pNewNode->locations->lNext;
	pNewNode->locations=lNode;
	pNewNode->pKey =pkey; 
	pNewNode->pNext=NULL;
	pNewNode->frequency=1;
	pNewNode->neighbors=NULL;
	pNewNode->neighnum=0;
	++hash_table_size;
	//pNewNode->location[pNewNode->times++]=position;
	 //pNewNode->nHashA=nHashA;
	 //pNewNode->nHashB=nHashB;

	//printf("\npkey=%lld\n",pNewNode->pKey);
	
	if(hashTable[pos]==NULL)
	   {
			//pNewNode->pNext = hashTable[pos];
			hashTable[pos] = pNewNode;
	   }
	else
	   {
	   	    tmpHead->pNext=pNewNode;
	   }
}


//remove key-value frome the hash table
void hash_table_remove(char * ReferStr,char* pkey,char *revStr,unsigned int k)
{
	 unsigned long int start = RSHash(pkey,k);
	 unsigned long int npos=start % HASH_TABLE_MAX_SIZE;
	 unsigned long int pos=npos;

	 Hashtable pHead = hashTable[pos];
	 Hashtable pLast = NULL;
	 Hashtable pRemove = NULL;
	 //int i;
	 while(pHead)
	 {
		if(pstrcmp(pkey,pHead->pKey,k) == 0)
		{
			pRemove = pHead;
		    //printf("Find the kmer!!\n");
			break;
		}
		pLast = pHead;
		pHead = pHead->pNext;
	 }
	 if(pRemove==NULL)
	 {
		revSeq(pkey,revStr,k);
		start=RSHash(revStr,k);
		pos=start % HASH_TABLE_MAX_SIZE;
		pHead=hashTable[pos];
		while(pHead)
		{
			if(pstrcmp(revStr,pHead->pKey,k)==0)
			{
				pRemove=pHead;
				//printf("find the kmer!\n");
				break;
			}
			pLast=pHead;
			pHead=pHead->pNext;
	    }
	  }
	  if(pRemove)
	  {
		if(pLast)
			pLast->pNext = pRemove->pNext;
		else
			hashTable[pos] = pRemove->pNext;
		 //printf("Remove pkey:%lld\n",pRemove->pKey-ReferStr);
		 free(pRemove->locations);
		 free(pRemove);
		 hash_table_size--;
	  }
}

//lookup a key in the hash table
Hashtable  hash_table_lookup(char* pkey,unsigned int k)
{
	unsigned long int npos = RSHash(pkey,k);
	unsigned long int mpos=npos % HASH_TABLE_MAX_SIZE;
	unsigned long int pos=mpos;
	if(hashTable[pos])
	{
			Hashtable pHead = hashTable[pos];
			while(pHead)
			{
				if(pstrcmp(pHead->pKey,pkey,k) == 0)
				{
					return pHead;
				}
				pHead = pHead->pNext;
			}	
	}
	return NULL;
}

//print the content in the hash table
void hash_table_print(char *ReferStr)
{
	printf("===========content of hash table=================\n");
	unsigned long long int i;
	for(i = 0; i < HASH_TABLE_MAX_SIZE; ++i)
		if(hashTable[i])
		{
			Hashtable pHead = hashTable[i];
			while(pHead)
			{
				printf("\n%lld=>\nposValue=%ld:frequency=%ld\n",i,pHead->pKey-ReferStr,pHead->frequency);
				pHead = pHead->pNext;
			}
			printf("\n");
		}
}

//free the memory of the hash table
void hash_table_release()
{
	long long int i;
	for(i = 0; i < HASH_TABLE_MAX_SIZE; ++i)
	{
		if(hashTable[i])
		{
			Hashtable pHead = hashTable[i];
			while(pHead)
			{
				Hashtable pTemp = pHead;
				pHead = pHead->pNext;
				if(pTemp)
				{
					//free(pTemp->neighbors);
					free(pTemp->locations);
					free(pTemp);
				}

			}
		}
		if(orihashTable[i])
		{
			Hashtable pHead = orihashTable[i];
			while(pHead)
			{
				Hashtable pTemp = pHead;
				pHead = pHead->pNext;
				if(pTemp)
				{
					//free(pTemp->neighbors);
					//free(pTemp->locations);
					free(pTemp);
				}

			}
		}
	}
}

/* ===============================hash table end=========================*/

/* ============================test function ============================*/
/*#define MAX_STR_LEN 20
#define MIN_STR_LEN 10
void rand_str(char r[])
{
	int i;
	int len = MIN_STR_LEN + rand() % (MAX_STR_LEN - MIN_STR_LEN);
	for(i = 0; i < len - 1; ++i)
		r[i] = 'a' + rand() % ( 'z' - 'a');
	r[len - 1] = '\0';
}*/

/*int main(int argc, char** argv)
{
	//srand(time(NULL));
	printf("Beginning....\n");
	hash_table_init();
	printf("insert testing.........\n");
	int n = 10;
	const char *key1 = "aaammdgggggcccccccgggggg";
	const char *key2 = "xzzyymmmmmmmmsssssaaaaaa";
	const char *key3 = "cdcdedaddsfsffsfdgfdgdfb";

	prepareCryptTable();
	
	hash_table_insert(key1,1);
	hash_table_insert(key1,7);
	hash_table_insert(key2,13);
	hash_table_insert(key2,19);
	hash_table_insert(key3,25);
	
	char str[MAX_STR_LEN + 1];
	while(n--)
	{
		rand_str(str);
		hash_table_insert(str,n);
	}*/
   

	//hash_table_print();

	/*printf("\nlookup testing..........\n");
	Hashtable pNode = hash_table_lookup(key1);
	printf("\nlookup result:%d\n", pNode->idValue);
	pNode = hash_table_lookup(key2);
	printf("\nlookup result:%d\n", pNode->idValue);

	printf("\nremove testing..........\n");
	printf("\nbefore remove %s:\n", key3);
	hash_table_print();
	hash_table_remove(key3);
	printf("\nafter remove:\n");
	hash_table_print();
	//hash_table_release();
	
	//system("pause");
	return 0;
}*/
