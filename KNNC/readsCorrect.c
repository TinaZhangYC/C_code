#include "readsCorrect.h"

Hashtable hashTable[HASH_TABLE_MAX_SIZE]; //array of hash table 
//Hashtable orihashTable[HASH_TABLE_MAX_SIZE];
unsigned long long int hash_table_size; //the number of elements in the hash table
//unsigned long long int cryptTable[0x500];
unsigned int readLen;
unsigned long long int readsnumber;
char nucleotide[]="ACGT";
ReadLine acceptReads;
ReadLine discardReads;
ReadLine corrigibleReads;
ReadLine prediscardReads;
ReadLine correctedReads;
ReadLine precorrectedReads;
ReadLine subacceptReads;

long long int power(int a,int n)
{
	long long int  result=1;
	int i;
	for(i=1;i<=n;i++)
	 result*=a;
	
	return result;
}

long long int  getId(char *seq)
{
	long long int k=strlen(seq);
	long long int i;
	long long int id=1;
	for(i=0;i<k;i++)
	{
		switch(seq[i])
		{
			case 'A':
					id+=0*power(4,k-i-1);
					break;
			case 'C':
					id+=1*power(4,k-i-1);
					break;
			case 'G':
					id+=2*power(4,k-i-1);
					break;
			case 'T':
					id+=3*power(4,k-i-1);
					break;
			default:
					break;
		};
	}
	 	return id;
}

char*  substr(char* str,long int start,long  int lenght)
{
	if(!str)	return   NULL;
	if(start>=strlen(str))	return   NULL;
	if(lenght==0 || lenght>(strlen(str)-start))
		lenght=strlen(str)-start;
	char* tmp=(char*)malloc((lenght+1)*sizeof(char));
	memset(tmp,0,lenght+1);
	memcpy(tmp,str+start,lenght);
	tmp[lenght]='\0';
	return   tmp;
}


/**************************************************/
/********** The kmerNet-creation module ***********/
/**************************************************/
int createKmerNet(char * oriReferStr,unsigned int k,char *revStr) 
{
	unsigned long long  int i=0,j=0;
	char *pseq;
	unsigned int p=0;  //number of neighbors
	unsigned int m;
    Hashtable locate;
	char *position;
	char origine;
	Hashtable pHead;
//	unsigned int lookflag=0; //0: hashTable 1: orihashTable
	//char *revStr=(char *)malloc((k+1)*sizeof(char));
	unsigned long long int site=0;

	/*************** Initiation ******************/
	hash_table_size=0;
	printf("Scanning.....\n");
	//printf("The len of refer:%lld\n",strlen(ReferStr));
	
	//fprintf(kmerNetOut,">The k-value:%d\n",k);
	//fprintf(kmerNetOut,">Theoretical Max theoryKmer-num:%lld\n",(readLen-k+1)*readnumber);
    //prepareCryptTable();
	pseq=oriReferStr;
    printf("The reads length:%d\n",readLen);
	printf("The reads numbers:%lld\n",readsnumber);
	//fprintf(kmerNetOut,">The read-length:%d\n",readLen);
	//fprintf(kmerNetOut,">The reads-num:%lld\n",readnumber)

	/**************construct hash table ************/
	for(i=0;i<readsnumber;i++)
	{
		for(j=0;j<=readLen-k;j++)
	    {
			site=pseq-oriReferStr; //relevate path
		    hash_table_insert(oriReferStr,pseq,revStr,k,site);
		    pseq++;
        }
		pseq--;
		pseq+=k;
	 }
   //printf("Handle %lld times!\n",N);
   printf("There are %lld kmers!\n",hash_table_size);
   //fprintf(kmerNetOut,">The realKmer-num:%lld\n",hash_table_size);

   /**************construct genotype network*****************/
   unsigned long long int number=0;

   //unsigned long long int count=0;

   for(i=0;i<HASH_TABLE_MAX_SIZE;i++)
	 {
	 	pHead=hashTable[i];
	 	while(pHead)
		{
			//find neighbors
			p=0;
			position=pHead->pKey;//absolute path
			for(j=0;j<k;j++) 
			{
				origine=*position;
				for(m=0;m<4;m++)
				{
					if(nucleotide[m]==origine)
						continue;
					else
					{
						(*position)=nucleotide[m];
						revSeq(pHead->pKey,revStr,k);//the reverse complement of the kmer
						if((locate=hash_table_lookup(pHead->pKey,k))!=NULL)
						    {
						     	Neighbor nNode;
							    while(!(nNode=calloc(1,sizeof(struct neighborNode))))
						        {
						           printf("malloc memory error for neighbor node! try again!!\n");
							    }
							    //memset(nNode,0,sizeof(struct neighborNode));
								nNode->neighbor=locate->pKey;
								nNode->frequency=locate->frequency;
								nNode->nNext=NULL;
								if(p==0)
								{
									pHead->neighbors=nNode;
								}
								else
								{
							    	nNode->nNext=pHead->neighbors->nNext;
							 		pHead->neighbors->nNext=nNode;
								}
								p++;
						    }//if
						  else if ((locate=hash_table_lookup(revStr,k))!=NULL)
						    {
						     	Neighbor nNode;
							    while(!(nNode=calloc(1,sizeof(struct neighborNode))))
						        {
						           printf("malloc memory error for neighbor node! try again!!\n");
							    }
							    //memset(nNode,0,sizeof(struct neighborNode));
								nNode->neighbor=locate->pKey;
								nNode->frequency=locate->frequency;
								nNode->nNext=NULL;
								if(p==0)
								{
									pHead->neighbors=nNode;
								}
								else
								{
							    	nNode->nNext=pHead->neighbors->nNext;
							 		pHead->neighbors->nNext=nNode;
								}
								p++;
						    }//else if: locate
					}//else nucleotide[m]!=origine
				}//for:m
				 (*position)=origine;
				   position++;
			 }//for:j<k
			//printf("kmer:%s\n",pHead->pKey);
			pHead->neighnum=p;
			number++;
			pHead=pHead->pNext;
		 }//while (pHead)
    }//for:i
    printf("The calculate number:%lld\n",number);

	//test
	/*Neighbor neigh;
	char *tmpneigh;
	Hashtable neighNode;
	number=0;
	for(i=0;i<HASH_TABLE_MAX_SIZE;i++)
	{
		pHead=hashTable[i];
		while(pHead)
		{
			number++;
			neigh=pHead->neighbors;
			while(neigh)
			{
				tmpneigh=neigh->neighbor;
				neighNode=hash_table_lookup(ReferStr,tmpneigh,k);
				if(neighNode==NULL)
				{
					revSeq(tmpneigh,revStr,k);
					neighNode=hash_table_lookup(ReferStr,revStr,k);
				}
				if(neighNode==NULL)
				{
					printf("The %lld kmer is error!\n",number);
					printf("Neighbors error in create kmerNet......\n");
				}
				neigh=neigh->nNext;
			}
			pHead=pHead->pNext;
		}
	}

   printf("The real kmer number:%lld\n",number);*/
	//hash_table_release();
	return 0;
}

/*************************************************************/
/**************** The filtering module **********************/
/*************************************************************/
unsigned long long int insertArray(unsigned long int *kmerFreq,unsigned long int freq,unsigned long long  int n)
{
	unsigned long long int i=0;
	for(i=0;i<n;i++)
	{
		if(kmerFreq[i]==freq)
		{
		    return i;
		}
	}
	kmerFreq[n]=freq;
	return n;
}
int cmp ( const void *a , const void *b)
{
     return *( int *)a - *(int *)b;
}

void countFivenum(unsigned long int *kmerFreq,unsigned long int *fiveNum)
{	
	unsigned long long int n,tmpn;
	//unsigned long long int *tmpN;
	unsigned long long int i=0;
	unsigned long int Q1,Median,Q3,IQR;
	Hashtable pHead;
	n=0;
	for(i=0;i<HASH_TABLE_MAX_SIZE;i++)
	{
		pHead=hashTable[i];
		while(pHead)
		{
			tmpn=insertArray(kmerFreq,pHead->frequency,n);
			if(tmpn==n)
			{
				n++;
			}
			pHead=pHead->pNext;
		}
	}
	qsort(kmerFreq,n,sizeof(kmerFreq[0]),cmp); //sort kmer freq
	
	//
	/*printf("The soretd freqs:\n");
	for(i=0;i<n;i++)
	{
		printf("%ld\t",kmerFreq[i]);
	}
	printf("\n");*/
	//

	Q1=kmerFreq[(n+1)/4];
	Median=kmerFreq[(n+1)/2];
	Q3=kmerFreq[3*((n+1)/4)];
	IQR=Q3-Q1;
	fiveNum[0]=kmerFreq[0]; //min
	fiveNum[1]=Q1;
	fiveNum[2]=Median;
	fiveNum[3]=Q3;
	fiveNum[4]=kmerFreq[n-1];//max
	fiveNum[5]=IQR;

	printf("Print the fivenum:\n");
	printf("Min:%ld\tQ1:%ld\tMedian:%ld\tQ3:%ld\tMax:%ld\n",fiveNum[0],fiveNum[1],fiveNum[2],fiveNum[3],fiveNum[4]);

}
int prefilter(char * oriReferStr,unsigned int k,unsigned long int hfreqThre,unsigned long  int lfreqThre,unsigned int acceptThre,unsigned int discardThre)
{
	unsigned long long int i=0;
	unsigned long int j=0;
	unsigned long long int accnum=0,cornum=0,disnum=0;
	char *pseq;
	Hashtable kmerNode;
	unsigned long long int atmpnum=0,dtmpnum=0;
	unsigned long int freq=0;
	char *revStr;
	//unsigned int lookflag=0; //0: hashTable 1:orihashTable
	while(!(revStr=(char *)malloc((k+1)*sizeof(char))))
	{
		printf("malloc memory error for revStr-array!Try again!\n");
	}
	memset(revStr,'\0',(k+1));
	//
	printf("start flitering.......\n");
	printf("readLen:%d\nreadsnumber:%lld\n",readLen,readsnumber);
	//
	pseq=oriReferStr;
	for(i=0;i<readsnumber;i++)
	{
		atmpnum=0;
		dtmpnum=0;
		for(j=0;j<=readLen-k;j++)
		{
			kmerNode=hash_table_lookup(pseq,k);
			if(kmerNode==NULL)
			{
				revSeq(pseq,revStr,k);//the reverse complement of the kmer
				kmerNode=hash_table_lookup(revStr,k);
			}
			if(kmerNode==NULL)
			{
				printf("Cannot find the kmer in filtering!\n");
				pseq++;
				continue;
			}
			freq=kmerNode->frequency;
			if(freq>=hfreqThre)
				atmpnum++;
			if(freq<=lfreqThre)
				dtmpnum++;
			pseq++;
		}//for: j<=readLen-k
		//accept
		if(atmpnum>=acceptThre && dtmpnum<discardThre)
		{
			ReadLine accread;
			while(!(accread=calloc(1,sizeof(struct readNode))))
			{
				printf("malloc memory error for accread node!Try again!\n");
			}
			accread->readId=i;
			accread->rNext=NULL;
			if(acceptReads==NULL)
			{
				acceptReads=accread;
			}
			else
			{
				accread->rNext=acceptReads->rNext;
				acceptReads->rNext=accread;
			}

			//
			//printf("acceptReads Id:%lld\n",accread->readId);
			//
			pseq--;
			pseq+=k;
			accnum++;
			continue; //back to i<readsnumber
		}//if: atmpnum
		 
		 //discard
		if(atmpnum<acceptThre && dtmpnum>=discardThre)
		{
			ReadLine disread;
			while(!(disread=calloc(1,sizeof(struct readNode))))
			{
				printf("malloc memory error for diread node!Try again!\n");
			}
			disread->readId=i;
			disread->rNext=NULL;
			if(discardReads==NULL)
			{
				discardReads=disread;
			}
			else
			{
				disread->rNext=discardReads->rNext;
				discardReads->rNext=disread;
			}
			
			//
			//printf("discardReads Id:%lld\n",disread->readId);
			//
			pseq--;
			pseq+=k;
			disnum++;
			continue; //back to i<readsnumber
		 }//if: dtmpnum  
	   
		   //corrigible
		ReadLine corread;
		while(!(corread=calloc(1,sizeof(struct readNode))))
		{
			printf("malloc memory error for corread node! Try again!\n");
		}
		corread->readId=i;
		corread->rNext=NULL;
		if(corrigibleReads==NULL)
		{
			corrigibleReads=corread;
		}
		else
		{
			corread->rNext=corrigibleReads->rNext;
			corrigibleReads->rNext=corread;
		}	
		pseq--;
		pseq+=k;
		cornum++;
		continue; //back to i<readsnumber
	}//for i<readsnumber
	printf("The infintely accepted reads number:%lld\n",accnum);
	printf("The infintely discarded reads number:%lld\n",disnum);
	printf("The corrigible reads number:%lld\n",cornum);
	free(revStr);
	return 0;
}//prefilter

int filter(char * oriReferStr,unsigned int k,char * revStr,unsigned long int *kmerFreq,unsigned long int *fiveNum)
{
	unsigned long int hfreqThre,lfreqThre;
	unsigned int acceptThre,discardThre;

	//Calculate the fivenum array
	countFivenum(kmerFreq,fiveNum);

	//filter and classify reads
	hfreqThre=fiveNum[1]; //Q1
	lfreqThre=1;
	acceptThre=readLen-k+1;
	discardThre=2;
	prefilter(oriReferStr,k,hfreqThre,lfreqThre,acceptThre,discardThre);
	return 0;
}

/*******************************************************/
/****************** The fixing module ******************/
/*******************************************************/
unsigned int pcmp(char *p1,char *p2,unsigned int k)
{
	unsigned int i,diff; 
	char *tmp1,*tmp2;
	tmp1=p1;
	tmp2=p2;
	diff=0;
	for(i=0;i<k;i++)
	{
		if(*tmp1!=*tmp2)
		{
			diff++;
		}
		tmp1++;
		tmp2++;
	}
	return diff;
}

unsigned int fixPos(char* p1,char * p2,unsigned int k)
{
	unsigned int i=0;
	char *tmp1,*tmp2;
	tmp1=p1;
	tmp2=p2;
	for(i=0;i<k;i++)
	{
		if(*tmp1!=*tmp2)
		  return i;
		tmp1++;
		tmp2++;
	}
	return k;
}

int correct(char *ReferStr,unsigned int k,unsigned int fixlevel,unsigned long int lfreqThre,unsigned long int hfreqThre)
{
	//unsigned long int *postFreq;
	unsigned int diffnum=0,diffpos=0;  //sign the different poses and position of two kmers
	unsigned long int lfreqnum=0,hfreqnum=0; //record the neigh-freqs after correct to determin whether to correct
	unsigned long long int i=0,j=0;
	unsigned long int maxNfreq=0,freq=0;
	unsigned long int id=0;
	unsigned long long int read=0;//realtive path
	unsigned int corrflag=0; //sign the correction of read: 0: no correction,1: correction
	char * localpos; //absolute path

	/*char *oread; 
	while(!(oread=(char *)malloc((readLen+1)*sizeof(char))))
	{
		printf("malloc memory error for origine read-array!Try again!\n");
	}
	memset(oread,'\0',(readLen+1));*/
	
	char *orikmer;
	while(!(orikmer=(char *)malloc((k+1)*sizeof(char))))
	{
		printf("malloc memory error for orikmer-array! Try again!\n");
	}
	memset(orikmer,'\0',(k+1));
	/*char *tmporikmer;
	while(!(tmporikmer=(char *)malloc((k+1)*sizeof(char))))
	{
		printf("malloc memory error for tmporikmer-array! Try again!\n");
	}
	memset(tmporikmer,'\0',(k+1));*/
	char *revStr;
	while(!(revStr=(char *)malloc((k+1)*sizeof(char))))
	{
		printf("malloc memory error for revStr-array! Try again!\n");
	}
	memset(revStr,'\0',(k+1));
  /*char *neighstr;
	while(!(neighstr=(char *)malloc((k+1)*sizeof(char))))
	{
		printf("malloc memory error for neighstr-array!Try again!\n");
	}
	memset(neighstr,'\0',(k+1));*/

	char *prefer,*neighpos; //absolute path
	//char *porireads;
	Hashtable tmpKmer=NULL;
	char *tmp; //absolute path
	char *pmemcpy;
	unsigned int startpos=0,endpos=0;
	unsigned long long int cornum=0;
	ReadLine pread=NULL;
	Hashtable kmer=NULL;
	Neighbor neigh=NULL;
	pread=corrigibleReads;
	while(pread)
	{
		//unsigned long long int m;
		prefer=ReferStr;
		id=pread->readId;
		read=id*readLen;
		prefer+=read;//absolute path
		
		corrflag=0;
		lfreqnum=0;
		hfreqnum=0;
		//localpos=prefer;
		for(i=0;i<=readLen-k;i++) //i: relative path
		{
			localpos=prefer+i; //absolute path
			kmer=hash_table_lookup(localpos,k); //changed kmer network
			if(kmer==NULL)
			{
				revSeq(localpos,revStr,k);//the reverse completement of the kmer
				kmer=hash_table_lookup(revStr,k); //origine kmer network
			}
			if(kmer==NULL)
			{
				printf("Cannot find the kmer in correcting:%ld,maybe it is silenced!\n",localpos-ReferStr);
				//localpos++;
				lfreqnum=readLen-k+1;
				ReadLine disread;
				while(!(disread=calloc(1,sizeof(struct readNode))))
				{
					printf("malloc memory error for diread node in the process of correction!Try again!\n");
				}
				disread->readId=pread->readId;
				disread->rNext=NULL;
				if(prediscardReads==NULL)
				{
					prediscardReads=disread;
				}
				else
				{
					disread->rNext=prediscardReads->rNext;
					prediscardReads->rNext=disread;
				}
				pread=pread->rNext;
				break; //go to while: pread*/
				//lfreqnum++;
				//continue; //back to i<=readLen-k
			 }
			freq=kmer->frequency;
			if(freq<=lfreqThre)  //should be corrected
			{
				//unsigned long long int num; //relative path
				//Hashtable tmpneigh;
				neigh=kmer->neighbors;
				if(neigh==NULL)  //acnode,cannot be corrected!
				{
					//printf("Cannot find neighbor!It's an acnode!\n");
					//localpos++;
					continue; //back to for: i<=readLen-k
				}
				maxNfreq=neigh->frequency;
				neighpos=neigh->neighbor;
				neigh=neigh->nNext;
				while(neigh)
				{
					if(maxNfreq<neigh->frequency)
					{
						neighpos=neigh->neighbor;
						maxNfreq=neigh->frequency;
					}
					neigh=neigh->nNext;
				}//while: neigh				
			/*	tmpneigh=hash_table_lookup(neighpos,k);  //origine kmer network
				if(tmpneigh==NULL)
				{
					revSeq(neighpos,revStr,k);
					tmpneigh=hash_table_lookup(revStr,k);  //origine kmer network
				}
				if(tmpneigh==NULL)
				{
					printf("This neighbor doesn't exist!,The neighbor kmer is corrected!\n");
					//localpos++;
					continue; //back to for : i<=readLen-k
				}
				maxNfreq=tmpneigh->frequency;  */
				if(maxNfreq<=lfreqThre) //cannot find proper neighbor,cannot correct the kmer
				 {	
				 	printf("Cannot find highfreq neighbor!\n");
					//localpos++;
					continue; //back to for: i<=readLen-k
				 }
				pmemcpy=memcpy(orikmer,localpos,k); //store the origine kmer
				if(pmemcpy==NULL)
				{
					printf("Localpos to orikmer memory copy failed in the correcting.....\n");
					continue; //back to for: i<=readLen-k
				}
				//memcpy(neighstr,tmpneigh->pKey,k);
				
				//change the origine refer!!
				diffnum=0;
				diffnum=pcmp(localpos,neighpos,k);//whether to use neighbor or its' reverse completement
				if(diffnum!=1) //neighbor reverse completement
				{
					//unsigned int tmpdiffnum=0;
					revSeq(neighpos,revStr,k);
				    /*tmpdiffnum=pcmp(revStr,k); //get rid of fake neighbor
					if(tmpdiffnum!=1)
					{
						printf("This is fake neighboring relationship!\n");
						continue; //back to i<=readLen-k
					}*/
					diffpos=fixPos(localpos,revStr,k); //pos on the kmer
					if(diffpos>=k)
					{
						printf("Cannot find different pos between two k-mers in the correcting process!\n");
						//localpos++;
						continue; //back to for: i<=readLen-k
					}
					pmemcpy=memcpy(localpos,revStr,k);  //change the original refer!!
					if(pmemcpy==NULL)
					{
						printf("revStr to localpos memory copy failed in the correcting......\n");
						continue; //back to for: i<=readLen-k
					}
				 }
				if(diffnum==1) //neighbor
				{
					diffpos=fixPos(localpos,neighpos,k);
					if(diffpos>=k)
					{
						printf("Cannot find different pos between two k-mers in the correcting process!\n");
						//localpos++;
						continue; //back to for: i<=readLen-k
					}
					pmemcpy=memcpy(localpos,neighpos,k); //correct the kmer  change the original refer!!
					if(pmemcpy==NULL)
					{
						printf("neighstr to localpos memory copy failed in the correcting......\n");
						continue; //back to for: i<=readLen-k
					}
				}
				int tmppos=0;
				diffpos+=i; //pos on the read
				tmppos=diffpos-k+1;
				startpos=0>tmppos?0:tmppos;
				/*if(tmppos<0)
				{
					startpos=0;
				}
				else
				{
					startpos=tmppos;
				}*/
				if(diffpos<readLen-k)
				{
					endpos=diffpos;
				}
				else
				{
					endpos=readLen-k;
				}
				//Hashtable tmpKmer;
				//char *tmp; //absolute path
				lfreqnum=0;
				hfreqnum=0;
				for(j=startpos;j<=endpos;j++)  //check the neighbor kmers' freq effection
				{
					if(j==i)
					{
						continue;
					}
					tmp=prefer+j; //absolute path
				    tmpKmer=hash_table_lookup(tmp,k); //
					if(tmpKmer==NULL)
					{
						revSeq(tmp,revStr,k);
						tmpKmer=hash_table_lookup(revStr,k);
					}
					if(tmpKmer==NULL)
					{
						printf("Cannot find the effected adjacent node!It's a new kmer!\n");
						memcpy(localpos,orikmer,k); //restore
						lfreqnum=readLen-k+1;
						break; //break for j<i
				//hash_table_insert(ReferStr,tmp,revStr,k,tmp-ReferStr); //insert the new kmer,but its' neighbor haven't fixed yet!
						//continue; //back to for: j<i
					}
					if((tmpKmer->frequency)>=hfreqThre)
					{
						hfreqnum++;
				    }
					if((tmpKmer->frequency)<=lfreqThre)
					{
						lfreqnum++;
					}
				}//for: j [num,i]
				if(lfreqnum<=(int)((fixlevel/100.0)*(endpos-startpos))) //to correct the kmer
				{
					corrflag=1;
					//fix the kmer-net
					//memcpy(tmporikmer,localpos,k);
					//memcpy(localpos,orikmer,k);
					//printf("Remove the kmer:%s\n",orikmer);
					//hash_table_remove(ReferStr,orikmer,revStr,k);//remove the origine kmer from the kmer-network,but its neighbors cannot be fixed yet and the kmerNet isn't the newest!
					
					//memcpy(localpos,tmporikmer,k);
					//increase the frequency of neighbor
					
					//memcpy(localpos,orikmer,k);
			/*		for(j=startpos;j<=endpos;j++)
					{
						tmp=prefer+j;
						tmpKmer=hash_table_lookup(lookflag,tmp,k);
						if(tmpKmer==NULL)
						{
							revSeq(tmp,revStr,k);
							tmpKmer=hash_table_lookup(lookflag,revStr,k);
						}
						if(tmpKmer==NULL)
						{
							printf("Cannot find the changed kmer!\n");
							memcpy(localpos,orikmer,k);
							break;
						}
						tmpKmer->frequency++;
						tmpKmer->neighnum--;
					}
			*/
					/*tmpKmer=hash_table_lookup(ReferStr,localpos,k);
					if(tmpKmer==NULL)
					{
						revSeq(localpos,revStr,k);
						tmpKmer=hash_table_lookup(ReferStr,revStr,k);
					}
					if(tmpKmer==NULL)
					{
						printf("Cannot find the correted kmer!\n");
						pmemcpy=memcpy(localpos,orikmer,k); //cannot fix the kmerNet
						if(pmemcpy==NULL)
						{
							printf("orikmer to localpos memory copy failed in the correcting......\n");
							continue; //back to for: i<=readLen-k
						}
						//localpos++;
						continue;  //back to for: i<=readLen-k
					}*/

					/*Location lRNode;
					while(!(lRNode=calloc(1,sizeof(struct locateNode))))
					{
			            printf("malloc memory error for location node! try again!!\n");
					}
					memset(lRNode,0,sizeof(struct locateNode));
					lRNode->location=localpos-ReferStr; //absolute path
					if(tmpKmer->locations==NULL)
					{
						tmpKmer->locations=lRNode;
					}
					else
					{
						lRNode->lNext=tmpKmer->locations->lNext;
						tmpKmer->locations->lNext=lRNode;
					}*/
					/*tmpKmer->frequency+=freq;
					//tmpKmer->frequency++;
					tmpKmer->neighnum--; //but its neighbors haven't fiexed yet!*/
				}//if :lfreqnum<=...
				else  //not correct
				{
					pmemcpy=memcpy(localpos,orikmer,k);
					if(pmemcpy==NULL)
					{
						printf("localpos to orikmer memory copy failed in the correcting......\n");
						continue; //back to for: i<=readLen-k
					}
				}
			}//if: freq<lfreqThre
			//localpos++;
		}//for : i<=readLen-k
		if(corrflag==1)  //this read is corrected
		{
			ReadLine corread;
			while(!(corread=calloc(1,sizeof(struct readNode))))
			{
				printf("malloc memory error for corread node in the process of correction!Try again!\n");
			}
            corread->readId=pread->readId;
			corread->rNext=NULL;
	        if(precorrectedReads==NULL)
	        {
	           precorrectedReads=corread;
			}
	       else
	       	{
	           corread->rNext=precorrectedReads->rNext;
               precorrectedReads->rNext=corread;
	   		}
			//memcpy(porireads,oread,readLen); //store the origine of corrected reads
			
			/*char *tmpori=porireads;
			unsigned int m=0;
			printf("Correct!The original read:\n");
			for(m=0;m<readLen;m++)
			  printf("%c",*tmpori++);
			printf("\n");*/

			//porireads+=readLen;
			cornum++;
		}//if:corrflag
		else //if(corrflag==0)  //discard the reads
		{
			ReadLine disread;
			while(!(disread=calloc(1,sizeof(struct readNode))))
			{
				printf("malloc memory error for diread node in the process of correction!Try again!\n");
			}
			disread->readId=pread->readId;
			disread->rNext=NULL;
			if(prediscardReads==NULL)
			{
				prediscardReads=disread;
			}
			else
			{
				disread->rNext=prediscardReads->rNext;
				prediscardReads->rNext=disread;
			}
			//memcpy(prefer,oread,readLen); //restore original read	
			//printf("Not correct oriread:%s\n",oread);

		 }//else if : corrflag==0
		pread=pread->rNext;
	}//while : pread :corrigibleReads

	//free memory
	printf("The precorrected reads number:%lld\n",cornum);
	free(revStr);
	//free(tmporikmer);
	free(orikmer);
	return 0;
}//correct function


/**************************************************************/
/********************** The refiltering module ****************/
/**************************************************************/
int refilter(char *ReferStr,char *oriReferStr,unsigned int k,unsigned int refilterLevel,unsigned long int lfreqThre)
{
	unsigned long int i=0,lfreqnum=0;
	unsigned long int freq=0;
	unsigned long int id=0;
	unsigned long long int read=0;
	unsigned long long int correctednum=0;
	Hashtable kmer;
	ReadLine pRead;
	char *prefer,*localpos;
	char *revStr;
	while(!(revStr=(char *)malloc((k+1)*sizeof(char))))
	{
		printf("malloc memory error for revStr-array! Try again!\n");
	}
	memset(revStr,'\0',(k+1));
	//char *porireads;
	//porireads=orireads;
	pRead=precorrectedReads;
	while(pRead)
	{
		prefer=ReferStr;
    	id=pRead->readId;
    	read=id*readLen;
		
		//for(i=0;i<read;i++)
		  //prefer++;

		prefer+=read;
		lfreqnum=0;
		//localpos=prefer;
      	for(i=0;i<=readLen-k;i++) //i: relative path
      	{
			localpos=prefer+i;
			kmer=hash_table_lookup(localpos,k);
			if(kmer==NULL)
			{
				revSeq(localpos,revStr,k);//the reverse complement of the kmer
				kmer=hash_table_lookup(revStr,k);
     		}
			if(kmer==NULL)
			{
				printf("refiltering precorrectedReads....cannot find the kmer!\n");
				//lfreqnum=readLen-k+1;
				//localpos++;
				//continue; //back to for :i<=readLen-k
				//break; //go to :discard
				lfreqnum++;
				continue;
			}
     		freq=kmer->frequency;
			//printf("The kmer freq:%ld\n",freq);
			if(freq<=lfreqThre)
			{
				lfreqnum++;
			}
			//localpos++;
		 }//for :readLen-k
		if(lfreqnum>refilterLevel) //correction failed, discard reads
		{
			ReadLine disread;
			while(!(disread=calloc(1,sizeof(struct readNode))))
			{
				printf("malloc memory error for diread node in the process of refiltering!Try again!\n");
			}
			disread->readId=pRead->readId;
			disread->rNext=NULL;
			if(discardReads==NULL)
			{
				discardReads=disread;
			}
			else
			{
				disread->rNext=discardReads->rNext;
				discardReads->rNext=disread;
			}
			//memcpy(prefer,porireads,readLen); //cannot corrected
		 }//discard corrected reads
		else //accept corrected reads
		{
			ReadLine corread;
			while(!(corread=calloc(1,sizeof(struct readNode))))
			{
				printf("malloc memory error for corread node in the process of refiltering!Try again!\n");
			}
			corread->readId=pRead->readId;
			corread->rNext=NULL;
			if(correctedReads==NULL)
			{
				correctedReads=corread;
			}
			else
			{
				corread->rNext=correctedReads->rNext;
				correctedReads->rNext=corread;
			}
			correctednum++;
		 }//accept the corrected read
		 //porireads+=readLen;
	     pRead=pRead->rNext;
	}//while: pRead
   
	pRead=prediscardReads;
	//pRead=corrigibleReads;
	while(pRead)
	{
		prefer=oriReferStr;
    	id=pRead->readId;
    	read=id*readLen;
    	
		//for(i=0;i<read;i++)
		  //prefer++;
		
		prefer+=read;
		lfreqnum=0;
		//localpos=prefer;
      	for(i=0;i<=readLen-k;i++) //i: relative path
      	{
			localpos=prefer+i;
			kmer=hash_table_lookup(localpos,k);
			if(kmer==NULL)
			{
				revSeq(localpos,revStr,k);//the reverse complement of the kmer
				kmer=hash_table_lookup(revStr,k);
     		}
			if(kmer==NULL)
			{
				printf("refiltering prediscarded....cannot find the kmer!\n");
				//lfreqnum=readLen-k+1;
				lfreqnum++;
				continue; //back to for :i<=readLen-k
				//break;
			}
     		freq=kmer->frequency;
			if(freq<=lfreqThre)
			{
				lfreqnum++;
			}
			//localpos++;
		 }//for :readLen-k
		if(lfreqnum>refilterLevel) //discard reads
		{
			ReadLine disread;
			while(!(disread=calloc(1,sizeof(struct readNode))))
			{
				printf("malloc memory error for diread node in the process of refiltering!Try again!\n");
			}
			disread->readId=pRead->readId;
			disread->rNext=NULL;
			if(discardReads==NULL)
			{
				discardReads=disread;
			}
			else
			{
				disread->rNext=discardReads->rNext;
				discardReads->rNext=disread;
			}
		 }
		 else  //accept the reads
		  {
		  	 ReadLine accread;
			 while(!(accread=calloc(1,sizeof(struct readNode))))
			 {
			 	printf("malloc memory error for accept node in the porcess of refiltering!Try again!\n");
			 }
			 accread->readId=pRead->readId;
			 accread->rNext=NULL;
			 if(subacceptReads==NULL)
			 {
			 	subacceptReads=accread;
			 }
			 else
			 {
			 	accread->rNext=subacceptReads->rNext;
				subacceptReads->rNext=accread;
			 }
		  }//else
		  pRead=pRead->rNext;
	   }//while: pRead
	printf("The final corrected reads number:%lld\n",correctednum);
	free(revStr);
	return 0;
}

/*************************************************************/
/******************** The output module **********************/
/*************************************************************/
void printKmerNet(char * oriReferStr,unsigned int k,FILE *kmerNetOut)
{
	unsigned long long int i,j;
	char *tmpos;
	Hashtable pHead;
	unsigned int m;
	//printf("Output the origine kmer Net......\n");
	fprintf(kmerNetOut,">The read-length:%d\n",readLen);  
	fprintf(kmerNetOut,">The reads-num:%lld\n",readsnumber);
	fprintf(kmerNetOut,">The k-value:%d\n",k);
	fprintf(kmerNetOut,">Theoretical Max theoryKmer-num:%lld\n",(readLen-k+1)*readsnumber);
	fprintf(kmerNetOut,">The realKmer-num:%lld\n",hash_table_size);
	for(i=0;i<HASH_TABLE_MAX_SIZE;i++)
	{
		pHead=hashTable[i];
		while(pHead)
		{
			fprintf(kmerNetOut,">Id=%ld;frequency=%ld\nkmer:",pHead->pKey-oriReferStr,pHead->frequency);
			tmpos=pHead->pKey;
			for(j=0;j<k;j++)
			  fprintf(kmerNetOut,"%c",*tmpos++);
			fprintf(kmerNetOut,"\nLocations:");
			Location tmplocate;
			tmplocate=pHead->locations;
			j=0;
			while(tmplocate)
			{
				fprintf(kmerNetOut,"%lld;",tmplocate->location);
				j++;
				if(j%80==0)
				 fprintf(kmerNetOut,"\n");
				tmplocate=tmplocate->lNext;
			}
			fprintf(kmerNetOut,"\nNeighbors:%d\n",pHead->neighnum);
			Neighbor tmpNeigh;
			//Hashtable tmp;
			tmpNeigh=pHead->neighbors;
			while(tmpNeigh)
			{ 
			   	tmpos=tmpNeigh->neighbor;
			   	/*tmp=hash_table_lookup(ReferStr,tmpos,k);
				if(tmp==NULL)
				{
					revSeq(tmpos,revStr,k);
					tmp=hash_table_lookup(ReferStr,revStr,k);
				}
				if(tmp==NULL)
				{
					tmpNeigh=tmpNeigh->nNext;
					continue;
				}*/
				fprintf(kmerNetOut,"#Id=%ld;frequency=%ld\n",tmpNeigh->neighbor-oriReferStr,tmpNeigh->frequency);
			   	for(m=0;m<k;m++)
			     	fprintf(kmerNetOut,"%c",*tmpos++);
			   	fprintf(kmerNetOut,"\n");
				tmpNeigh=tmpNeigh->nNext;
			}//while :tmpNeigh
			fprintf(kmerNetOut,"\n");
			pHead=pHead->pNext;
		}//while (pHead)
	}//for:i<HASH_MAX_SIZE
	printf("Finish outputing kmerNet!\n");
}
void printAcceptReads(char * oriReferStr,unsigned int k,FILE *accReadsOut)
{
	unsigned long long int scale;
	unsigned long int i;
	unsigned long long int accnum=0;
	char *prefer;
	printf("Output the acceptReads......\n");
	ReadLine accpread;
	accpread=acceptReads;
	while(accpread)
	{
		scale=(accpread->readId)*readLen;
		prefer=oriReferStr+scale;
		fprintf(accReadsOut,">%ld\n",accpread->readId);
		for(i=0;i<readLen;i++)
		{
			fprintf(accReadsOut,"%c",*prefer++);
		}
		fprintf(accReadsOut,"\n");
		accpread=accpread->rNext;
		accnum++;
	}
	printf("The accepted reads number:%lld\n",accnum);
	printf("Finish outputing acceptReads!\n");
}
void printsubAcceptReads(char * oriReferStr,unsigned int k,FILE *subaccReadsOut)
{
	unsigned long long int scale;
	unsigned long int i=0;
	unsigned long long int accnum=0;
	char *prefer;
	printf("Output the subacceptReads......\n");
	ReadLine accpread;
	accpread=subacceptReads;
	while(accpread)
	{
		scale=(accpread->readId)*readLen;
		prefer=oriReferStr+scale;
		fprintf(subaccReadsOut,">%ld\n",accpread->readId);
		for(i=0;i<readLen;i++)
		{
			fprintf(subaccReadsOut,"%c",*prefer++);
		}
		fprintf(subaccReadsOut,"\n");
		accpread=accpread->rNext;
		accnum++;
	}
	printf("The sub accepted reads number:%lld\n",accnum);
	printf("Finish outputing subacceptReads!\n");
}
void printDiscardReads(char * oriReferStr,unsigned int k,FILE *disReadsOut)
{
	unsigned long long int scale;
	unsigned long int i;
	unsigned long long int disnum=0;
	char *prefer;
	printf("Output the discardReads......\n");
	ReadLine disread;
	disread=discardReads;
	while(disread)
	{
		scale=(disread->readId)*readLen;
		prefer=oriReferStr+scale;
		fprintf(disReadsOut,">%ld\n",disread->readId);
		for(i=0;i<readLen;i++)
		{
			fprintf(disReadsOut,"%c",*prefer++);
		}
		fprintf(disReadsOut,"\n");
		disread=disread->rNext;
		disnum++;
	}
	printf("The dicarded reads number:%lld\n",disnum);
	printf("Finish outputing the discardReads!\n");
}
void printprecorReads(char * ReferStr,unsigned int k,FILE *precorReadsOut)
{
	unsigned long long int scale;
	unsigned long int i;
	char *prefer;
	unsigned long long int precornum=0;
	ReadLine precorread;
	printf("Output the precorReads......\n");
	precorread=precorrectedReads;
	while(precorread)
	{
		scale=(precorread->readId)*readLen;
		prefer=ReferStr+scale;
		fprintf(precorReadsOut,">%ld\n",precorread->readId);
		for(i=0;i<readLen;i++)
		{
			fprintf(precorReadsOut,"%c",*prefer++);
		}
		fprintf(precorReadsOut,"\n");
		precorread=precorread->rNext;
		precornum++;
	}
	printf("The precorrected reads number:%lld\n",precornum);
	printf("Finish outputing the precorReads!\n");
}
void printCorrectReads(char * ReferStr,unsigned int k,FILE *corReadsOut)
{
	unsigned long long int scale;
	unsigned long int i;
	unsigned long long int cornum=0;
	char *prefer;
	ReadLine corread;
	printf("Output the correctedReads......\n");
	corread=correctedReads;
	while(corread)
	{
		scale=(corread->readId)*readLen;
		prefer=ReferStr+scale;
		fprintf(corReadsOut,">%ld\n",corread->readId);
		for(i=0;i<readLen;i++)
		{
			fprintf(corReadsOut,"%c",*prefer++);
		}
		fprintf(corReadsOut,"\n");
		corread=corread->rNext;
		cornum++;
	}
	printf("The corrected reads number:%lld\n",cornum);
	printf("Finish outputing the correctedReads!\n");
}
/*****************************************************/
/**************** Initiation Module********************/
/*****************************************************/
unsigned long long int Initiation(char * ReferStr,unsigned int k,char *readsLibName)
{
    
	unsigned long long int  i=0,n=0;
	//unsigned long long int j=0;
	//char *ReferStr=(char *)malloc(MAX_REFER_LEN*sizeof(char));
	char *tmpStr;
	while(!(tmpStr=(char *)malloc(MAX_LEN*sizeof(char))))
	{
		printf("malloc memory error for tmpStr-array!Try again!\n");
	}
    //memset(tmpStr,'\0',MAX_LEN);
	//unsigned long long int *seqNo=(unsigned long long int*)malloc(MAX_SEQ_NUM*sizeof(unsigned long long int));
	//unsigned long long int readsnumber=0;
	//unsigned int readLen=0;
	FILE *readsIn;
    char *pstr;
	char *pmemcpy;
	//char *opmemcpy;

	hash_table_size=0;
	//open the file
	if((readsIn=fopen(readsLibName,"r"))==NULL)
	{
		printf("Cannot open the reads library file!\n");
		exit(1);
	}
	i=0;
	n=0;
	pstr=ReferStr;
	//opstr=oriReferStr;
	readsnumber=0;
	readLen=0;
	//char *tmp;
	while(!feof(readsIn))
	{
        //fscanf(referIn,"%s",tmpStr);
		fgets(tmpStr,NAME_LEN,readsIn);
		tmpStr[strlen(tmpStr)]='\0';
		//printf("seq:%s",tmpStr);
		if(tmpStr[0]=='>')
         {
		 	  readLen=0;
              //memset(tmpStr,'\0',strlen(tmpStr));
			//printf("The sequence:%s\n",tmpStr);
            fscanf(readsIn,"%s",tmpStr);
			//fgets(tmpStr,MAX_LEN,readsIn);	
			//printf("The real seq:%s",tmpStr);
			while(tmpStr[0]!='#')
				{
                	i=strlen(tmpStr);
					//printf("The line length:%ld\n",i);
					tmpStr[i]='\0';
					//printf("The real seq:%s\n",tmpStr);
					pmemcpy=memcpy(pstr,tmpStr,i);
					if(pmemcpy==NULL)
					{
						printf("tmpStr to pstr memory copy failed in the initiating......\n");
						continue; 
					}

					/*opmemcpy=memcpy(opstr,tmpStr,i);
					if(opmemcpy==NULL)
					{
						printf("tmpStr to opstr memory copy failed in the initiating.....\n");
						continue;
					}*/

                 	//memset(tmpStr,'\0',MAX_LEN);
					/*tmp=pstr;
					printf("The refer segment:\n");
					for(j=0;j<i;j++)
						printf("%c",*tmp++);*/
					pstr+=i;
					//opstr+=i;
					n+=i;
                 	readLen+=i;
					fscanf(readsIn,"%s",tmpStr);
					//fgets(tmpStr,MAX_LEN,readsIn);
            	 }//while:tmpStr
				//seqNo[seqnumber++]=seqLen;
		        //printf("The seqlen of sequence:%d\n",readLen);
				readsnumber++;
		   }//if:tmpStr
		   //printf("The seqlen of sequence:%lld\n",seqLen);
	  }//while

		ReferStr[n]='\0';
		//oriReferStr[n]='\0';
	 for(i=0;ReferStr[i]!='\0';i++)
	   ReferStr[i]=toupper(ReferStr[i]);
	 //for(i=0;oriReferStr[i]!='\0';i++)
	   //oriReferStr[i]=toupper(oriReferStr[i]);
     hash_table_init();
 	
	//printf("The reference length is:%lld\n",n);
 	//printf("The reference is:\n%s\n",ReferStr);
	
	 free(tmpStr);
	 fclose(readsIn);
	 return n;
}

/*******************************************************/
/****************** Scheduling Module ******************/
/*******************************************************/
int schedule(unsigned int k,char *readsLibName,unsigned int lfreq,unsigned int hfreq,unsigned int accThre,unsigned int disThre,unsigned int fixscale,unsigned int refilterscale)
{
	/************defination and allocate memory*******/
	printf("Scheduling.......\n");
		//define the output files' names
	char prekmerNetName[256],postkmerNetName[256],accReadsName[256],disReadsName[256],distmpReadsName[256],corReadsName[256],precorReadsName[256],subacceptReadsName[256];
	char prekmersuffix[]=".prekmerNet";
	char postkmersuffix[]=".postkmerNet";
	char accsuffix[]=".acceptReads";
	char dissuffix[]=".discardReads";
	char distmpsuffix[]=".discardReadstmp";
	char corsuffix[]=".correctReads";
	char precorsuffix[]=".precorReads";
	char subaccsuffix[]=".subaccReads";
	strcpy(prekmerNetName,readsLibName);
	strcpy(postkmerNetName,readsLibName);
	strcpy(accReadsName,readsLibName);
	strcpy(disReadsName,readsLibName);
	strcpy(distmpReadsName,readsLibName);
	strcpy(corReadsName,readsLibName);
	strcpy(precorReadsName,readsLibName);
	strcpy(subacceptReadsName,readsLibName);

	strcat(prekmerNetName,prekmersuffix);
	strcat(postkmerNetName,postkmersuffix);
	strcat(accReadsName,accsuffix);
	strcat(disReadsName,dissuffix);
	strcat(distmpReadsName,distmpsuffix);
	strcat(corReadsName,corsuffix);
	strcat(precorReadsName,precorsuffix);
	strcat(subacceptReadsName,subaccsuffix);

		//definition of scales
	long int tmplfreq;
	unsigned long int hfreqThre,lfreqThre;
	unsigned int fixlevel,acceptThre,discardThre;
	unsigned int refilterLevel;

		//allocate memory
	char *ReferStr;
	while(!(ReferStr=(char*)malloc(MAX_REFER_LEN*sizeof(char))))
	{
		printf("malloc memory error for ReferStr-array!Try again!\n");
	}
	
	char *revStr;
	while(!(revStr=(char *)malloc((k+1)*sizeof(char))))
	{
		printf("malloc memory error for revStr-array!Try again!\n");
	}
	//unsigned int readLen;
	//unsigned long long int readsnumber;
	unsigned long int *kmerFreq;
	unsigned long int *fiveNum;
	while(!(fiveNum=(unsigned long int *)malloc(8*sizeof(unsigned long int))))
	{
		printf("malloc memory error for fiveNum-array!Try again!\n");
	}
	char *oriReferStr;
	while(!(oriReferStr=(char *)malloc(MAX_REFER_LEN*sizeof(char))))
	{
		printf("malloc memory error for orireferStr-array!Try again!\n");
	}
	unsigned long long int referLen=0;

		//define file pointers and open files
	FILE *prekmerNetOut,*postkmerNetOut,*accReadsOut,*disReadsOut,*distmpReadsOut,*corReadsOut,*precorReadsOut,*subaccReadsOut;
	//printf("Open files!\n");
	if(!(prekmerNetOut=fopen(prekmerNetName,"w+")))
	{
		printf("Cannot open the prekmerNet file !!\n");
		exit(1);
	}
	/*if(!(postkmerNetOut=fopen(postkmerNetName,"w+")))
	{
		printf("Cannot open the postkmerNet file !!\n");
		exit(1);
	}*/
	if(!(accReadsOut=fopen(accReadsName,"w+")))
	{
		printf("Cannot open the acceptReads file !!\n");
		exit(1);
	}
	if(!(disReadsOut=fopen(disReadsName,"w+")))
	{
		printf("Cannot open the discardReads file !!\n");
		exit(1);
	}
/*	if(!(distmpReadsOut=fopen(distmpReadsName,"w+")))
	{
		printf("Cannot open the discardReads file !!\n");
		exit(1);
	}
*/
	if(!(corReadsOut=fopen(corReadsName,"w+")))
	{
		printf("Cannot open the correctReads file !!\n");
		exit(1);
	}
/*	if(!(precorReadsOut=fopen(precorReadsName,"w+")))
	{
		printf("Cannot open the precorrectReads file !!\n");
		exit(1);
	}
*/
	if(!(subaccReadsOut=fopen(subacceptReadsName,"w+")))
	{
		printf("Cannot open the subacceptReads file !!\n");
		exit(1);
	}

	/********** Initiation *************/
	//char *pmemcpy;
	printf("Starting the Initiation Module......\n");
	referLen=Initiation(ReferStr,k,readsLibName);
    strcpy(oriReferStr,ReferStr);	
	/*pmemcpy=memcpy(oriReferStr,ReferStr,(referLen+1));
	if(pmemcpy==NULL)
	{
		printf("ReferStr to oriReferStr memory copy failed in the scheduling.....\n");
		return 1;
	}*/
	printf("Finishing Initiation!\n");

	/******** KmerNet creation *********/
	printf("Starting the KmerNet creation Module......\n ");
	createKmerNet(oriReferStr,k,revStr);
	printf("Output the original kmer network......\n");
	printKmerNet(oriReferStr,k,prekmerNetOut); //output the origine kmerNet
	//backup hashTable
/*	unsigned long long int i=0;
	char *oriStr;
	oriStr=oriReferStr;
	Hashtable pHead;
	Hashtable kmerNode;
	Neighbor newneigh,oldneigh;
	unsigned int neighnum=0; 
	for(i=0;i<HASH_TABLE_MAX_SIZE;i++)
	{
		pHead=hashTable[i];
		while(pHead)
		{
			while(!(kmerNode=calloc(1,sizeof(struct HashNode))))
			{
				printf("Memory error for hashNode in scheduing.....Try again!\n");
			}
			kmerNode->pKey=oriStr+(pHead->pKey-ReferStr);
			kmerNode->frequency=pHead->frequency;
			neighnum=pHead->neighnum;
			kmerNode->neighnum=neighnum;
			oldneigh=pHead->neighbors;
			while(neighnum)
			{
				while(!(newneigh=calloc(1,sizeof(struct neighborNode))))
				{
					printf("Memory error for neigh in scheduind.....Try again!\n");
				}
				newneigh->neighbor=oriStr+(oldneigh->neighbor-ReferStr);
				newneigh->frequency=oldneigh->frequency;
				if(kmerNode->neighbors==NULL)
				{
					kmerNode->neighbors=newneigh;
				}
				else
				{
					newneigh->nNext=kmerNode->neighbors->nNext;
					kmerNode->neighbors->nNext=newneigh;
				}
				neighnum--;
				oldneigh=oldneigh->nNext;
			}
			kmerNode->pNext=NULL;
			if(orihashTable[i]==NULL)
			    orihashTable[i]=kmerNode;
			else
			{
				kmerNode->pNext=orihashTable[i]->pNext;
				orihashTable[i]->pNext=kmerNode;
			}
			pHead=pHead->pNext;
		}
     }
*/

	//test
/*	unsigned long long int t=0;
	  for(i=0;i<HASH_TABLE_MAX_SIZE;i++)
	  {
	  	pHead=orihashTable[i];
		while(pHead)
		{
			t++;
			pHead=pHead->pNext;
		}
	  }
	 printf("The origine kmer number:%lld\n",t);
*/	//fclose(prekmerNetOut);
	printf("Finishing KmerNet creation!\n");

	/******* Correction Process ********/
	printf("Starting Correction process......\n");
	 	//filtering
	printf("Filtering......\n");
	while(!(kmerFreq=(unsigned long int *)malloc(hash_table_size*sizeof(unsigned long int))))
	{
		printf("malloc memory error for kmerFreq-array!Try again!\n");
	}
	countFivenum(kmerFreq,fiveNum);
	hfreqThre=(int)(fiveNum[1]/hfreq*1.0); //Q1/hfreq
	tmplfreq=0>(fiveNum[1]-(int)(1.5*fiveNum[5]))?0:(fiveNum[1]-(int)(1.5*fiveNum[5]));
	lfreqThre=lfreq>tmplfreq?lfreq:tmplfreq;
	//lfreqThre=lfreq;
	acceptThre=(int)((accThre/100.0)*(readLen-k+1));
	discardThre=(int)((disThre/100.0)*(readLen-k+1));
	refilterLevel=(int)((refilterscale/100.0)*(readLen-k+1));
	fixlevel=fixscale;//when the fixlevel persentage neighbors kmer freqs on the reads <lfreqThre,then correct the kmer

	printf("k-value:%d\n",k);
	printf("lfreThre:%ld\n",lfreqThre);
	printf("hfreqThre:%ld\n",hfreqThre);
	printf("acceptThre:%d\n",acceptThre);
	printf("fixlvel-weight:%d\n",fixlevel);
	printf("refilterNum:%d\n",refilterLevel);

	prefilter(oriReferStr,k,hfreqThre,lfreqThre,acceptThre,discardThre);
	printAcceptReads(oriReferStr,k,accReadsOut);
	//fclose(accReadsOut);
//	printDiscardReads(oriReferStr,k,distmpReadsOut);
	//fclose(distmpReadsOut);
	printf("Finishing filtering module!\n");
		//fixing ,change ReferStr
	//lfreqThre=1;
	//hfreqThre=fiveNum[1]; //Q1
	//fixlevel=fixscale;//when the fixlevel persentage neighbors kmer freqs on the reads <lfreqThre,then correct the kmer
	printf("Correcting......\n");
	correct(ReferStr,k,fixlevel,lfreqThre,hfreqThre);
//	printprecorReads(ReferStr,k,precorReadsOut);
	//fclose(precorReadsOut);
	printf("Finishing correction module!\n");
	   //refilteing	
	//refilterLevel=(int)((refilterscale/100.0)*(readLen-k+1));
	
	//if(!(strcmp(ReferStr,oriReferStr)))
	// printf("The oriReferStr is effected!!\n");


	printf("Refiltering......\n");
	refilter(ReferStr,oriReferStr,k,refilterLevel,lfreqThre);
	printCorrectReads(ReferStr,k,corReadsOut);
	printDiscardReads(oriReferStr,k,disReadsOut);
	printsubAcceptReads(oriReferStr,k,subaccReadsOut);
	//fclose(corReadsOut);
	//fclose(disReadsOut);
	printf("Finishing refiltering module!\n");


	/************ Output kmerNet **********/
	//printf("Output the final kmerNet......\n");
	//printKmerNet(ReferStr,k,revStr,postkmerNetOut);
	//fclose(postkmerNetOut);
/*	unsigned long long int t=0;
	  for(i=0;i<HASH_TABLE_MAX_SIZE;i++)
	  {
	  	pHead=orihashTable[i];
		while(pHead)
		{
			t++;
			pHead=pHead->pNext;
		}
	  }
	printf("The final original kmer-network kmers number:%lld\n",t);
*/
//free memory and close files
	free(ReferStr);
	free(revStr);
	free(oriReferStr);
	free(kmerFreq);
	hash_table_release();
	fclose(prekmerNetOut);
	//fclose(postkmerNetOut);
	fclose(accReadsOut);
	//fclose(distmpReadsOut);
	fclose(disReadsOut);
	//fclose(precorReadsOut);
	fclose(corReadsOut);
	return 0;
}
int main(int argc,char **argv)
{
    
	 int c;   // c: the option  
	unsigned int k;  //k: define the length of k-mer c: the option  
	unsigned int lfreq; //user-specified lowfreq [1,2....]
	unsigned int hfreq; //hfreq weight : Q1/hfreq [1,2....]
	unsigned int accThre; //acceptThre weight: (accThre/10)*(readLen-k+1) [1...10]
	unsigned int disThre; //discardThre weight: (disThre/10)*(readLen-k+1)[1...10]
	unsigned int fixlevel; //fix weight: (fixlevel/10)*(endpos-startpos) [1...10]
	unsigned int refilter; //refilter weight:(refilter/10)*(readLen-k+1) [1...10]
	genome_net_opt  opt;
	char readsName[256];
 
    
	strcpy(readsName,argv[15]);
	
	printf("The input file:%s\n",readsName);


	if((c=getopt(argc,argv,"k:l:h:a:d:f:r:"))<0)
	{
		printf("Usage:\n	./readsCorrect -k 3 -l 1 -h 2 -a 10 -d 10 -f 1 -r 2  infile\n");
		printf(" 	            	-k define the length of k-mer c: the option\n");
		printf("			-l user-specified lowfreq [1,2....]\n");
		printf("			-h hfreq weight : Q1/hfreq [1,2....]\n");
		printf("			-a acceptThre weight: (accThre/100)*(readLen-k+1) [1...100]\n");
		printf("			-d discardThre weight: (disThre/100)*(readLen-k+1)[1...100]\n");
		printf("			-f fix weight: (fixlevel/100)*(endpos-startpos) [1...100]\n");
		printf("	    		-r refilter weight:(refilter/100)*(readLen-k+1) [1...100]\n");
		return 0;
	}
    while(c>=0)
	{
		switch(c)
		{
			case 'k':opt.k=atoi(optarg); 
			 			break;
			case 'l':opt.l=atoi(optarg);
						break;
			case 'h':opt.h=atoi(optarg);
						break;
			case 'a':opt.a=atoi(optarg);
						break;
			case 'd':opt.d=atoi(optarg);
						break;
			case 'f':opt.f=atoi(optarg);
						break;
			case 'r':opt.r=atoi(optarg);
						break;
		};
		c=getopt(argc,argv,"k:l:h:a:d:f:r:");
	}
    k=opt.k;
	lfreq=opt.l;
	hfreq=opt.h;
	accThre=opt.a;
	disThre=opt.d;
	fixlevel=opt.f;
	refilter=opt.r;
	printf("Starting......\n");
    schedule(k,readsName,lfreq,hfreq,accThre,disThre,fixlevel,refilter); //schedule the core modules
	printf("finished!!\n");
	return 0;
}
