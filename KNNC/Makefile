CC=gcc
CFILE1=hashTable
CFILE2=KNNC
$(CFILE2):		$(CFILE1).o		$(CFILE2).o
	$(CC)	-Wall -mcmodel=medium  -g  $(CFILE1).o	$(CFILE2).o    -o $(CFILE2)
$(CFILE1).o:	$(CFILE1).c		$(CFILE1).h
	$(CC)	-Wall -mcmodel=medium  -o  $(CFILE1).o  -c $(CFILE1).c
$(CFILE2).o:	$(CFILE2).c		$(CFILE2).h
	$(CC) 	-Wall -mcmodel=medium  -o  $(CFILE2).o	-c  $(CFILE2).c
clean:
	rm	-f *.o
cleanall:
	rm -f *.o
	rm -f $(CFILE2)
