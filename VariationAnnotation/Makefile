#
# Makefile for program design
#
# Turn on debugging

# icc -w -vec-report0 -DMKL_VML_MIXED -I/opt/intel/mkl/10.0.4.023/include GeneComparison.cpp -o map.exe -Wl,--start-group /opt/intel/mkl/10.0.4.023/lib/em64t/libmkl_intel_lp64.a /opt/intel/mkl/10.0.4.023/lib/em64t/libmkl_intel_thread.a /opt/intel/mkl/10.0.4.023/lib/em64t/libmkl_core.a -Wl,--end-group /opt/intel/mkl/10.0.4.023/lib/em64t/libguide.a -lpthread -lm

# CC=g++34
# CFLAGS=-g -O3 -I/usr/include/mysql -Wno-deprecated

# CC=LC_ALL=C icpc

CC=gcc
CFLAGS=-w  -lpthread -lm -I/pub/include/mysql  -std=gnu99 

LIBAGS=-L/usr/lib64 -L/pub/lib64/mysql -lmysqlclient

objects = main.o variant.o  mysql.o

GA : $(objects)
	$(CC) $(CFLAGS) $(LIBAGS) -g -o ./annVar  $(objects)

###############################################
# modules
main.o : variant.h main.c
	$(CC) $(CFLAGS) -g  -c main.c

variant.o : variant.h variant.c
	$(CC) $(CFLAGS) -g -c variant.c

mysql.o : variant.h mysql.c
	$(CC) $(CFLAGS) -g -c mysql.c

###############################################
