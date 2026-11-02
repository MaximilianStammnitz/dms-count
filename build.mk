# dms_count - C/htslib implementation of deepmutscan's count_variants.py (byte-identical output).
# Needs htslib headers and library (1.17+ for the thread-pool API; tested with 1.21):
#   make -f build.mk HTSLIB=/path/to/htslib   # a built htslib source tree (libhts.a inside)
#   make -f build.mk HTSLIB_PREFIX=/usr/local # an installed htslib (include/htslib, lib/libhts.so)
CC      ?= gcc
CFLAGS  ?= -O2 -Wall -Wno-misleading-indentation
ifdef HTSLIB
  INC = -I$(HTSLIB)
  LIB = $(HTSLIB)/libhts.a -lz -lbz2 -llzma -lpthread -lm
else
  HTSLIB_PREFIX ?= /usr/local
  INC = -I$(HTSLIB_PREFIX)/include
  LIB = -L$(HTSLIB_PREFIX)/lib -lhts -lpthread -lm
endif
dms_count: dms_count.c
	$(CC) $(CFLAGS) -o $@ $< $(INC) $(LIB)
clean:
	rm -f dms_count
