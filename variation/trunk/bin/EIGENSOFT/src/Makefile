HOMEL=$(PWD)
DEBUG_OPTIONS= -g -pg
BIN=$(HOMEL)/bin
# "make eigenstrat" to make eigenstrat program (ditto for other programs)
# "mv eigenstrat ../bin" (or "make install" to copy all C programs to ../bin)
# "make clean" to clean up extra files in this directory
# "make clobber" to clobber all files and subdirectories except source files
#    so as to enable recompiling from scratch.
.PHONY: nicklib dirs rmdirs clobber

PROGS= smartpca convertf twstats eigenstrat eigenstratQTL pca
NLIB=$(HOMEL)/smartlib/nicklib.a
all: $(NLIB) $(PROGS)

IDIR=$(HOMEL)/smartinclude
VPATH=.:nicksrc
BLAS = blas-3
# may need to change to BLAS = blas (depends on blas/lapack installation)

CFLAGS= -c -O -I$(IDIR) -Wimplicit
OBJ=strsubs.o sortit.o vsubs.o statsubs.o linsubs.o getpars.o xsearch.o gauss.o	gds.o
TWTAB=\"$(HOMEL)/smarttables/twtable\"

statsubs.o:     nicksrc/statsubs.c
	$(CC)  $(CFLAGS) -DTWTAB=$(TWTAB) -o statsubs.o nicksrc/statsubs.c

M1=smartpca
M1O=smartpca.o  twsubs.o mcio.o admutils.o egsubs.o eigsubs.o  eigx.o 

M2=convertf
M2O=convertf.o  mcio.o admutils.o

M3=twstats
M3O=twstats.o  

M4=eigenstrat
M4O=eigenstrat.o

M5=eigenstratQTL
M5O=eigenstratQTL.o

M6=pca
M6O=pca.o eigsubs.o eigx.o

M7=pcatoy
M7O=pcatoy.o eigsubs.o eigx.o

CC=gcc
       

install:	all
	cp $(PROGS) $(HOMEL)/bin

uninstall:   
	rm -f $(BIN)/smartpca
	rm -f $(BIN)/convertf
	rm -f $(BIN)/twstats
	rm -f $(BIN)/eigenstrat
	rm -f $(BIN)/eigenstratQTL
	rm -f $(BIN)/pca
	rm -f $(NLIB)

$(M1): $(NLIB) $(M1O)
	$(CC) -O -I$(IDIR) $(DEBUG_OPTIONS) -o $(M1) $(M1O) $(NLIB) /usr/lib/libf2c.a  -lm -llapack -l$(BLAS) -Wimplicit # -lf2c
	#2008-12-02 -lf2c doesn't work. ERROR: "/usr/bin/../lib/libf2c.so: undefined reference to `MAIN__'"

$(M2): $(NLIB) $(M2O)
	$(CC) -O -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M2) $(M2O) $(NLIB) -Wimplicit

$(M3): $(NLIB) $(M3O)
	$(CC) -O -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M3) $(M3O) $(NLIB) -Wimplicit

$(M4): $(NLIB) $(M4O)
	$(CC) -O -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M4) $(M4O) $(NLIB) -Wimplicit

$(M5): $(NLIB) $(M5O)
	$(CC) -O -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M5) $(M5O) $(NLIB) -Wimplicit

$(M6): $(NLIB) $(M6O)
	$(CC) -O -I$(IDIR) $(DEBUG_OPTIONS) -o $(M6) $(M6O) $(NLIB) /usr/lib/libf2c.a -lm -llapack -l$(BLAS) -Wimplicit #-lf2c 

$(M7): $(NLIB) $(M7O)
	$(CC) -O -I$(IDIR) $(DEBUG_OPTIONS) -o $(M7) $(M7O) $(NLIB) -lm -llapack -l$(BLAS) -lf2c -Wimplicit

$(NLIB):	$(OBJ)	#2008-12-02 phony target dirs triggers compilation no matter whether any source file has been modified. removed.
	mkdir -p  $(HOMEL)/smartlib
	mkdir -p  $(HOMEL)/smarttables
	mkdir -p  $(HOMEL)/smartinclude
	mkdir -p  $(BIN)
	cp twtable  $(HOMEL)/smarttables
	cp nicksrc/*.h  $(IDIR)	ar -r libnick.a $(OBJ)
	cp libnick.a  $(NLIB)

#nicklib:	dirs   libnick.a 	#2008-12-02 phony target nicklib is totally eliminated by using $(NLIB).

#dirs:	#2008-12-02 phony target dirs triggers compilation no matter whether any source file has been modified. removed.


clean: 
	rm -f *.o 
	rm -f core
	rm -f libnick.a
	rm -f $(PROGS)

clobber: clean rmdirs uninstall


rmdirs: 
	rm -rf $(HOMEL)/smartlib 
	rm -rf $(HOMEL)/smarttables 
	rm -rf $(HOMEL)/smartinclude 

