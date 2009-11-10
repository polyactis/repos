
BoostPythonTarget	= GADA.so
BoostPythonLib  = -lboost_python

CFLAGS	= -O3
CXXFLAGS	= -O3 -fPIC
SharedLibFlags	= -shared -fPIC
IncludeDirs = `python-config --includes`
CXXCOMPILER = g++

all:	GADA $(BoostPythonTarget)

#for boost.python target(shared libraries)
$(BoostPythonTarget):	%.so:	%.o
	-$(CXXCOMPILER) $(BoostPythonLib) $(SharedLibFlags) -lm $< BaseGADA.o -o $@


.cc.o:
	-$(CXXCOMPILER) $(CFLAGS) $(IncludeDirs) $(SharedLibFlags) -c $< -o $@

.c.o:
	-$(CXXCOMPILER) $(CFLAGS) $(IncludeDirs) $(SharedLibFlags) -c $< -o $@

GADA:	GADA.o BaseGADA.o BaseGADA.h consoledefines.h
	-$(CXXCOMPILER) $(CFLAGS) $(IncludeDirs) $(BoostPythonLib) $(SharedLibFlags) -o GADA GADA.o BaseGADA.o -lm

clean:
	-rm GADA *.o $(BoostPythonTarget)