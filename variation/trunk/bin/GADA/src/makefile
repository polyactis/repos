
BoostPythonTarget	= GADA.so
BoostPythonLib  = -lboost_python

CFLAGS	= -O3
CXXFLAGS	= -fPIC $(CFLAGS)
SharedLibFlags	= -shared -fPIC
IncludeDirs = `python-config --includes`
CXXCOMPILER = g++

all:	GADA $(BoostPythonTarget)

#for boost.python target(shared libraries)
$(BoostPythonTarget):	%.so:	%.o
	-$(CXXCOMPILER) $(BoostPythonLib) $(SharedLibFlags) -lm $< BaseGADA.o -o $@


.cc.o:
	-$(CXXCOMPILER) $(CXXFLAGS) $(IncludeDirs) $(SharedLibFlags) -c $< -o $@

.c.o:
	-$(CXXCOMPILER) $(CXXFLAGS) $(IncludeDirs) $(SharedLibFlags) -c $< -o $@

GADA:	GADA.o BaseGADA.o BaseGADA.h consoledefines.h
	-$(CXXCOMPILER) $(CXXFLAGS) $(IncludeDirs) $(BoostPythonLib) `python-config --libs` -o GADA GADA.o BaseGADA.o -lm

clean:
	-rm GADA *.o $(BoostPythonTarget)