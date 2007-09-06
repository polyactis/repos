c++tg=

CFLAGS=-O3

all:	$(c++tg) forwsim


.c.o:
	-gcc $(CFLAGS) -c $< -o $@

.cpp.o:
	-g++ $(CFLAGS) -c $< -o $@

forwsim:	forwsim.o mtrand.o
	-g++ $(CFLAGS) -o $@ $^

$(c++tg): 	%:	%.o
	-g++ $(CFLAGS) $(IDIR) $< -o $@ 

$(ctg):	%:	%.o
	-gcc $(CFLAGS)  $< -o $@




clean:
	-rm -f *~  *.bb *.bbg *.da *.gcov *.o forwsim $(ctg) $(c++tg)
