all:	PCApolya

.cpp.o:
	-g++ $(CFLAGS) -c $< -o $@

PCApolya:	%:	%.o
	-g++ $(CFLAGS)  -o $@ $<

clean:
	-rm PCApolya *.o
