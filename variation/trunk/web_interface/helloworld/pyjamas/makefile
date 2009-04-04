dir=../helloworld/public/Accession/

all: $(dir)Accession.html

.py.html:
	-pyjsbuild $< -o $(dir) $(FLAGS)

$(dir)Accession.html:	Accession.py
	-pyjsbuild $< -o $(dir) $(FLAGS)

clean:
	-rm $(dir)/*
