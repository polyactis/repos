#
# Copyright (c) 2006 Hin-Tak Leung
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use, copy,
# modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

EXE = cel_qnorm_pass1 cel_p1readtest gtype_cel_to_pq cel_p1merge

THE_REST = cdfimage cdfimage+ cdfimage+snp celimage

BOOST_OBJS = file_descriptor.o mapped_file.o zlib.o bzip2.o 

CXXFLAGS ?=
LDFLAGS ?=
CXXFLAGS += -I. -DHAVE_BOOST -O2 -Wall -pipe -fno-strict-aliasing
LDFLAGS += -lz -lbz2


.PHONY : default clean

default: $(EXE)

all: $(EXE) $(THE_REST)

cel_qnorm_pass1 : cel_qnorm_pass1.o stable_compare.o CELFileData.o FileIO.o $(BOOST_OBJS) 
	$(CXX) -o $@ $^ $(LDFLAGS) 

cel_p1readtest : cel_p1readtest.o
	$(CXX) -o $@ $^ $(LDFLAGS)

cel_p1merge : cel_p1merge.o
	$(CXX) -o $@ $^ $(LDFLAGS)

gtype_cel_to_pq : gtype_cel_to_pq.o CDFFileData.o stable_compare_rank.o hash.o log_average.o CELFileData.o FileIO.o $(BOOST_OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

cdfimage : cdfimage.o CDFFileData.o FileIO.o 
	$(CXX) -o $@ $^ $(LDFLAGS)

cdfimage+ : cdfimage+cel.o CDFFileData.o CELFileData.o FileIO.o $(BOOST_OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

cdfimage+snp : cdfimage+snp.o hash.o CDFFileData.o FileIO.o $(BOOST_OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

celimage : celimage.o stable_compare_rank.o stable_compare_rank.o CELFileData.o FileIO.o $(BOOST_OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean :
	$(RM) *.o *~ $(EXE) $(THE_REST)

#.ps.pdf is the old fashion way

%.pdf: %.ps
	ps2pdf -sPAPERSIZE=a4 -dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode -dPDFsettings=/prepress $^ $@

%.ps: %.dvi
	dvips -Ppdf -G0 $^ -o $@

%.dvi: %.tex
	latex $^
	latex $^
