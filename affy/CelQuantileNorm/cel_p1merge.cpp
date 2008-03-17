/**
 * Copyright (c) 2006 Hin-Tak Leung
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 **/

#include <iostream>

//2008-03-15 yh: CEL_INTENSITY_ARRAY_SIZE defined in common.h
#include "common.h"
//#define CEL_INTENSITY_ARRAY_SIZE 6553600

using namespace std;

int main(int argc, char **argv)
{
  double *sum = (double *)calloc(CEL_INTENSITY_ARRAY_SIZE, sizeof(double));
  double *current = (double *)calloc(CEL_INTENSITY_ARRAY_SIZE, sizeof(double));
  FILE *descfile = NULL; 
  FILE *outfile = NULL; 
  int writehelp = 0;
  int sumsample = 0;

  extern int optind;

  optind = 1;
  while (1) 
    {
      int c = getopt(argc, argv,"i:o:h");
      if (c == -1)
      break;
      switch (c)
	{
	case 'i':
	  descfile = fopen(optarg, "rb");
	  break;		  
	case 'o':
	  outfile = fopen(optarg, "wb");
	  break;
	case 'h':
	  writehelp = 1;
	  break;		  
	default:
	  break;
	}
      
    }

  if ((descfile == NULL) || (outfile == NULL) || writehelp)
    {
      cout << "Usage:" << endl;
      cout << "     cel_p1merge -o outfile -i desfile" << endl;
      exit(1);
    }

  char filename[1024];
  int weight;
  while (fscanf(descfile, "%s %i",filename, &weight) == 2) {
    printf("reading %s with weight %i\n", filename, weight);
    FILE *infile = fopen(filename, "rb");    
    size_t readsize = fread(current, sizeof(double), CEL_INTENSITY_ARRAY_SIZE, infile);
    cout << "read " << readsize << " of " << CEL_INTENSITY_ARRAY_SIZE << endl;
    fclose(infile);
    sumsample += weight;
    for (int i = 0; i < CEL_INTENSITY_ARRAY_SIZE ; i++) {
      sum[i] += current[i] * weight; 
    }
  }

  for (int i = 0; i < CEL_INTENSITY_ARRAY_SIZE ; i++) {
    sum[i] /= sumsample; 
  }

  size_t writesize = fwrite(sum, sizeof(double), CEL_INTENSITY_ARRAY_SIZE, outfile);
  cout << "written " << writesize << " of " << CEL_INTENSITY_ARRAY_SIZE << endl;
  fclose(outfile);
  return 0;
}
