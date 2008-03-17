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
  FILE *infile = NULL;

  extern int optind;

  optind = 1;
  while (1) 
    {
      int c = getopt(argc, argv,"i:");
      if (c == -1)
      break;
      switch (c)
	{
	case 'i':
	  infile = fopen(optarg, "rb");
	  break;		  
	default:
	  break;
	}
      
    }

  if (infile == NULL)
    {
      //print_usage(argv);
      cout << "hey, needs some input files" << endl;
      exit(1);
    }
  
  size_t readsize = fread(sum, sizeof(double), CEL_INTENSITY_ARRAY_SIZE, infile);
  cout << "read " << readsize << " of " << CEL_INTENSITY_ARRAY_SIZE << endl;
  fclose(infile);

  for (int icel=0; icel< CEL_INTENSITY_ARRAY_SIZE; icel++)
    cout << icel << ": " << sum[icel] << endl;

  return 0;
}
