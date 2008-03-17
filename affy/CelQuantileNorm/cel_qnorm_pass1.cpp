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
//#include <cstdlib>
//#include <algorithm>
#include "CELFileData.h"
#include "stable_compare.h"
//2008-03-15 yh: CEL_INTENSITY_ARRAY_SIZE defined in common.h
#include "common.h"
//#define CEL_INTENSITY_ARRAY_SIZE 6553600

using namespace affxcel;
using namespace std;

int main(int argc, char **argv)
{
	float *intensity = (float *)calloc(CEL_INTENSITY_ARRAY_SIZE, sizeof(float));
	double *sum = (double *)calloc(CEL_INTENSITY_ARRAY_SIZE, sizeof(double));
	int read_count = 0;

	FILE *outfile = NULL;

	extern int optind;

	optind = 1;
	while (1)
	{
		int c = getopt(argc, argv,"o:");
		if (c == -1)
			break;
		switch (c)
		{
		case 'o':
			outfile = fopen(optarg, "wb");
			break;
		default:
			break;
		}

	}

	if ((optind >= argc) || (outfile == NULL))
	{
		//print_usage(argv);
		cout << "hey, needs some input and output files" << endl;
		cout << "e.g. runs with: program -o intfile *NSP* " << endl;
		exit(1);
	}

	while (optind < argc)
	{
		CCELFileData cel;
		cel.SetFileName(argv[optind]);

		cerr << "Processing " << cel.GetFileName() << endl;
		if (cel.Exists())
		{
			if(!cel.Read())
			{
				cerr << " read failed" << endl;
				exit(1);
			}

			if (cel.GetNumCells() != CEL_INTENSITY_ARRAY_SIZE)
			{
				cerr << "array sizes don't agree: " << cel.GetNumCells() << endl;
				exit(1);
			}

			for (int icel=0; icel< CEL_INTENSITY_ARRAY_SIZE; icel++)
			{
				intensity[icel] = cel.GetIntensity(icel);
			}
		}
		else
		{
			cerr << " doesn't exist" << endl;
			exit(1);
		}
		//2008/03/16 yh: sort the intensity from one array
		qsort(intensity, CEL_INTENSITY_ARRAY_SIZE, sizeof(float), stable_compare_float);
		//2008/03/16 yh: add the intensity of the same rank from each array to a sum, to take average later
		for (int icel=0; icel< CEL_INTENSITY_ARRAY_SIZE; icel++)
			sum[icel] += intensity[icel];

		int neg_count = 0;
		while (intensity[neg_count] <= 0.0)
		{
			neg_count++;
		}
		if (neg_count)
		{
			cerr << cel.GetFileName() << " contains " << neg_count << " non-positive probe intensities." << endl;
		}
		cerr << cel.GetFileName() << " low: " << intensity[0] << " , high: " << intensity[CEL_INTENSITY_ARRAY_SIZE-1] << endl;

		read_count++;
		cel.Close();

		optind++;
	}

	if (read_count)
	{
		//2008/03/16 yh: take average for each summed&ranked intensity
		for (int icel=0; icel< CEL_INTENSITY_ARRAY_SIZE; icel++)
			sum[icel] /= read_count;
		size_t writesize = fwrite(sum, sizeof(double), CEL_INTENSITY_ARRAY_SIZE, outfile);
		cout << "written " << writesize << " of " << CEL_INTENSITY_ARRAY_SIZE << endl;
	}
	fclose(outfile);

	return 0;
}

