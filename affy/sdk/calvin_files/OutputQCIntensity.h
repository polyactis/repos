#ifndef OUTPUTQCINTENSITY_H_
#define OUTPUTQCINTENSITY_H_

#endif /*OUTPUTQCINTENSITY_H_*/

/*
*
*2008-07-04
* given cel & cdf files, output the intensity of QC probes
*/
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <boost/tokenizer.hpp>	//for tokenizer, parse input file
#include <fstream>
#include <boost/array.hpp>
#include <string>
#include <boost/python.hpp>

#include <map>		//for std::map
#include <boost/utility.hpp>             // for boost::tie
#include <boost/multi_array.hpp>	//for boost::multi_array

#include <getopt.h>	//to parse program options.

#include "fusion/src/FusionCELData.h"
#include "fusion/src/FusionCDFData.h"
/*#include "file/CDFFileData.h"*/

using namespace affymetrix_fusion_io;

using namespace boost;
using namespace boost::python;
using namespace std;

//typedef boost::multi_array<float, 2> intensity_matrix_type;	//for an indigenous multi-dimension array besides a boost::python::numeric::array
//typedef intensity_matrix_type::index index;
	
class OutputQCIntensity
{
	public:
		FusionCELData cel;
		FusionCDFData cdf;
		char* celFileName;
		char* cdfFileName;
		char* outputFileName;
		
		OutputQCIntensity();
		OutputQCIntensity(char* _celFileName, char* _cdfFileName, char* _outputFileName);
		~OutputQCIntensity();
		void _outputQCIntensity(char* celFileName, char* cdfFileName, char* outputFileName);
		void run();
};
