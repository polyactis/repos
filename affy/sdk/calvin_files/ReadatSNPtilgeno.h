/*
*
*2008-04-07
* a boost.python module to read atSNPtil_geno cel file and return intensity matrix.
*      
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

typedef boost::multi_array<float, 2> intensity_matrix_type;	//for an indigenous multi-dimension array besides a boost::python::numeric::array
//typedef intensity_matrix_type::index index;
	
class ReadatSNPtilgeno
{
	public:
		FusionCELData cel;
		FusionCDFData cdf;
		int ecotypeid;
		char* probeSubsetFname;
		char* celFileName;
		char* cdfFileName;
		char* outputFileName;
		
		std::map<std::string, int > snpID2index;
		std::map<std::string, std::string> snpID2alleles;
		//boost::python::numeric::array data_matrix;
		
		ReadatSNPtilgeno(char* _celFileName, char* _cdfFileName, char* _probeSubsetFname, char* _outputFileName, int _ecotypeid);
		ReadatSNPtilgeno(int _ecotypeid);
		~ReadatSNPtilgeno();
		void constructsnpID2index(char* probeSubsetFname);
		//boost::python::dict readCEL(char* celFileName, char* cdfFileName, boost::python::numeric::array& data_matrix);
		void _readCEL(char* celFileName, char* cdfFileName, intensity_matrix_type& data_matrix);
		//void outputCEL(char* outputFileName, boost::python::numeric::array& data_matrix, boost::python::dict& snpID2index, boost::python::dict& snpID2alleles);
		void _outputCEL(char* outputFileName, intensity_matrix_type& data_matrix, std::map<std::string, int >& snpID2index, std::map<std::string, std::string>& snpID2alleles);
		void run();
};