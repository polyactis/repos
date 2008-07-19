/*
*
*2008-04-07
* a boost.python module to read atSNPtil_geno cel file and return intensity matrix.
*/


#include "ReadatSNPtilgeno.h"

ReadatSNPtilgeno::ReadatSNPtilgeno(int _ecotypeid)\
	:ecotypeid(_ecotypeid)
{
	#if defined(DEBUG)
		std::cerr<<"ReadatSNPtilgeno Begins"<<std::endl;
		std::cerr<<ecotypeid<<std::endl;
	#endif
}

ReadatSNPtilgeno::ReadatSNPtilgeno(char* _celFileName, char* _cdfFileName, char* _probeSubsetFname, char* _outputFileName, int _ecotypeid, int _debug)\
	:celFileName(_celFileName), cdfFileName(_cdfFileName), probeSubsetFname(_probeSubsetFname), outputFileName(_outputFileName), ecotypeid(_ecotypeid), debug(_debug)
{
	#if defined(DEBUG)
		std::cerr<<"ReadatSNPtilgeno Begins"<<std::endl;
		std::cerr << celFileName << '\t' << cdfFileName << '\t' << probeSubsetFname << '\t' << ecotypeid << std::endl;
	#endif
}



ReadatSNPtilgeno::~ReadatSNPtilgeno()
{
	#if defined(DEBUG)
		std::cerr<<"ReadatSNPtilgeno Exits"<<std::endl;
	#endif
}


void ReadatSNPtilgeno::constructsnpID2index(char* probeSubsetFname)
{
	#if defined(DEBUG)
		std::cerr<<"ReadatSNPtilgeno constructsnpID2index ..."<<std::endl;
	#endif
	std::ifstream probeSubsetF(probeSubsetFname);
	//boost::char_separator<char> sep(" \t()");		//blank, '\t' or '(' or ')' is the separator
	//typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	int row_index = 0;
	for (std::string line; std::getline(probeSubsetF, line);)
	{
		/*
		tokenizer line_toks(line, sep);
		for (tokenizer::iterator tokenizer_iter = line_toks.begin(); tokenizer_iter!=line_toks.end();++tokenizer_iter)
		{
			#if defined(DEBUG0)
				std::cerr<<*tokenizer_iter<<'\t';
			#endif
			signature_frequency_vector.push_back(atoi((*tokenizer_iter).c_str()));
		}
		
		#if defined(DEBUG)
			std::cerr<<line<< "\t" << row_index << "\t" << line.c_str() << std::endl;
		#endif
		*/
		snpID2index[line] = row_index;
		row_index++;
		/*
		#if defined(DEBUG)
			std::cerr<<line<< "\t" << row_index << std::endl;
		#endif
		*/
	}
	probeSubsetF.close();
	#if defined(DEBUG)
		std::cerr<<"ReadatSNPtilgeno constructsnpID2index exits."<<std::endl;
	#endif
}
/*
boost::python::dict ReadatSNPtilgeno::readCEL(char* celFileName, char* cdfFileName, boost::python::numeric::array& data_matrix)
{
	cel.SetFileName(celFileName);
	if (cel.Read() == false)
	{
		cout << "Failed to read the file." << endl;
		exit(4);
	}
	
	cdf.SetFileName(cdfFileName);
	if (cdf.Read() == false)
	{
		std::cerr << "Failed to read the CDF file." << endl;
		exit(4);
	}
	
	FusionCDFFileHeader cdfheader;
	cdfheader = cdf.GetHeader();
	
	boost::python::dict snpID2alleles;
	std::string probeSetName;
	FusionCDFProbeSetInformation cdfProbeSetInfo;
	FusionCDFProbeGroupInformation group;
	FusionCDFProbeInformation probe;
	float intensity = 0;
	int row_index;
	for (int i=0; i<cdfheader.GetNumProbeSets(); i++)
	{
		probeSetName = cdf.GetProbeSetName(i);
		//probeSetType = cdf.GetProbeSetType(i);
		cdf.GetProbeSetInformation(i, cdfProbeSetInfo);
		if (snpID2index.has_key(probeSetName))
		{
			row_index = snpID2index[probeSetName];
			std::string alleles("");
			int ngroups = cdfProbeSetInfo.GetNumGroups();
			for (int igroup=0; igroup<ngroups; igroup++)
			{
				cdfProbeSetInfo.GetGroupInformation(igroup, group);
				alleles += group.GetName();
				group.GetCell(0, probe);
				intensity = cel.GetIntensity(probe.GetX(), probe.GetY());
				data_matrix[boost::make_tuple(row_index,igroup)] = intensity;
			}
			snpID2alleles[probeSetName] = alleles;
		}
	}
	return snpID2alleles;
}
*/

void ReadatSNPtilgeno::_readCEL(char* celFileName, char* cdfFileName, intensity_matrix_type& data_matrix)
{
	#if defined(DEBUG)
		std::cerr<<"ReadatSNPtilgeno _readCEL ..."<<std::endl;
	#endif
	cel.SetFileName(celFileName);
	if (cel.Read() == false)
	{
		std::cerr << "Failed to read the file." << endl;
		exit(4);
	}
	
	cdf.SetFileName(cdfFileName);
	if (cdf.Read() == false)
	{
		std::cerr << "Failed to read the CDF file." << endl;
		exit(4);
	}
	
	FusionCDFFileHeader cdfheader;
	cdfheader = cdf.GetHeader();
	
	std::map<std::string, int >::iterator snpID2index_it;
	std::string probeSetName;
	FusionCDFProbeSetInformation cdfProbeSetInfo;
	FusionCDFProbeGroupInformation group;
	FusionCDFProbeInformation probe;
	float intensity = 0;
	int row_index;
	for (int i=0; i<cdfheader.GetNumProbeSets(); i++)
	{
		probeSetName = cdf.GetProbeSetName(i);
		//probeSetType = cdf.GetProbeSetType(i);
		cdf.GetProbeSetInformation(i, cdfProbeSetInfo);
		snpID2index_it = snpID2index.find(probeSetName);
		if (snpID2index_it!=snpID2index.end())
		{
			row_index = snpID2index[probeSetName];
			std::string alleles("");
			int ngroups = cdfProbeSetInfo.GetNumGroups();
			for (int igroup=0; igroup<ngroups; igroup++)
			{
				cdfProbeSetInfo.GetGroupInformation(igroup, group);
				alleles += group.GetName();
				if (debug!=0)
					std::cerr<< probeSetName << " group " << igroup << " " << group.GetName() << " " << alleles <<std::endl;
				group.GetCell(0, probe);
				intensity = cel.GetIntensity(probe.GetX(), probe.GetY());
				data_matrix[row_index][igroup] = intensity;
			}
			snpID2alleles[probeSetName] = alleles;
			if (debug!=0)
				std::cerr<< (probeSetName==(*snpID2index_it).first) <<std::endl;
		}
	}
	#if defined(DEBUG)
		std::cerr<<"ReadatSNPtilgeno _readCEL exits."<<std::endl;
	#endif
}

/*
void ReadatSNPtilgeno::outputCEL(char* outputFileName, boost::python::numeric::array& data_matrix, boost::python::dict& snpID2index, boost::python::dict& snpID2alleles)
{
	ofstream out;
	out.open(outputFileName);
	boost::python::list snpIDkeys = snpID2alleles.keys();
	std::string snpID;
	boost::python::tuple shape = boost::python::extract<boost::python::tuple>(data_matrix.attr("shape"));
	int numberOfRows;
	int numberOfColumns;
	numberOfRows = boost::python::extract<int>(shape[0]);
	if (numberOfRows == 0)
	{
		cerr<<"No data"<<endl;
		numberOfColumns = 0;
	}
	else
		numberOfColumns = boost::python::extract<int>(shape[1]);
	
	out << "SNP ID\t" << "sense1" << "\t" <<"sense2" << "\t" << "antisense1\t" << "antisense2\t" << std::endl; 
	int row_index;
	for (int i = 0; i < len(snpIDkeys); ++i)
	{
		snpID = boost::python::extract<std::string>(snpIDkeys[i]);
		row_index = boost::python::extract<int>(snpID2index[snpID]);
		snpID += "_";
		snpID += boost::python::extract<std::string>(snpID2alleles[snpID]);
		out << snpID;
		for (int j=0; j<numberOfColumns; j++)
		{
			out << "\t";
			out << boost::python::extract<float>(data_matrix[row_index][j]);
		}
		out << std::endl;
	}
	out.close();
}
*/

void ReadatSNPtilgeno::_outputCEL(char* outputFileName, intensity_matrix_type& data_matrix, \
	std::map<std::string, int >& snpID2index, std::map<std::string, std::string>& snpID2alleles)
{
	#if defined(DEBUG)
		std::cerr<<"ReadatSNPtilgeno _outputCEL ... " << snpID2alleles.size() << " snps in snpID2alleles."<<std::endl;
	#endif
	ofstream out;
	out.open(outputFileName);
	std::string snpID;
	//boost::python::tuple shape = boost::python::extract<tuple>(data_matrix.attr("shape"));
	int numberOfRows;
	int numberOfColumns;
	numberOfRows = data_matrix.shape()[0];
	if (numberOfRows == 0)
	{
		cerr<<"No data"<<endl;
		numberOfColumns = 0;
	}
	else
		numberOfColumns = data_matrix.shape()[1];
	
	out << "SNP ID\t" << ecotypeid << "_sense1\t" << ecotypeid << "_sense2\t" << ecotypeid << "_antisense1\t" << ecotypeid << "_antisense2" << std::endl; 
	int row_index;
	std::map<std::string, int >::iterator snpID2index_it;
	std::map<std::string, std::string >::iterator snpID2alleles_it;
	for ( snpID2index_it = snpID2index.begin(); snpID2index_it != snpID2index.end(); ++snpID2index_it)
	{
		snpID = (*snpID2index_it).first;
		row_index = (*snpID2index_it).second;
		snpID2alleles_it = snpID2alleles.find(snpID);
		snpID += "_";
		if (snpID2alleles_it!=snpID2alleles.end())
			snpID += (*snpID2alleles_it).second;
		if (debug!=0)
			std::cerr<< (*snpID2alleles_it).second <<std::endl;
		out << snpID;
		for (int j=0; j<numberOfColumns; j++)
		{
			out << "\t";
			out << data_matrix[row_index][j];
		}
		out << std::endl;
	}
	#if defined(DEBUG)
		std::cerr<<"ReadatSNPtilgeno _outputCEL exits."<<std::endl;
	#endif
	out.close();
}
void ReadatSNPtilgeno::run()
{
	constructsnpID2index(probeSubsetFname);
	
	int no_of_snps = snpID2index.size();
	intensity_matrix_type data_matrix(boost::extents[no_of_snps][4]);
	_readCEL(celFileName, cdfFileName, data_matrix);
	_outputCEL(outputFileName, data_matrix, snpID2index, snpID2alleles);
}


BOOST_PYTHON_MODULE(ReadatSNPtilgeno)
{
	boost::python::class_<ReadatSNPtilgeno>("ReadatSNPtilgeno", init<int>())
		//.def("readCDF", &ReadatSNPtilgeno::readCDF)
		//.def("readCEL", &ReadatSNPtilgeno::readCEL)
		.def("constructsnpID2index", &ReadatSNPtilgeno::constructsnpID2index)
		//.def("outputCEL", &ReadatSNPtilgeno::outputCEL)
		.def_readonly("ecotypeid", &ReadatSNPtilgeno::ecotypeid)
		.def_readonly("snpID2index", &ReadatSNPtilgeno::snpID2index)
	;
	def("set_module_and_type", &numeric::array::set_module_and_type);
}

void print_usage(FILE* stream,int exit_code)
{
	assert(stream !=NULL);
    fprintf(stream,"Usage: ReadatSNPtilgeno options -i celFileName -d cdfFileName -e ecotypeid -p probeSubsetFname\n");
	fprintf(stream,"\t-h  --help	Display the usage infomation.\n"\
		"\t-o ..., --output=...	Write output to file, ReadatSNPtilgeno.output(default)\n"\
		"\t-i ..., --input=...	celFileName\n"\
		"\t-d ...,	cdfFileName\n"\
		"\t-e ...,	ecotypeid\n"\
		"\t-p ...,	probeSubsetFname\n"\
		"\t-b,	toggle debug mode, more output\n"\
		"\tFor long option, = or ' '(blank) is same.\n"\
		"\tLine tokenizer is one tab\n");
	fprintf(stream, "Examples:\n");
	fprintf(stream, "\tReadatSNPtilgeno -i /Network/Data/250k/db/raw_data/160_raw_data.cel -d /Network/Data/250k/raw_data/atSNPtil_geno.cdf -e 1 -p ~/script/affy/250k_test/250kprobe_subset.txt  -o /tmp/160.tsv\n");
	exit(3);
}

int main(int argc, char* argv[])
{
	int next_option;
	const char* const short_options="ho:i:d:e:p:b";
	const struct option long_options[]={
	  {"help",0,NULL,'h'},
	  {"output",1,NULL,'o'},
	  {"input",1,NULL,'i'},
	  {NULL,0,NULL,0}
	};
	
	char* output_filename = "ReadatSNPtilgeno.output";
	char* celFileName = NULL;
	char* cdfFileName = NULL;
	char* probeSubsetFname = NULL;
	int ecotypeid = -1;
	int debug = 0;

	do
	{
		next_option=getopt_long(argc,argv,short_options,long_options,NULL);
		switch(next_option)
		{
		case 'h':
			print_usage(stdout,0);
	  		exit(1);
		case 'o':
			output_filename=optarg;
			break;
		case 'i':
			celFileName = optarg;
			break;
		case 'd':
			cdfFileName = optarg;
			break;
		case 'e':
			ecotypeid = atoi(optarg);
			break;
		case 'p':
			probeSubsetFname = optarg;
			break;
		case 'b':
			debug = 1;
			break;
		case '?':
			print_usage(stderr,-1);
		case -1:
			break;
		default:
			abort();
		}
	}while(next_option!=-1);

	if (celFileName!=NULL && cdfFileName!=NULL && ecotypeid!=-1 && probeSubsetFname!=NULL)
	{
		
		ReadatSNPtilgeno instance(celFileName, cdfFileName, probeSubsetFname, output_filename, ecotypeid, debug);
		instance.run();
		
	}
	else
		print_usage(stderr, 1);
}
