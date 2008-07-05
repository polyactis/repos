/*
*
*2008-07-04
* given cel & cdf files, output the intensity of QC probes
*/


#include "OutputQCIntensity.h"

OutputQCIntensity::OutputQCIntensity()
{
	#if defined(DEBUG)
		std::cerr<<"OutputQCIntensity Begins"<<std::endl;
	#endif
}

OutputQCIntensity::OutputQCIntensity(char* _celFileName, char* _cdfFileName, char* _outputFileName)\
	:celFileName(_celFileName), cdfFileName(_cdfFileName), outputFileName(_outputFileName)
{
	#if defined(DEBUG)
		std::cerr<<"OutputQCIntensity Begins"<<std::endl;
		std::cerr << celFileName << '\t' << cdfFileName << std::endl;
	#endif
}



OutputQCIntensity::~OutputQCIntensity()
{
	#if defined(DEBUG)
		std::cerr<<"OutputQCIntensity Exits"<<std::endl;
	#endif
}

void OutputQCIntensity::_outputQCIntensity(char* celFileName, char* cdfFileName, char* outputFileName)
{
	#if defined(DEBUG)
		std::cerr<<"OutputQCIntensity _outputQCIntensity ..."<<std::endl;
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
	
	ofstream out;
	out.open(outputFileName);
	out << "QCProbeSet\tX\tY\tintensity\n";
	
	FusionCDFFileHeader cdfheader;
	cdfheader = cdf.GetHeader();
	
	std::string probeSetName;
	FusionCDFProbeSetInformation cdfProbeSetInfo;
	FusionCDFProbeGroupInformation group;
	FusionCDFProbeInformation probe;
	
	float intensity = 0;
	
	FusionCDFQCProbeSetInformation qcProbeSetInfo;
	FusionCDFQCProbeInformation qcProbeInfo;
	for (int i=0; i<cdfheader.GetNumQCProbeSets(); i++)	//2008-07-04 output QC probeset information
	{
		cdf.GetQCProbeSetInformation(i, qcProbeSetInfo);
		cout << "\tQCProbeSet: " << i << " has " << qcProbeSetInfo.GetNumCells() << " cells/probes." << endl;
		for (int icell=0; icell<qcProbeSetInfo.GetNumCells(); icell++)
		{
			qcProbeSetInfo.GetProbeInformation(icell, qcProbeInfo);
			/*
			cout <<"\t\tProbe: " << icell << " X=" << qcProbeInfo.GetX() \
			<< " Y="<<qcProbeInfo.GetY() << " Plen=" << qcProbeInfo.GetPLen() \
			<< " IsPerfectMatchProbe=" << qcProbeInfo.IsPerfectMatchProbe()\
			<< " IsBackgroundProbe=" << qcProbeInfo.IsBackgroundProbe() << endl;
			*/
			intensity = cel.GetIntensity(qcProbeInfo.GetX(), qcProbeInfo.GetY());
			out << i << "\t" << qcProbeInfo.GetX() << "\t" << qcProbeInfo.GetY() << "\t" << intensity << endl;
		}
	}
	
	#if defined(DEBUG)
		std::cerr<<"OutputQCIntensity _outputQCIntensity exits."<<std::endl;
	#endif
	out.close();
}


void OutputQCIntensity::run()
{
	_outputQCIntensity(celFileName, cdfFileName, outputFileName);
}

/*
BOOST_PYTHON_MODULE(OutputQCIntensity)
{
	boost::python::class_<OutputQCIntensity>("OutputQCIntensity")
		//.def("readCDF", &ReadatSNPtilgeno::readCDF)
		//.def("readCEL", &ReadatSNPtilgeno::readCEL)
		.def("_outputQCIntensity", &OutputQCIntensity::_outputQCIntensity)
	;
	def("set_module_and_type", &numeric::array::set_module_and_type);
}
*/

void print_usage(FILE* stream,int exit_code)
{
	assert(stream !=NULL);
        fprintf(stream,"Usage: OutputQCIntensity options -i celFileName -d cdfFileName\n");
	fprintf(stream,"\t-h  --help	Display the usage infomation.\n"\
		"\t-i ..., --input=...	celFileName\n"\
		"\t-d ...,	cdfFileName\n"\
		"\t-o ..., --output=...	Write output to file, OutputQCIntensity.output(default)\n"\
		"\tFor long option, = or ' '(blank) is same.\n"\
		"\tLine tokenizer is one tab\n");
	exit(3);
}

int main(int argc, char* argv[])
{
	int next_option;
	const char* const short_options="ho:i:d:e:p:";
	const struct option long_options[]={
	  {"help",0,NULL,'h'},
	  {"output",1,NULL,'o'},
	  {"input",1,NULL,'i'},
	  {NULL,0,NULL,0}
	};
	
	char* output_filename = "OutputQCIntensity.output";
	char* celFileName = NULL;
	char* cdfFileName = NULL;

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
		case '?':
			print_usage(stderr,-1);
		case -1:
			break;
		default:
			abort();
		}
	}while(next_option!=-1);

	if (celFileName!=NULL && cdfFileName!=NULL)
	{
		
		OutputQCIntensity instance(celFileName, cdfFileName, output_filename);
		instance.run();
		
	}
	else
		print_usage(stderr, 1);
}
