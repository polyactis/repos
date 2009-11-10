/*
	GADA v1.0 Genome Alteration Detection Algorithm
    Copyright (C) 2008  Childrens Hospital of Los Angeles

	GADA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    GADA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GADA.  If not, see <http://www.gnu.org/licenses/>.

	Author:
		Roger Pique-Regi    piquereg@usc.edu

*/

#include "BaseGADA.h"


// Global variables with algorithm parameters

double T=5.0; //Backward elimination threshold
//double T2=5.0; //Segment collapse to base Non-Alteration level threshold
double BaseAmp=0.0;  //Base-level
double a=0.2; //SBL hyperprior parameter
double sigma2=-1; //Variance observed, if negative value, it will be estimated by the mean of the differences
			  // I would recommend to be estimated on all the chromosomes and as a trimmed mean.
int MinLen=0; //Lenght in number of probes for a CNA segment to be called significan.
int SelectClassifySegments=0; //Classify segment into altered state (1), otherwise 0
int SelectEstimateBaseAmp=1; //Estimate Neutral hybridization amplitude.
char *InputFile;
char *OutputFile;


#ifndef GADABIN
class GADA
{
	public:
		int M;
		int i;
		int *Iext;
		int K;
		double *tn;
		int *SegLen;
		double *SegAmp;
		double *SegState;
		double *Wext;

		FILE *fin,*fout;

		double T; //Backward elimination threshold
		//double T2=5.0; //Segment collapse to base Non-Alteration level threshold
		double BaseAmp;  //Base-level
		double a; //SBL hyperprior parameter
		double sigma2; //Variance observed, if negative value, it will be estimated by the mean of the differences
					  // I would recommend to be estimated on all the chromosomes and as a trimmed mean.
		int MinLen; //Lenght in number of probes for a CNA segment to be called significan.
		int SelectClassifySegments; //Classify segment into altered state (1), otherwise 0
		int SelectEstimateBaseAmp; //Estimate Neutral hybridization amplitude.
		char *InputFile;
		char *OutputFile;

		string input_fname;
		string output_fname;

		GADA();
		~GADA();

		void readInIntensity(boost::python::list intensity_list);
		boost::python::list run(boost::python::list intensity_list, double aAlpha, double TBackElim, int MinSegLen);
		void cleanupMemory();
};

GADA::GADA()
{
	T=5.0; //Backward elimination threshold
	//double T2=5.0; //Segment collapse to base Non-Alteration level threshold
	BaseAmp=0.0;  //Base-level
	a=0.2; //SBL hyperprior parameter
	sigma2=-1; //Variance observed, if negative value, it will be estimated by the mean of the differences
				  // I would recommend to be estimated on all the chromosomes and as a trimmed mean.
	MinLen=0; //Lenght in number of probes for a CNA segment to be called significan.
	SelectClassifySegments=0; //Classify segment into altered state (1), otherwise 0
	SelectEstimateBaseAmp=1; //Estimate Neutral hybridization amplitude.
}

GADA::~GADA()
{

}
void GADA::readInIntensity(boost::python::list intensity_list)
{
#if defined(DEBUG)
	cerr<< boost::format("# Start reading ... \n") ;
#endif
	int no_of_probes = boost::python::extract<int>(intensity_list.attr("__len__")());
	M = no_of_probes;
	tn= (double *) calloc(no_of_probes, sizeof(double));

	for(int i=0; i<no_of_probes-1;i++)
		tn[i] = boost::python::extract<double>(intensity_list[i]);

#if defined(DEBUG)
	cerr<< boost::format("# M=%1% probes in input file\n") % M ;
#endif
}

void GADA::cleanupMemory()
{
	free(tn);
	//free(Iext);
	free(SegLen);
	free(SegAmp);
	//free(SegState);
	//free(Wext);
}

boost::python::list GADA::run(boost::python::list intensity_list, double aAlpha, double TBackElim, int MinSegLen)
{
	readInIntensity(intensity_list);


	K = SBLandBE(tn, M, &sigma2, aAlpha, 0, 0, &Iext, &Wext);

#if defined(DEBUG)
	//K=SBLandBE(tn,M,&sigma2,a,T,MinLen,&Iext,&w);
    cerr<< boost::format("# Overall mean %1%.\n")%Wext[0];
	cerr<< boost::format("# Sigma^2=%1%.\n")%sigma2;
	cerr<< boost::format("# Found %1% breakpoints after SBL\n")%K;
#endif

	BEwTandMinLen(Wext,Iext,&K, sigma2, TBackElim, MinSegLen);

#if defined(DEBUG)
    cerr<< boost::format("# Kept %1% breakpoints after BE\n")%K;
#endif

	SegLen= (int*) calloc(K+1,sizeof(int));
	SegAmp= (double *) calloc(K+1,sizeof(double));
	IextToSegLen(Iext,SegLen,K);
	IextWextToSegAmp(Iext,Wext,SegAmp,K);
#if defined(DEBUG)
    cerr<< boost::format("# Making segments\n");
#endif

    boost::python::list return_ls;

	//fprintf(fout,"Start\tStop\tLength\tAmpl\n");
	for(i=0;i<K+1;i++)
	{
		boost::python::list d_row;
		d_row.append(Iext[i]+1);
		d_row.append(Iext[i+1]);
		d_row.append(SegLen[i]);
		d_row.append(SegAmp[i]);
		//cerr<< boost::format("%1% \t %2% \t %3% %4% \n")%Iext[i] % Iext[i+1] % SegLen[i] % SegAmp[i];
		return_ls.append(d_row);
	}

	cleanupMemory();
	return return_ls;

}

BOOST_PYTHON_MODULE(GADA)
{
	using namespace boost::python;
	class_<GADA>("GADA")
		.def("run", &GADA::run)
	;

}

#endif

void help_message(FILE *fd)
{
	//fprintf(fd,"# Welcome to GADA 1.0 \n");
	fprintf(fd,"# Usage:\n");
	fprintf(fd,"# \t GADA [-T 5] [-a 0.8] [-s -1] [-M 3] -i input.txt -o output.txt\n");
	fprintf(fd,"# \t input.txt is a single column text file with no header\n");
	fprintf(fd,"# Possible options:\n");
	fprintf(fd,"# \t -i\t Input file. Otherwise standard input assumed\n");
	fprintf(fd,"# \t -o\t Output file. Otherwise standard output assumed\n");
	fprintf(fd,"# \t -a\t is the SBL hyperprior parameter for a breakpoint. It is the \n");
	fprintf(fd,"# \t\t shape parameter of the Gamma distribution. Default value %g.\n",a);
	fprintf(fd,"# \t\t Higher (lower) value more (less) breakpoints expected a priori\n");
	fprintf(fd,"# \t -T\t is the backward elimination critical value for a breakpoint. \n");
	fprintf(fd,"# \t\t i.e. the minimum difference between the segment means divided\n");
	fprintf(fd,"# \t\t by sigma. The default value for T is %g\n",T);
	fprintf(fd,"# \t -M\t is the minimum size in number of probes for a segment to be \n");
	fprintf(fd,"# \t\t deemed significant. Default value is %d.\n",MinLen);
	fprintf(fd,"# \t -s\t The variance estimate for the noise. If not provided or \n");
	fprintf(fd,"# \t\t negative it will be estimated from the provided data. \n");
	fprintf(fd,"# \t\t We recomend to estimate this on the entire data and not \n");
    fprintf(fd,"# \t\t separately for each chromosome.\n");
	fprintf(fd,"# \t -c\t Classify segments into altered state (L)oss -- (N)eutral \n");
	fprintf(fd,"# \t\t (G)ain). If c option is not specified, segments are returned \n");
	fprintf(fd,"# \t\t with their mean.\n");
	fprintf(fd,"# \t -b\t Mean amplitude associated to the Neutral state. If not \n");
	fprintf(fd,"# \t\t provided, and c option is used, then it is estimated as the\n");
	fprintf(fd,"# \t\t median value of all probes hybridization values after running \n");
	fprintf(fd,"# \t\t the algorithm. We recomend to estimate this on chromosomes\n");
	fprintf(fd,"# \t\t that are known to have a Neutral state on most areas. In some\n");
	fprintf(fd,"# \t\t cases this value may be known if we have been applied some \n");
	fprintf(fd,"# \t\t normalization, preprocessing or using another sample as ref.\n");
	fprintf(fd,"# \t -h \t Prints this help message.\n");
}

void help_and_exit(FILE *fd,int code)
{
	fprintf(fd, "Invalid syntax. Use GADA -h for help\n");
	exit(code);
}


void Configure(int ac, char *av[], string &input_fname, string &output_fname)
{
#if defined(DEBUG)
    cerr<< boost::format("# Getting commandline arguments ...");
#endif
    int CLcount,i;
	FILE *fd;

	// Parse the command line
	CLcount=1;


	fd=stdout;
	fprintf(fd,"# NumArgs = %d \n",ac);
	fprintf(fd,"# CallStr = ");

	for(i=0;i<ac;i++)
		fprintf(fd,"%s ",av[i]);
	fprintf(fd,"\n# Parsing Arguments: \n");


	while (CLcount < ac)
	{
		if (0 == strncmp (av[CLcount], "-h", 2))
		{
			help_message(stderr);
			exit(0);
		}

		else if (0 == strncmp (av[CLcount], "-c", 2))
		{
			SelectClassifySegments = 1;
			CLcount+=1;
			fprintf(fd,"# -c option activated to classify segments into altered states\n");
		}
		else if (0 == strncmp (av[CLcount], "-i", 2)) //! Input file
		{
			if(0 == strncmp (av[CLcount+1], "-", 1))
				help_and_exit(stderr,1);
			InputFile = (char*)malloc(strlen(av[CLcount+1]));
			//strcpy(InputFile, av[CLcount+1]);
			input_fname = av[CLcount+1];
			cout << boost::format("# Input file: %1%")%input_fname << std::endl;
			InputFile = av[CLcount+1];
			CLcount += 2;
			//fprintf(fd,"# Input file: %s \n",InputFile);
		}
		else if (0 == strncmp (av[CLcount], "-o", 2)) //! Output File
		{
			if(0 == strncmp (av[CLcount+1], "-", 1))
				help_and_exit(stderr,1);
			OutputFile = (char*)malloc(strlen(av[CLcount+1]));
			strcpy(OutputFile,av[CLcount+1]);
			output_fname = av[CLcount+1];
			//OutputFile = av[CLcount+1];
			CLcount += 2;
			fprintf(fd,"# Output file: %s \n",OutputFile);
		}
		else if (0 == strncmp (av[CLcount], "-a", 2)) //! a parameter
		{
			if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
			sscanf (av[CLcount+1], "%lf", &a);
			CLcount += 2;
			fprintf(fd,"# a= %g \n",a);
		}
		else if (0 == strncmp (av[CLcount], "-T", 2)) //! T parameter
		{
			if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
			sscanf (av[CLcount+1], "%lf", &T);
			CLcount += 2;
			fprintf(fd,"# T= %g \n",T);
		}
		else if (0 == strncmp (av[CLcount], "-M", 2)) //! MinLen parameter
		{
			if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
			sscanf (av[CLcount+1], "%d", &MinLen);
			CLcount += 2;
			fprintf(fd,"# MinLen= %d \n",MinLen);
		}
		else if (0 == strncmp (av[CLcount], "-b", 2)) //! BaseAmp parameter
		{
			//if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
			sscanf (av[CLcount+1], "%lf", &BaseAmp);
			CLcount += 2;
			fprintf(fd,"# BaseAmp= %g \n",BaseAmp);
			SelectEstimateBaseAmp=0;
		}
		else if (0 == strncmp (av[CLcount], "-s", 3)) //! sigma2 parameter
		{
			//if(0 == strncmp (av[CLcount+1], "-", 1))help_and_exit(stderr,1);
			sscanf (av[CLcount+1], "%lf", &sigma2);
			CLcount += 2;
			fprintf(fd,"# sigma2= %g \n",sigma2);
		}
		else
		{
			help_and_exit(stderr,1);
		}
	}
#if defined(DEBUG)
	cerr<< boost::format("Done.\n");
#endif
}

int main(int argc, char *argv[])
{
#if defined(DEBUG)
	cerr<< boost::format("# Starting ....\n");
#endif
	int M=1000;
	int i;
	int *Iext;
	int K;
	double *tn;
	int *SegLen;
	double *SegAmp;
	double *SegState;
	double *Wext;

	FILE *fin,*fout;

	string input_fname;
	string output_fname;

	tn= (double *) calloc(M,sizeof(double));

	Configure(argc,argv, input_fname, output_fname);

	// fin = stdin;
	fin = fopen(InputFile, "r");
	ifstream in;
	in.open(InputFile);
	// fout = stdout;
	fout = fopen(OutputFile, "w");
	//ofstream out;
	//out.open(output_fname.c_str());

	fprintf(stdout,"# GADA v1.0 Genome Alteration Detection Algorithm\n");
	fprintf(stdout,"# Copyright (C) 2008  Childrens Hospital of Los Angeles\n");
	fprintf(stdout,"# author: Roger Pique-Regi piquereg@usc.edu\n");


	cout << boost::format("# Parameter setting: a=%1%,T=%2%,MinLen=%3%,sigma2=%4%,BaseAmp=%5%.") % a % T % MinLen % sigma2 % BaseAmp << std::endl;
	//fprintf(fout,"#Parameter setting: a=%g,T=%g,MinLen=%d,sigma2=%g,BaseAmp=%g\n",a,T,MinLen,sigma2,BaseAmp);

	/*
	#if defined(DEBUG)
		std::cerr<<"Read in the data...";
	#endif
	string line;
	while(getline(in, line))
	{
		cout << line[0];
	}
	 */
	i=0;
	while(!feof(fin))
	{
		fscanf(fin,"%lf",&tn[i++]);
		if(i>=M){
			M=M+1000;
			tn=(double *) realloc(tn,M*sizeof(double));
		}
	}
	M=i-1;
	tn= (double*) realloc(tn,M*sizeof(double));

	fprintf(fout,"# Reading M=%d probes in input file\n",M);

	K=SBLandBE(tn,M,&sigma2,a,0,0,&Iext,&Wext);

	//K=SBLandBE(tn,M,&sigma2,a,T,MinLen,&Iext,&w);
    fprintf(fout,"# Overall mean %g\n",Wext[0]);
	fprintf(fout,"# Sigma^2=%g\n",sigma2);
	fprintf(fout,"# Found %d breakpoints after SBL\n",K);


	BEwTandMinLen(Wext,Iext,&K,sigma2,T,MinLen);
    fprintf(fout,"# Kept %d breakpoints after BE\n",K);

	SegLen= (int*) calloc(K+1,sizeof(int));
	SegAmp= (double *) calloc(K+1,sizeof(double));
	IextToSegLen(Iext,SegLen,K);
	IextWextToSegAmp(Iext,Wext,SegAmp,K);
    fprintf(fout,"# Making segments\n");

	//Collapse Segments
	if(SelectClassifySegments==1)
	{
		if(SelectEstimateBaseAmp==1)
		{
			BaseAmp=CompBaseAmpMedianMethod(SegLen,SegAmp,K);
			fprintf(fout,"# Estimating BaseAmp\n");
		}
		fprintf(fout,"# BaseAmp=%g \n",BaseAmp);
		fprintf(fout,"# Classify Segments \n",BaseAmp);


		SegState= (double *) calloc(K+1,sizeof(double));
		for(i=0;i<=K;i++)SegState[i]=SegAmp[i];
		CollapseAmpTtest(SegState,SegLen,K,BaseAmp,sigma2,T);
	}


	if(SelectClassifySegments==0)
	{
    	fprintf(fout,"Start\tStop\tLength\tAmpl\n");
		for(i=0;i<K+1;i++)
			fprintf(fout,"%d\t%d\t%d\t%g\n",Iext[i]+1,Iext[i+1],SegLen[i],SegAmp[i]);
	}
	else if(SelectClassifySegments==1)
	{
	   	fprintf(fout,"Start\tStop\tLenght\tAmpl\tState\n");
		for(i=0;i<K+1;i++)
		{
			fprintf(fout,"%d\t%d\t%d\t%g\t",Iext[i]+1,Iext[i+1],SegLen[i],SegAmp[i]);
			if(SegState[i]>BaseAmp)
				fprintf(fout,"G");
			else if(SegState[i]<BaseAmp)
				fprintf(fout,"L");
			else
				fprintf(fout,"N");
			fprintf(fout,"\n");
		}


	}




}
