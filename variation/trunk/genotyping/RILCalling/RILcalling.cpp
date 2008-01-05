// This program is a base-calling algorithm for recombinant inbred lines in
// an Arabidopsis context
#include<cmath>
// NB. Bases are coded as follows:
// 0 = A
// 1 = C
// 2 = G
// 3 = T
// or 0 (reference) and 1 (mutant) for SFPs/SNPs
 
// options
const bool bEstimateLocations=0; // =1 to estimate locations, =0 to keep them fixed
const bool bJohnsFormat=1;
bool bEstLocations=bEstimateLocations;
const int iNoSNPs=10; // the number of SNPs
const int iNoLines=200;
int iNoOfAlleles=4;
int iNoLinesTemp=iNoLines;
int iNoSNPsTemp=iNoSNPs;
int **giIncidenceMx; // stores incidence of SNPs when we are estimating locations.
double **gdDistanceMx; // stores the pairwise distance between loci

#include "MCMCFunction.h"
#include <iostream>
#include <algorithm>
using namespace std;
#include<string>
#include<fstream>
#include<math.h>
#include <cassert>
#include <time.h>
#include <iostream.h>
//#include <stdio.h>
#include <malloc.h>
#include <memory>
#include<stdlib.h>
#include<search.h>

ofstream ofile;
ifstream infile1("input.txt");// | ios::nocreate);

//using namespace ranlib; 

char **cRILarray,*cParentArray1,*cParentArray2; // hold the SNP info
double ***RILprior,***RILposterior;
int *iSNPIndices,*giOldSNPIndices;//,*iCopy; // stores the indices of the SNPs, if known
double *gdSNPLocations,*gdOldSNPLocations; // stores the locations of the SNPs along the chromosome 
int giBase; // stores the identity of the base which has its state changed in the MCMC proposal
int giLine; // stores the identity of the line which has its state changed in the MCMC proposal

int NUMBER_OF_MCMC_ITERATIONS=10000;//000;
int NUMBER_OF_ITERATIONS_BETWEEN_SAMPLES=1000; 
int BURN_IN_PERIOD=0;//50000;
double dRhoParameter=1e-1;

// NCNC proposal flags
bool gbLocationsSwitched=0; // indicates whether the SNP order was changed this iteration
bool gbLocationsMoved=0; // indicates whether a SNP locus has been moved (along the sequence)

long int ranseed=-78;
const double dEPSILON=1e-8;
int iMCscheme=1; //=2 for pure rejection
double gdRejectionMethodTolerance=1;
long int iNormalCount,iMutantCount;
double gdDist; // the distance form the target in a rejection method

// some declarations to make the MCMC calculations easier
int giOldNoRecs,giNewNoOfRecs;
double gdNewLnProductOfPriorsTerm,gdOldLnProductOfPriorsTerm;
double gdNewLogProbOfData,gdOldLogProbOfData;


#define PI 3.141592654



//************************************************
/*  Random number generator from Numerical recipes (Ed 2, ran1) */
#define IA 16807
#define IM 2147483647      
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double rand01(long *idum) 
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;
	if(*idum<=0 || !iy) {
		if(-(*idum) < 1) *idum=1;
		else *idum= -(*idum);
		for(j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ) -IR*k;
			if(*idum < 0) *idum += IM;
			if(j <NTAB) iv[j]= *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if(*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j]; 
	iv[j]= *idum;
	if((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#ifdef TRAND
main() {
	int i,bin[10];
	long dummy= -9284322;
	for(i=0;i<10;i++) bin[i]=0;
	for(i=0;i<20000;i++) 
		bin[(int)10*rand01(&dummy)]++;
	for(i=0;i<10;i++)
		printf("\n\t%d %lf",i,bin[i]/20000.0);
}
#endif

long int factorial(int n)
{
	long int dFact;
	if (n==0)
		dFact=1;
	else
		dFact=n*factorial(n-1);
	return dFact;
}

//***************further prob functions********************
double PoissonProb(double dParameter,int iNo)
{
	return pow(dParameter,iNo)*exp(-1*dParameter)/factorial(iNo);
}

//***********************************************************

void PrintArray(char *cArray)
{
	ofile<<endl;
	for (int i=0;i<iNoSNPs;i++)
		ofile<<cArray[i];
}

void PrintdArray(double *dArray)
{
	ofile<<"  ";
	for (int i=0;i<iNoSNPs;i++)
		ofile<<dArray[i]<<" ";
}

void PrintiArray(int *iArray)
{
	ofile<<"  ";
	for (int i=0;i<iNoSNPs;i++)
		ofile<<iArray[i]<<" ";
}

void SetUpHeaderLineOfOutputFile() // Writes the column headings for the output file
{
}

void AcceptNewState(int iProposalType, int iIt, int iTypeOfMCMC)
{
	// update global record-keepers
	int i,m;
	giOldNoRecs=giNewNoOfRecs;
	gdOldLnProductOfPriorsTerm=gdNewLnProductOfPriorsTerm;
	gdOldLogProbOfData=gdNewLogProbOfData;
	for (m=0;m<iNoLines;m++)
		for (i=0;i<iNoSNPs;i++)
		{
			if (iIt>BURN_IN_PERIOD)
			{
				switch(cRILarray[m][i])
				{
				case 'A':
					RILposterior[m][i][0]++;
					break;
				case 'C':
					RILposterior[m][i][1]++;
					break;
				case 'G':
					RILposterior[m][i][2]++;
					break;
				case 'T':
					RILposterior[m][i][3]++;
					break;
				}
			}
		}
}

void RejectNewState(int iProposalType,int iIterationCount,int iTypeOfMCMC)
{
	// The parameters indicates the type of proposal we made and the current iteration number
	int i,m;
	// We simply return the base to what it was before
	if (cRILarray[giLine][giBase]==cParentArray1[giBase])
	{
		cRILarray[giLine][giBase]=cParentArray2[giBase];
	}
	else
	{
		cRILarray[giLine][giBase]=cParentArray1[giBase];
	}

	//need to update one global record-keeper
	for (m=0;m<iNoLines;m++)
		for (i=0;i<iNoSNPs;i++)
		{
			if (iIterationCount>BURN_IN_PERIOD)
			{
				switch(cRILarray[m][i])
				{
				case 'A':
					RILposterior[m][i][0]++;
					break;
				case 'C':
					RILposterior[m][i][1]++;
					break;
				case 'G':
					RILposterior[m][i][2]++;
					break;
				case 'T':
					RILposterior[m][i][3]++;
					break;
				}
			}
		}

		// return SNPs to previous locations if they were moved
		if (gbLocationsSwitched)
		{
			for (i=0;i<iNoSNPs;i++)
				iSNPIndices[i]=giOldSNPIndices[i];
			gbLocationsSwitched=0;
		}

		// return SNPs to previous locations if they were moved
		if (gbLocationsMoved)
		{
			for (i=0;i<iNoSNPs;i++)
				gdSNPLocations[i]=gdOldSNPLocations[i];
			gbLocationsMoved=0;
		}
}

void SwitchSNPLocations() // switches the ordering of two loci (leaving their locations the same, but switched)
{
	// Pick two SNPs and switch them
	double p;
	int i,j,k;
	for (i=0;i<iNoSNPs;i++)
		giOldSNPIndices[i]=iSNPIndices[i];
	p=rand01(&ranseed);
	i=(int)floor(p*iNoSNPs);
	p=rand01(&ranseed);
	j=(int)floor(p*(iNoSNPs-1));
	if (j>=i)
		j++;
	k=iSNPIndices[i];
	iSNPIndices[i]=iSNPIndices[j];
	iSNPIndices[j]=k;
}

void MoveSNPLocations() // changes the location of a single locus
{
	// Pick a SNP and switch them
	double p,dLow,dHigh;
	int i;
	for (i=0;i<iNoSNPs;i++)
		gdOldSNPLocations[i]=gdSNPLocations[i];
	p=rand01(&ranseed);
	i=(int)floor(p*iNoSNPs);
	switch (i)
	{
	case 0: 
		dLow=0;
		dHigh=gdSNPLocations[1];
		break;
	case iNoSNPs-1:
		dLow=gdSNPLocations[i-1];
		dHigh=1;
		break;
	default:
		dLow=gdSNPLocations[i-1];
		dHigh=gdSNPLocations[i+1];
	}
	p=rand01(&ranseed);
	gdSNPLocations[i]=dLow+p*(dHigh-dLow);
}

int ProposeNewState(int iIterationCount,int iTypeOfMCMC) // Returns an integer indicating the type of proposal made.
{
	double p;
	if (bEstimateLocations)
	{
		// sometime we re-order the bases
		p=rand01(&ranseed);
		if (p<0.5)
		{
			SwitchSNPLocations();
			gbLocationsSwitched=1;
		}
		else
			gbLocationsSwitched=0;

		// sometimes we move some of the loci
		p=rand01(&ranseed);
		if (p<0.5)
		{
			MoveSNPLocations();
			gbLocationsMoved=1;
		}
		else
			gbLocationsMoved=0;

	}


	// pick a random base, on a random line and change its state
	p=rand01(&ranseed);
	giLine=(int)floor(p*iNoLines);
	p=rand01(&ranseed);
	giBase=(int)floor(p*iNoSNPs); // this is the base we change
	if (cRILarray[giLine][giBase]==cParentArray1[giBase])
	{
		cRILarray[giLine][giBase]=cParentArray2[giBase];
	}
	else
	{
		cRILarray[giLine][giBase]=cParentArray1[giBase];
	}	
	return 1;


}

double CalculateLnOfPriorTerm(int iTypeOfProposal,int iMCMCtype) // The parameter indicates the type of proposal made
{
	return 0;
}

double CalculateDistanceToData()
{
	return -9;
}

void GenerateSingleRejectionSample()
{
	

}


int FindIndexOfIthSNP(int Ith) // returns the index of the ith SNP, where the lefmost SNP is SNP 0.
{
	int i;
	bool bFound=0;
	//double dIthLocation;
	//for (i=0;i<iNoSNPs;i++)
	//	iCopy[i]=iSNPIndices[i];
	//ofile<<"\nUnsorted= ";
	//for (i=0;i<iNoSNPs;i++)	
	//ofile<<iCopy[i]<<" ";
	//sort(iCopy,iCopy+iNoSNPs);
	//ofile<<"\nsorted= ";
	//for (i=0;i<iNoSNPs;i++)
	//ofile<<iCopy[i]<<" ";
	//qsort((void*) dSNPLocations[0],iNoSNPs+1,sizeof (int),compare);
	//dIthLocation=iCopy[Ith];
	i=0;
	if (iSNPIndices[i]==Ith)
		bFound=1;
	while (bFound==0)
	{
		i++;
		if (iSNPIndices[i]==Ith)
			bFound=1;
		if (i>=iNoSNPs)
		{
			ofile<<"\nFailed to find ith SNP";
			//ofile<<"\ndLoc="<<dIthLocation;
			PrintiArray(iSNPIndices);
			exit(0);
		}
	}
	return i;
}

void UpdateIncidenceMx()
{
	int i,j,k;
	for (i=0;i<iNoSNPs;i++)
	{
		k=iSNPIndices[i];
		//ofile<<"\nSNP "<<i<<" is next to SNPs ";
		for (j=0;j<iNoSNPs;j++)
		{
			if ((iSNPIndices[j]==k-1)||(iSNPIndices[j]==k+1))
			{
				giIncidenceMx[i][j]++;
				giIncidenceMx[j][i]++;
				//ofile<<j<<"  ";
			}
		}
	}
}

void UpdateDistanceMx()
{
	int i,j,k,m;
	double dDist;
	//ofile<<endl;
	//PrintiArray(iSNPIndices);
	//PrintdArray(gdSNPLocations);
	for (i=0;i<iNoSNPs;i++)
	{
		k=iSNPIndices[i];
		//ofile<<"\nSNP "<<i<<" is next to SNPs ";
		for (j=0;j<iNoSNPs;j++)
		{
			if (j!=i)
			{
				m=iSNPIndices[j];
				if (k<m)
					dDist=gdSNPLocations[m]-gdSNPLocations[k];
				else
				{
					if (k>m)
						dDist=gdSNPLocations[k]-gdSNPLocations[m];
					else
					{
						ofile<<"\nError updating distance mx.";
						exit(0);
					}
				}
				//ofile<<"\n"<<dDist;
				if (dDist<=0)
				{
					ofile<<"\nError finding distance.";
					exit(0);
				}
				gdDistanceMx[i][j]+=dDist;
				gdDistanceMx[j][i]+=dDist;
			}
		}
	}
}

int CalculateNumberOfRecombinationsInRIL()
{
	// returns the minimum number of recombinations that must have occured in the RILs
	bool bParentSet=0,bMatch1,bMatch2;
	int iCurrentParent=-9;
	int i,iNoRecs=0,j,m;
	for (m=0;m<iNoLines;m++)
	{
		bParentSet=0;
		//ofile<<endl<<iNoRecs<<" ";
		for (i=0;i<iNoSNPs;i++)
		{
			//ofile<<cRILarray[m][i];
			// does the ith element match parent 1 and 2?
			switch (bEstimateLocations)
			{
			case 0: // SNPs are listed 'in order'	
				bMatch1=(cParentArray1[i]==cRILarray[m][i]);
				bMatch2=(cParentArray2[i]==cRILarray[m][i]);
				break;
			case 1: // SNPs may not be in order
				j=FindIndexOfIthSNP(i);
				//ofile<<"\nj="<<j;
				bMatch1=(cParentArray1[j]==cRILarray[m][j]);
				bMatch2=(cParentArray2[j]==cRILarray[m][j]);
				break;
			}
			//ofile<<endl<<cRILarray[i]<<cParentArray1[i]<<cParentArray2[i];
			if (bMatch1==bMatch2)
			{
				// no need to change the parent
			}
			else
			{
				if (bMatch1==1)
				{
					if (bParentSet==0)
					{
						iCurrentParent=1;
						bParentSet=1;
					}
					else
					{
						if (iCurrentParent==2)
						{
							iNoRecs++;
							iCurrentParent=1;
						}
					}
				}
				else
				{
					if (bParentSet==0)
					{
						iCurrentParent=2;
						bParentSet=1;
					}
					else
					{
						if (iCurrentParent==1)
						{
							iNoRecs++;
							iCurrentParent=2;
						}
					}
				}
				if ((bMatch1==0)&&(bMatch2==0))
				{
					ofile<<"\nTrying to calculate number of recs when the RIL doesn't match either parent. Exit.";
					exit(0);
				}
			}

		}
	}
	return iNoRecs;
}

double 	CalculateLnProductOfPriorTermsForRIL()
{
	double dProb=0;
	for (int m=0;m<iNoLines;m++)
	for (int i=0;i<iNoSNPs;i++)
	{
		switch (cRILarray[m][i])
		{
		case 'A':
			dProb+=log(RILprior[m][i][0]);
			//ofile<<"\nA p="<<dProb;
			break;
		case 'C':
			dProb+=log(RILprior[m][i][1]);
			//ofile<<"\nC p="<<dProb;
			break;
		case 'G':
			dProb+=log(RILprior[m][i][2]);
			//ofile<<"\nG p="<<dProb;
			break;
		case 'T':
			dProb+=log(RILprior[m][i][3]);
			//ofile<<"\nT p="<<dProb;
			break;
		default:
			ofile<<"\nDefault case in CalculateLnProductOfPriorTermsForRIL(). Exit.";exit(0);
		}
	}

	return dProb;
}


double CalculateLnProbOfDataTerm(int iIterationNumber, bool bLastProposalAccepted,int iMCscheme) // the parameters indicate the iteration number and  whether the last proposal was accepted 
{

	double dP1,dP2=1,dP3;
	int iRec,iExtraRecs;
	switch (iMCscheme)
	{
	case 1: // MCMC
		iRec=CalculateNumberOfRecombinationsInRIL();
		giNewNoOfRecs=iRec;
		iExtraRecs=iRec-giOldNoRecs;
		if (iExtraRecs==0)
			return 0;
		dP1=pow(dRhoParameter,iExtraRecs);
		if (iExtraRecs>0)
		{
			for (int i=giOldNoRecs+1;i<giOldNoRecs+iExtraRecs+1;i++)
				dP2=dP2*i;
			dP1=log(dP1/dP2);
		}
		else
		{
			for (int i=giOldNoRecs;i>giOldNoRecs+iExtraRecs;i--)
				dP2=dP2*i;
			dP1=log(dP2/dP1);
		}
		//ofile<<"\nProb terms: "<<iExtraRecs<<" "<<giOldNoRecs<<" "<<dP1<<" "<<dP2<<" "<<gdNewLnProductOfPriorsTerm<<" "<<gdOldLnProductOfPriorsTerm;
		gdNewLnProductOfPriorsTerm=CalculateLnProductOfPriorTermsForRIL();
		gdNewLogProbOfData=dP1+gdNewLnProductOfPriorsTerm;
		//ofile<<"\nProb terms: "<<iRec<<" "<<dP1<<" "<<dP2<<" "<<gdNewLnProductOfPriorsTerm<<" "<<gdOldLnProductOfPriorsTerm;
		return dP1+gdNewLnProductOfPriorsTerm-gdOldLnProductOfPriorsTerm;
		break;
	case 2: // rejection
		ofile<<"\nYet to write CalculateLnProbOfDataTerm() case for rejection method. Exit.";exit(0);
		break;
	case 3: // no-likelihoods
		ofile<<"\nYet to write CalculateLnProbOfDataTerm() case for 'no-l/hood' method. Exit.";exit(0);
		break;

	}
	return 0;
}

double CalculateLnTransitionKernelTerm(int iTypeOfProposal, int iMCMCtype)
{
	return 0;
}

double CalculateLnJacobianTerm(int iTypeOfProposal, int iMCMCtype) // the parameter indiactes the type of proposal made
{
	return 0;
}

void OutputStatisticsOfInterest(int iIterationCount, double dLnHastingsRatio,bool bAccept, double dMeanAcceptanceRate, int iProposalType) // The function for outputting the state of the process during the MCMC run.
{
	PrintArray(cRILarray[0]);
	ofile<<"  #Recs="<<CalculateNumberOfRecombinationsInRIL()<<" "<<giBase<<" "
		<<dLnHastingsRatio<<" "<<bAccept<<" "<<gdOldLogProbOfData;
	if (bEstimateLocations)
	{
		PrintiArray(iSNPIndices);
		PrintdArray(gdSNPLocations);
	}
}

void DeclareMemory()
{
	int i,j,m,k;
	cRILarray= new char*[iNoLines];
	for (i=0;i<iNoLines;i++)
		cRILarray[i]= new char[iNoSNPs];
	cParentArray1= new char[iNoSNPs];
	cParentArray2= new char[iNoSNPs];
	iSNPIndices= new int[iNoSNPs];
	giOldSNPIndices= new int[iNoSNPs];
	//iCopy= new int[iNoSNPs];
	RILprior= new double**[iNoLines];
	for (i=0;i<iNoLines;i++)
	{
		RILprior[i]= new double*[iNoSNPs];
		for (j=0;j<iNoSNPs;j++)
			RILprior[i][j]= new double[iNoOfAlleles];
	}
	RILposterior= new double**[iNoLines];
	for (i=0;i<iNoLines;i++)
		RILposterior[i]= new double*[iNoSNPs];
	for (j=0;j<iNoLines;j++)
	{
		for (i=0;i<iNoSNPs;i++)
		{
			RILposterior[j][i]= new double[iNoOfAlleles];
			for (k=0;k<iNoOfAlleles;k++)
				RILposterior[j][i][k]=0;
		}
	}

	if (bEstimateLocations)
	{
		gdSNPLocations=new double [iNoSNPs];
		gdOldSNPLocations=new double [iNoSNPs];
		giIncidenceMx= new int*[iNoSNPs];
		for (i=0;i<iNoSNPs;i++)
			giIncidenceMx[i]= new int[iNoSNPs];
		gdDistanceMx= new double*[iNoSNPs];
		for (i=0;i<iNoSNPs;i++)
			gdDistanceMx[i]= new double[iNoSNPs];
		for (i=0;i<iNoSNPs;i++)
			for (j=0;j<iNoSNPs;j++)
			{
				giIncidenceMx[i][j]=0;
				gdDistanceMx[i][j]=0;
			}

		
	}


	
	for (i=0;i<iNoSNPs;i++)
	{
		for (m=0;m<iNoLines;m++)
			cRILarray[m][i]='x';
		cParentArray1[i]='A';
		cParentArray2[i]='T';
	}
}

void DeleteMemory()
{
	delete cRILarray;
	delete cParentArray1;
	delete cParentArray2;
	delete iSNPIndices;
	delete giOldSNPIndices;
	delete gdSNPLocations;
	delete gdOldSNPLocations;
	delete gdDistanceMx;
	delete giIncidenceMx;

	//delete iCopy;
	delete RILprior;
	delete RILposterior;
}

void RandomizeLocations()
{
	double p;
	int i,j,k;
	for (i=iNoSNPs-1;i>0;i--)
	{
		p=rand01(&ranseed);
		j=(int)floor(p*i);
		k=iSNPIndices[j];
		iSNPIndices[j]=iSNPIndices[i];
		iSNPIndices[i]=k;
	}
	ofile<<"\nAfter randomizing:";
	for (i=0;i<iNoSNPs;i++)
		ofile<<iSNPIndices[i]<<" ";

}

void ReadData()
{
	int iTemp;
	double dJohnsPriorCertainty=0.90;
	switch(bJohnsFormat)
	{
	case 0:
		for (int m=0;m<iNoLines;m++)
		{
			ofile<<"\nPrior for line "<<m<<endl;
			for (int i=0;i<iNoSNPs;i++)
			{
				infile1>>RILprior[m][i][0];
				infile1>>RILprior[m][i][1];
				infile1>>RILprior[m][i][2];
				infile1>>RILprior[m][i][3];
				ofile<<endl<<RILprior[m][i][0]<<" "<<RILprior[m][i][1]<<" "<<RILprior[m][i][2]<<" "<<RILprior[m][i][3];
			}
		}
		for (int i=0;i<iNoSNPs;i++)
			iSNPIndices[i]=i;
		if (bEstimateLocations)
		{
			for (int i=0;i<iNoSNPs;i++)
				infile1>>gdSNPLocations[i];
		}
		break;
	case 1:
		for (int m=0;m<iNoLines;m++)
		{
			for (int i=0;i<iNoSNPs;i++)
			{
				infile1>>iTemp;
				switch (iTemp)
				{
				case 0:
					RILprior[m][i][0]=dJohnsPriorCertainty;
					RILprior[m][i][1]=(1-dJohnsPriorCertainty)/3.0;
					RILprior[m][i][2]=(1-dJohnsPriorCertainty)/3.0;
					RILprior[m][i][3]=(1-dJohnsPriorCertainty)/3.0;
					break;
				case 2:
					RILprior[m][i][0]=(1-dJohnsPriorCertainty)/3.0;
					RILprior[m][i][1]=(1-dJohnsPriorCertainty)/3.0;
					RILprior[m][i][2]=(1-dJohnsPriorCertainty)/3.0;
					RILprior[m][i][3]=dJohnsPriorCertainty;
					break;
				default:
					ofile<<"\nRead something other than 0 or 2 in ReadData(). Exit.";
					exit(0);
				}
			}
		}

		for (int i=0;i<iNoSNPs;i++)
			iSNPIndices[i]=i;
		if (bEstimateLocations)
		{
			for (int i=0;i<iNoSNPs;i++)
			{
				if (i<5)
					gdSNPLocations[i]=0.09*i;
				else
					gdSNPLocations[i]=0.09*i+0.09;
			}
		}
		break;
	}

	// store the most likely base call in the input data
	/*for (int m=0;m<iNoLines;m++)
	{
		for (int i=0;i<iNoSNPs;i++)
		{
			j=FindBiggest(RILprior[m][i][0],RILprior[m][i][1],RILprior[m][i][2],RILprior[m][i][3]);
			InitialBestCall[m][i]=j;
		}
	}*/



	//ofile<<"\nSNP locations: ";
	//PrintdArray(gdSNPLocations);
	//for (int i=0;i<iNoSNPs;i++)
	//	ofile<<iSNPIndices[i]<<" ";
	//	if (bEstimateLocations)
	//		RandomizeLocations();

}

void OutputPrior(double **dPrior)
{
}

void PrintMx(double **dMx, int iDim1, int iDim2)
{
	int i,j;
	for (i=0;i<iDim1;i++)
	{
		ofile<<endl;
		for (j=0;j<iDim2;j++)
			ofile<<dMx[i][j]<<" ";
	}
}

void SetState(char **cMx, int iElement, int iState, int iLine)
{
	// sets the iElement element of mx cMx equal to the base corresponding to iState
	switch (iState)
	{
	case 0:
		cMx[iLine][iElement]='A';
		break;
	case 1:
		cMx[iLine][iElement]='C';
		break;
	case 2:
		cMx[iLine][iElement]='G';
		break;
	case 3:
		cMx[iLine][iElement]='T';
		break;
	default:
		ofile<<"\nError. Default case in SetState(). Exit.";exit(0);
	}

}

void InitializeRILstate()
{
	// set each base to it's most likely state
	double dPMax;
	int i,j,k,m;
	for (m=0;m<iNoLines;m++)
		for (i=0;i<iNoSNPs;i++)
		{
			dPMax=RILprior[m][i][0];
			k=0;
			for (j=1;j<iNoOfAlleles;j++)
				if (dPMax<RILprior[m][i][j])
				{
					dPMax=RILprior[m][i][j];
					k=j;
				}
			SetState(cRILarray,i,k,m);
		}
		//ofile<<"\nInitial state:";
		//PrintArray(cRILarray[0]);
}

void OutputFlags()
{
	ofile<<"\nEstimating RILs.";
	if (bEstimateLocations==1)
		ofile<<"\nEstimating Locations.";
	else
		ofile<<"\nLocations Are Fixed";
}

void SortLocations()
{
}

int FindBiggest(double d0, double d1, double d2, double d3)
{
	int iBig=0;
	if (d1>d0)
	{
		d0=d1;
		iBig=1;
	}
	if (d2>d0)
	{
		d0=d2;
		iBig=2;
	}
	if (d3>d0)
	{
		d0=d3;
		iBig=3;
	}
	return iBig;
}

bool bDoPriorAndPosteriorCallsAgree(int iLine, int iSNP)
{
	int j,k;
	j=FindBiggest(RILprior[iLine][iSNP][0],RILprior[iLine][iSNP][1],RILprior[iLine][iSNP][2],RILprior[iLine][iSNP][3]);
	k=FindBiggest(RILposterior[iLine][iSNP][0],RILposterior[iLine][iSNP][1],RILposterior[iLine][iSNP][2],RILposterior[iLine][iSNP][3]);
	if (j==k)
		return 1;
	else 
		return 0;
}



main()
{

	ofile.open("output1.txt");
	ofile.precision(5);
	DeclareMemory();
	OutputFlags();
	ReadData();
	PrintArray(cParentArray1);
	PrintArray(cParentArray2);
	InitializeRILstate();
	//SortLocations();
	ofile<<"\nNoRecs="<<CalculateNumberOfRecombinationsInRIL();

	DoMCMC(NUMBER_OF_MCMC_ITERATIONS,iMCscheme);

	DeleteMemory();
}
