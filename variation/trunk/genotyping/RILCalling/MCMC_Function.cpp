//#include <math.h>

// Generic function for conducting MCMC routine. It calls various
// external functions, listed below, to fill in the details specific
// to a given application.
#include "MCMCFunction.h"

extern double rand01(long int*);
extern long int ranseed;
extern int NUMBER_OF_MCMC_ITERATIONS;
extern int BURN_IN_PERIOD;
extern int iMCscheme; // the type of MCMC scheme used:
	//			  =1 for pure MCMC
	//            =2 for rejection method
	//            =3 for no-likelihoods
extern 	int giOldNoRecs,giNewNoOfRecs;
extern double gdOldLnProductOfPriorsTerm,gdNewLnProductOfPriorsTerm;
extern double ***RILposterior,***RILprior;
extern int **giIncidenceMx;
extern double **gdDistanceMx;
extern bool bEstLocations;
extern int iNoSNPsTemp;
extern int iNoLinesTemp;

//extern int iNoSNPs;
extern void UpdateIncidenceMx();
extern void UpdateDistanceMx();
extern bool bDoPriorAndPosteriorCallsAgree(int,int);
extern int FindBiggest(double, double, double, double);



extern ofstream ofile;
extern void SetUpHeaderLineOfOutputFile(); // Writes thew column headings for the output file
extern void AcceptNewState(int,int,int);
extern void RejectNewState(int,int,int); // The parameters indicates the type of proposal we made and the current iteration number
extern int ProposeNewState(int,int); // Returns an integer indicating the type of proposal made.
extern double CalculateLnOfPriorTerm(int,int); // The parameterindicates the type of proposal made
extern double CalculateLnProbOfDataTerm(int,bool,int); // the parameters indicate the iteration number
// and  whether the last proposal was accepted 
extern double CalculateLnTransitionKernelTerm(int,int);
extern double CalculateLnJacobianTerm(int,int); // the parameter indiactes the type of proposal made
extern void OutputStatisticsOfInterest(int,double,bool,double,int); // The function for outputting the state of the process during the MCMC run.
extern void	GenerateSingleRejectionSample();
extern double CalculateDistanceToData();
extern double gdDist;
extern double gdRejectionMethodTolerance;

//Parameters
extern int NUMBER_OF_ITERATIONS_BETWEEN_SAMPLES;  //50//00//00//200//500//00//00

#define NUMBER_OF_TYPES_OF_PROPOSAL 20 // upper bound for the no. of different ways of proposing a new state

//Debugging flags
#define ACCEPT_HALF_OF_PROPOSALS 0 // =1 to accept half of the proposals - for debugging
//#define LONG_OUTPUT  
#define OUTPUT_HR_PROBS 0 // =1 to output the components of the HR


double CalculateLnHastingsRatio(int iIterationNumber, bool bLastProposalAccepted, int iTypeOfProposal, int iMCMCtype)
{
	//   iMCMCtype=1 for pure MCMC
	//            =2 for rejection method
	//            =3 for no-likelihoods

	double dLnPriorTerm,dLnProbOfDataTerm,dLnTransitionKernelTerm,dLnJacobianTerm;
	//if (iIterationNumber==BURN_IN_PERIOD)
	//	ofile<<"\nkjklj";
	dLnPriorTerm=CalculateLnOfPriorTerm(iTypeOfProposal,iMCMCtype);
	dLnProbOfDataTerm=CalculateLnProbOfDataTerm(iIterationNumber,bLastProposalAccepted,iMCMCtype);
	dLnTransitionKernelTerm=CalculateLnTransitionKernelTerm(iTypeOfProposal,iMCMCtype);
	dLnJacobianTerm=CalculateLnJacobianTerm(iTypeOfProposal,iMCMCtype);

	if (OUTPUT_HR_PROBS)
		ofile<<"\n* "<<dLnPriorTerm<<" "<<dLnProbOfDataTerm<<" "<<dLnTransitionKernelTerm<<" "<<dLnJacobianTerm;
	if (ACCEPT_HALF_OF_PROPOSALS)
		return log(0.5);
	else
	return dLnPriorTerm+dLnProbOfDataTerm+dLnTransitionKernelTerm
		+dLnJacobianTerm;
}

void *pFunctionPtr(); // might use this later


void DoMCMC(int iNumberOfIterations, int iTypeOfMCMC)
{
	// iTypeOfMCMC=1 for pure MCMC
	//            =2 for rejection method
	//            =3 for no-likelihoods
	int iIterationCount,iProposalType,i,j,m,kk;
	bool bAccept=0;
	double dLnProb,dLnHastingsRatio;
	double dMeanAcceptanceRate=1;
	double iaProposalFrequencies[NUMBER_OF_TYPES_OF_PROPOSAL];
	double iaAcceptanceFrequencies[NUMBER_OF_TYPES_OF_PROPOSAL];

	if (ACCEPT_HALF_OF_PROPOSALS)
		ofile<<"\nACCEPTING HALF OF THE PROPOSALS - FOR DEBUGGING";

	//Initialize everything
	// probs for hastings ratio calculations
	if (iTypeOfMCMC!=2)
		dLnHastingsRatio=CalculateLnHastingsRatio(0,0,1,iMCscheme);
	giOldNoRecs=giNewNoOfRecs;
	gdOldLnProductOfPriorsTerm=gdNewLnProductOfPriorsTerm;

	memset(iaProposalFrequencies,0,sizeof(iaProposalFrequencies));
	memset(iaAcceptanceFrequencies,0,sizeof(iaAcceptanceFrequencies));

	// Set up the header line of the output file
	SetUpHeaderLineOfOutputFile();
	// Do MCMC loop
	for (iIterationCount=0;iIterationCount<=iNumberOfIterations;iIterationCount++)
	{
		//		pFunctionPtr=&ProposeNewState;
		//Propose new state
		iProposalType=ProposeNewState(iIterationCount,iTypeOfMCMC);
		assert(iProposalType<NUMBER_OF_TYPES_OF_PROPOSAL);
		iaProposalFrequencies[iProposalType]++;

		//Calculate Hastings Ratio
		dLnHastingsRatio=CalculateLnHastingsRatio(iIterationCount,bAccept,iProposalType,iTypeOfMCMC);
#ifdef LONG_OUTPUT
		ofile<<"\nlogHR="<<dLnHastingsRatio;
#endif
		//Decide whether to accept and take appropriate action
		switch(iMCscheme)
		{
		case 1: // MCMC
			dLnProb=log(rand01(&ranseed));

			// Force the first iteration (and the iteration at the end of the burn in period to be accepted so that the program variables are seeded.
			if ((iMCscheme!=2)&&((iIterationCount==0)||(iIterationCount==BURN_IN_PERIOD)))
				bAccept=1;
			else
				bAccept=(dLnProb<dLnHastingsRatio);
			if (bAccept==1)
			{
				AcceptNewState(iProposalType,iIterationCount,iTypeOfMCMC);
				iaAcceptanceFrequencies[iProposalType]++;
			}
			else
				RejectNewState(iProposalType,iIterationCount,iTypeOfMCMC);
			break;
		case 2: // rejection
			GenerateSingleRejectionSample();
			gdDist=CalculateDistanceToData();
			if (gdDist<gdRejectionMethodTolerance)
			{
				bAccept=1;
				AcceptNewState(iProposalType,iIterationCount,iTypeOfMCMC);
				iaAcceptanceFrequencies[iProposalType]++;
			}
			else
			{
				bAccept=0;
				RejectNewState(iProposalType,iIterationCount,iTypeOfMCMC);
			}
			break;
		case 3: // no likelihoods
			ofile<<"\nNo code written for acceptance of iMCscheme=3. Exit.";
			exit(0);
			break;
		default:
			ofile<<"\nNo code written for acceptance of iMCscheme="<<iMCscheme<<". Exit.";
			exit(0);
		}

		//collect any statistics of interest
		dMeanAcceptanceRate=(dMeanAcceptanceRate*(double)iIterationCount+bAccept)
			/(double)(iIterationCount+1);

		if ((iIterationCount>BURN_IN_PERIOD)&&(bEstLocations))
		{
			UpdateIncidenceMx();
			UpdateDistanceMx();
		}

		//output statistics for every Nth sample
		switch (iMCscheme)
		{
		case 1: // pure MCMC
			if (iIterationCount % NUMBER_OF_ITERATIONS_BETWEEN_SAMPLES== 0)
				OutputStatisticsOfInterest(iIterationCount,dLnHastingsRatio,bAccept,dMeanAcceptanceRate,iProposalType);
			break;
		case 2: // rejection
			if (bAccept==1)
				OutputStatisticsOfInterest(iIterationCount,dLnHastingsRatio,bAccept,dMeanAcceptanceRate,iProposalType);
			break;
		case 3:// no-likelihoods
			break;
		}

	}

	//Print out summary statistics
	ofile<<"\n\n===============SUMMARY STATISTICS FOLLOW=============";	
	ofile<<"\nAcceptance statistics for each proposal type:";
	ofile<<"\nProposal type    No.of proposals      No.accepted     Acceptance rate";
	for (i=0;i<NUMBER_OF_TYPES_OF_PROPOSAL;i++)
		if (iaProposalFrequencies[i]>0)
			ofile<<"\n    "<<i<<"                 "<<iaProposalFrequencies[i]<<"                    "
			<<iaAcceptanceFrequencies[i]<<"          "<<iaAcceptanceFrequencies[i]/(double)iaProposalFrequencies[i];
	ofile<<"\n=====================================================";
	ofile<<"\nPrior distribution for RILs:";
	for (m=0;m<iNoLinesTemp;m++)
	{
		ofile<<"\nLine "<<m<<",  P(A): ";
		for (i=0;i<iNoSNPsTemp;i++)
			ofile<<RILprior[m][i][0]<<" ";
		ofile<<"\nP(C): ";
		for (i=0;i<iNoSNPsTemp;i++)
			ofile<<RILprior[m][i][1]<<" ";
		ofile<<"\nP(G): ";
		for (i=0;i<iNoSNPsTemp;i++)
			ofile<<RILprior[m][i][2]<<" ";
		ofile<<"\nP(T): ";
		for (i=0;i<iNoSNPsTemp;i++)
			ofile<<RILprior[m][i][3]<<" ";
	}
	ofile<<"\n=====================================================";
	ofile<<"\n=====================================================";
	ofile<<"\nFinal posterior distribution for RILs:";
	for (m=0;m<iNoLinesTemp;m++)
	{
		ofile<<"\nLine "<<m<<",  P(A): ";
		for (i=0;i<iNoSNPsTemp;i++)
			ofile<<RILposterior[m][i][0]/(iNumberOfIterations-BURN_IN_PERIOD)<<" ";
		ofile<<"\nP(C): ";
		for (i=0;i<iNoSNPsTemp;i++)
			ofile<<RILposterior[m][i][1]/(iNumberOfIterations-BURN_IN_PERIOD)<<" ";
		ofile<<"\nP(G): ";
		for (i=0;i<iNoSNPsTemp;i++)
			ofile<<RILposterior[m][i][2]/(iNumberOfIterations-BURN_IN_PERIOD)<<" ";
		ofile<<"\nP(T): ";
		for (i=0;i<iNoSNPsTemp;i++)
			ofile<<RILposterior[m][i][3]/(iNumberOfIterations-BURN_IN_PERIOD)<<" ";
	}
	ofile<<"\n=====================================================";
	if (bEstLocations)
	{
		ofile<<"\nFinal SNP incidence matrix:";
		for (i=0;i<iNoSNPsTemp;i++)
		{
			ofile<<endl;
			for (j=0;j<iNoSNPsTemp;j++)
				ofile<<0.5*(double)giIncidenceMx[i][j]/(NUMBER_OF_MCMC_ITERATIONS-BURN_IN_PERIOD)<<" ";
		}
	}
	ofile<<"\n=====================================================";
	if (bEstLocations)
	{
		ofile<<"\nFinal SNP distance matrix:";
		for (i=0;i<iNoSNPsTemp;i++)
		{
			ofile<<endl;
			for (j=0;j<iNoSNPsTemp;j++)
				ofile<<0.5*(double)gdDistanceMx[i][j]/(NUMBER_OF_MCMC_ITERATIONS-BURN_IN_PERIOD)<<" ";
		}
	}

	// check whether prior and posterior calls agree
	ofile<<"\n========================================";
	ofile<<"\n====Changed calls are listed below======";
	ofile<<"\n========================================";
	for (m=0;m<iNoLinesTemp;m++)
	{
		for (i=0;i<iNoSNPsTemp;i++)
		{
			if (bDoPriorAndPosteriorCallsAgree(m,i)==0)
			{
				ofile<<"\nChanged call on line "<<m<<" SNP "<<i<<" from ";
				kk=FindBiggest(RILprior[m][i][0],RILprior[m][i][1],RILprior[m][i][2], RILprior[m][i][3]);
				ofile<<kk<<" to ";
				kk=FindBiggest(RILposterior[m][i][0],RILposterior[m][i][1],RILposterior[m][i][2], RILposterior[m][i][3]);
				ofile<<kk;
			}
		}
	}
}
