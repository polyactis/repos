
#include "FusionCHPData.h"
#include "FusionCHPLegacyData.h"
#include <iostream>

using namespace std;
using namespace affymetrix_fusion_io;

int main(int argc, char **argv)
{
	const char* fileName = argv[1];
	try
	{
		FusionCHPData *chp = FusionCHPDataReg::Read(fileName);
		if (chp == NULL)
		{
			cout << "Error reading the file." << endl;
			return 0;
		}
		FusionCHPLegacyData *legchp = FusionCHPLegacyData::FromBase(chp);
		if (legchp == NULL)
		{
			cout << "The example is for expression CHP files only" << endl;
			return 0;
		}
		if (legchp->GetHeader().GetAssayType() != FusionExpression)
		{
			cout << "The example is for expression only." << endl;
			return 0;
		}

		FusionExpressionProbeSetResults psResults;
		float sum = 0;
		int n = legchp->GetHeader().GetNumProbeSets();
		for (int i = 0; i < n; i++)
		{
			legchp->GetExpressionResults(i, psResults);
			sum += psResults.GetSignal();
		}
		float avg = sum / n;
		cout << "The average signal is: " << avg << endl;

		// Delete the chp object.
		delete legchp;
	}
	catch (...)
	{
		cout << "Error reading the file." << endl;
	}
	return 0;
}
