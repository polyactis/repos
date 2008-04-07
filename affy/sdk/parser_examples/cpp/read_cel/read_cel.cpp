
#include "fusion/src/FusionCELData.h"
#include "fusion/src/FusionCDFData.h"
/*#include "file/CDFFileData.h"*/
#include <iostream>

using namespace std;
using namespace affymetrix_fusion_io;

void CheckCDF(FusionCDFData &cdf)
{
	FusionCDFFileHeader cdfheader;
	cdfheader = cdf.GetHeader();

	cout << "For " << cdf.GetFileName() << endl;
	cout << "\t ChipType = " << cdf.GetChipType() << endl;
	std::vector<std::string> chipTypes = cdf.GetChipTypes();
	cout << "\t ChipTypes = ";
	for (int i=0; i<chipTypes.size(); i++)
		cout << chipTypes[i] << " ";
	cout << endl;
	cout << "\t" << cdfheader.GetNumProbeSets() << " probesets;" << endl;
	cout << '\t' << cdfheader.GetNumQCProbeSets() << " QC probesets;\n";
	cout << '\t' << cdfheader.GetCols() << " columns;" << endl;
	cout << '\t' << cdfheader.GetRows() << " rows;" << endl;
	cout << endl;
	std::string name;
	FusionCDFProbeSetInformation cdfProbeSetInfo;
	FusionCDFQCProbeSetInformation qcProbeSetInfo;
	affxcdf::GeneChipProbeSetType probeSetType;
	affxcdf::DirectionType directionType;
	FusionCDFProbeGroupInformation group;
	FusionCDFProbeInformation probe;
	for (int i=0; i<cdfheader.GetNumProbeSets(); i++)
	{
		name = cdf.GetProbeSetName(i);
		probeSetType = cdf.GetProbeSetType(i);
		cdf.GetProbeSetInformation(i, cdfProbeSetInfo);
		cout << "\tProbeSet: " << name << ". type=" << probeSetType << " direction=" << cdfProbeSetInfo.GetDirection() \
		<< " NumLists=" << cdfProbeSetInfo.GetNumLists()\
		<< " NumGroups=" << cdfProbeSetInfo.GetNumGroups()\
		<< " NumCells=" << cdfProbeSetInfo.GetNumCells()\
		<< " NumCellsPerList=" << cdfProbeSetInfo.GetNumCellsPerList()\
		<< " ProbeSetNumber=" << cdfProbeSetInfo.GetProbeSetNumber() << endl;
		int ngroups = cdfProbeSetInfo.GetNumGroups();
		for (int igroup=0; igroup<ngroups; igroup++)
		{
			cdfProbeSetInfo.GetGroupInformation(igroup, group);
			int ncells = group.GetNumCells();
			cout << "\t\tGroup: " << group.GetName() << ". start=" << group.GetStart() << " stop=" << group.GetStop() \
			<< " direction=" << group.GetDirection()\
			<< " NumLists=" << group.GetNumLists()\
			<< " NumCells=" << group.GetNumCells()\
			<< " NumCellsPerList=" << group.GetNumCellsPerList() << endl;
			for (int icell=0; icell<ncells; icell++)
			{
				group.GetCell(icell, probe);
				cout << "\t\t\tProbe: " << probe.GetExpos() << ". ListIndex=" << probe.GetListIndex() << " PBase=" << probe.GetPBase() \
				<< " TBase=" << probe.GetTBase()\
				<< " X=" << probe.GetX()\
				<< " Y=" << probe.GetY() << endl;
			}
		}
	}
}


int main(int argc, char **argv)
{
	const char* celFileName = argv[1];
	const char* cdfFileName = argv[2];
	FusionCELData cel;
	FusionCDFData cdf;
	try
	{
		cel.SetFileName(celFileName);
		if (cel.Read() == false)
		{
			cout << "Failed to read the file." << endl;
			return 0;
		}
		int n = cel.GetNumCells();
		float sum = 0;
		for (int i = 0; i < n; i++)
		{
			sum += cel.GetIntensity(i);
		}
		float avg = sum / n;
		cout << "The average intensity is: " << avg << endl;

		cdf.SetFileName(cdfFileName);
		if (cdf.Read() == false)
		{
			cout << "Failed to read the CDF file." << endl;
			return 0;
		}

		int nsets = cdf.GetHeader().GetNumProbeSets();
		CheckCDF(cdf);


		std::string name;
		for (int iset=0; iset<nsets; iset++)
		{
			name = cdf.GetProbeSetName(iset);
			sum = 0;
			FusionCDFProbeSetInformation set;
			cdf.GetProbeSetInformation(iset, set);
			int ngroups = set.GetNumGroups();
			cout << "probe " << name << " has " << ngroups << " groups." << endl;
			for (int igroup=0; igroup<ngroups; igroup++)
			{
				FusionCDFProbeGroupInformation group;
				set.GetGroupInformation(igroup, group);
				int ncells = group.GetNumCells();
				cout << "probe " << name << " group " << igroup << " has " << ncells << " cells." << endl;
				for (int icell=0; icell<ncells; icell++)
				{
					FusionCDFProbeInformation probe;
					group.GetCell(icell, probe);
					sum += cel.GetIntensity(probe.GetX(), probe.GetY());
				}
			}
			avg = sum / set.GetNumCells();
			cout << "The average probe set intensity (" << name << ") is " << avg << endl;
		}
	}
	catch (...)
	{
		cout << "Error in reading the file.";
	}

	return 0;
}
