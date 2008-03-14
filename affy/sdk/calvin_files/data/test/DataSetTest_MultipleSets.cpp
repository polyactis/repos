////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////

//
#include <string>
//
#include "DataException.h"
#include "DataSetTest_MultipleSets.h"
#include "GenericData.h"
#include "GenericDataTypes.h"
#include "GenericFileReader.h"

using namespace std;
using namespace affymetrix_calvin_io;

#define TEST_DATA_DAT_FILE_WITH_GRIDS "../data/test.file.data_dat_with_grid"

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetTest_MultipleSets );

DataSetTest_MultipleSets::DataSetTest_MultipleSets()
{
	data = 0;
	dataPlane = 0;
}

DataSetTest_MultipleSets::~DataSetTest_MultipleSets()
{

}

void DataSetTest_MultipleSets::setUp()
{
	// Open one DataSet.
	data = new GenericData;
	GenericFileReader reader;
	std::string name = TEST_DATA_DAT_FILE_WITH_GRIDS;
	reader.SetFilename(name);
	reader.ReadHeader(*data);
	try
	{
		dataPlane = data->DataSet(0,0);	// acquired data is the first DataSet of the first DataGroup.
	}
	catch(affymetrix_calvin_exceptions::CalvinException&)
	{
		dataPlane = 0;
	}
}

void DataSetTest_MultipleSets::tearDown()
{
	if (dataPlane)
		dataPlane->Delete();
	delete data;
}

void DataSetTest_MultipleSets::testCreation()
{
	// Check that the DataSet was created.
	CPPUNIT_ASSERT(dataPlane);
}

void DataSetTest_MultipleSets::testCreateSecondDataSet()
{
	DataSet* gridDP=NULL;
	CPPUNIT_ASSERT_NO_THROW(gridDP = data->DataSet(L"",L"grid position"));
	CPPUNIT_ASSERT_NO_THROW(gridDP->Delete());
}


void DataSetTest_MultipleSets::testOpenTwoDataSets()
{
	DataSet* gridDP=NULL;
	CPPUNIT_ASSERT_NO_THROW(gridDP = data->DataSet(L"",L"grid position"));

	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT(gridDP->Open() == true);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
	CPPUNIT_ASSERT_NO_THROW(gridDP->Delete());	// also calls Close
}

void DataSetTest_MultipleSets::testAccessDataFromTwoDifferentDataSets()
{
	// Create second DataSet
	DataSet* gridDP=NULL;
	CPPUNIT_ASSERT_NO_THROW(gridDP = data->DataSet(L"",L"grid position"));

	// Open both DataSets
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT(gridDP->Open() == true);

	// Access the data from both DataSets

	// Access from the grid DataSet
	int32_t rows = gridDP->Rows();
	int32_t cols = gridDP->Cols();
	for (int32_t row = 0; row < rows; ++row)
	{
		for (int32_t col = 0; col < cols; ++col)
		{
			// Compute the expected value
			float expected = (float)(row*100 + col/2); // corner = col/2

			float value;
			CPPUNIT_ASSERT_NO_THROW(gridDP->GetData(row, col, value));
			CPPUNIT_ASSERT(value == expected);
		}
	}

	// Access the acquired data DataSet
	rows = dataPlane->Rows();
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		u_int16_t value;
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 0, value));
		CPPUNIT_ASSERT(value == expected);
	}

	// Access the grid DataSet again.
	rows = gridDP->Rows();
	cols = gridDP->Cols();
	for (int32_t row = 0; row < rows; ++row)
	{
		for (int32_t col = 0; col < cols; ++col)
		{
			// Compute the expected value
			float expected = (float)(row*100 + col/2); // corner = col/2

			float value;
			CPPUNIT_ASSERT_NO_THROW(gridDP->GetData(row, col, value));
			CPPUNIT_ASSERT(value == expected);
		}
	}

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
	CPPUNIT_ASSERT_NO_THROW(gridDP->Delete());	// also calls Close
}

void DataSetTest_MultipleSets::testAccessDataFromTwoInstancesOfSamePlane()
{
	// Create second instance of the acquired data DataSet
	DataSet* gridDP=NULL;
	CPPUNIT_ASSERT_NO_THROW(gridDP = data->DataSet(L"",L"acquired data"));

	// Open both DataSets
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT(gridDP->Open() == true);

	// Access the first instance DataSet
	int32_t rows = dataPlane->Rows();
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		u_int16_t value;
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 0, value));
		CPPUNIT_ASSERT(value == expected);
	}

	// Access the second instance of the acquired data DataSet
	rows = gridDP->Rows();
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		u_int16_t value;
		CPPUNIT_ASSERT_NO_THROW(gridDP->GetData(row, 0, value));
		CPPUNIT_ASSERT(value == expected);
	}

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
	CPPUNIT_ASSERT_NO_THROW(gridDP->Delete());	// also calls Close
}
