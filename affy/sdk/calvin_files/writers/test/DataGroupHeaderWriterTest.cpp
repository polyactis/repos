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

#include "DataGroupHeaderWriterTest.h"
#include "DataGroupHeaderWriter.h"
#include "ColumnInfo.h"

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( DataGroupHeaderWriterTest );

void DataGroupHeaderWriterTest::setUp()
{
    writer = new DataGroupHeaderWriter();
}

void DataGroupHeaderWriterTest::tearDown()
{
    delete writer;
}

void DataGroupHeaderWriterTest::testCreation()
{
	DataGroupHeaderWriter w;
	CPPUNIT_ASSERT(1);
}

void DataGroupHeaderWriterTest::WriteTest()
{
	std::ofstream os;
    std::string f = "data_dataGroup_header";
    os.open(f.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
	if (!os)
	{
		CPPUNIT_ASSERT(0);
	}
    else
    {
	    DataGroupHeader hdr(L"dataGroup");
	    writer->Write(os, hdr);
        CPPUNIT_ASSERT(1);
        os.close();
    }

    //TODO: use data dataGroup header reader to read header and verify data
}
