////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

#ifndef _CalvinCHPMultiDataFileBufferWriterTest_HEADER_
#define _CalvinCHPMultiDataFileBufferWriterTest_HEADER_

#include <cppunit/extensions/HelperMacros.h>

class CalvinCHPMultiDataFileBufferWriterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( CalvinCHPMultiDataFileBufferWriterTest);

	CPPUNIT_TEST ( testBuffer1 );
	CPPUNIT_TEST ( testBuffer2 );
	CPPUNIT_TEST ( testBuffer3 );

	CPPUNIT_TEST_SUITE_END();

    void UpdateFile1();
    void UpdateFile2(std::vector<std::string> &fileNames, std::vector<int> offset);
    void UpdateFile3(std::vector<std::string> &fileNames, std::vector<int> offset);
	void CreateReferenceFile1();
	void CreateReferenceFile2(const char *fileName, int offset);
	void CreateReferenceFile3(const char *fileName, int offset);
    void ChecFile2(const char *fileName, int offset);
    void ChecFile3(const char *fileName, int offset);

public:
	void setUp();
	void tearDown();

	void testBuffer1();
	void testBuffer2();
	void testBuffer3();
};

#endif
