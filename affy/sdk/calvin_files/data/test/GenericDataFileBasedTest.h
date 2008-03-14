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

#ifndef __GENERICDATAFILEBASEDTEST_H_
#define __GENERICDATAFILEBASEDTEST_H_

#include <cppunit/extensions/HelperMacros.h>
#include "GenericData.h"

/*
 * These are the general file-based GenericData tests.
 * Only tests that could not be done in GenericDataTest have been added here.
 */
class GenericDataFileBasedTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE (GenericDataFileBasedTest);

	CPPUNIT_TEST (CreationTest);
	CPPUNIT_TEST (DataGroupMMTest);
	CPPUNIT_TEST (DataGroupFSreamTest);
	CPPUNIT_TEST (DataGroupFStreamLoadEntireDataSetTest);

	CPPUNIT_TEST_SUITE_END();

public:
	GenericDataFileBasedTest();
	~GenericDataFileBasedTest();

	void setUp();
	void tearDown();

	void CreationTest();
	void DataGroupMMTest();
	void DataGroupFSreamTest();
	void DataGroupFStreamLoadEntireDataSetTest();
};

#endif // __GENERICDATAFILEBASEDTEST_H_
