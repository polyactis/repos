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

#ifndef __GENERICFILEUPDATERTEST_H_
#define __GENERICFILEUPDATERTEST_H_

#include "GenericFileUpdater.h"
#include <cppunit/extensions/HelperMacros.h>

using namespace affymetrix_calvin_io;

class GenericFileUpdaterTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( GenericFileUpdaterTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( UpdaterTest );
	CPPUNIT_TEST ( Updater2Test );
	CPPUNIT_TEST ( UpdateExistingDataSet1Test );	// TODO: make this test pass
	CPPUNIT_TEST ( UpdateExistingDataSet2Test );

	CPPUNIT_TEST_SUITE_END();

public:

	void setUp();
	void tearDown();
	void testCreation();
	void UpdaterTest();
	void Updater2Test();
	void UpdateExistingDataSet1Test();
	void UpdateExistingDataSet2Test();

private:

	bool Copy(std::string filename, std::string copyFilename);
	
};

#endif // __GENERICFILEUPDATERTEST_H_
