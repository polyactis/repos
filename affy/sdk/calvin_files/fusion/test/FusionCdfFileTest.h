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
#pragma once

#include "CDFFileData.h"
#include "FusionCDFData.h"

#include <cppunit/extensions/HelperMacros.h>

class FusionCdfFileTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( FusionCdfFileTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( testExpressionV3 );
	CPPUNIT_TEST ( testExpressionXDA );
	CPPUNIT_TEST ( testMissingFile );
	CPPUNIT_TEST ( testMappingV3 );
	CPPUNIT_TEST ( testMappingXDA );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();
	void testRegObjects(affxcdf::CCDFFileData &gcosCdf, affymetrix_fusion_io::FusionCDFData &fusionCdf);
	void testQCObjects(affxcdf::CCDFFileData &gcosCdf, affymetrix_fusion_io::FusionCDFData &fusionCdf);

	void testExpressionV3();
	void testExpressionXDA();
	void testMissingFile();
	void testMappingV3();
	void testMappingXDA();

	void testCreation();
};
