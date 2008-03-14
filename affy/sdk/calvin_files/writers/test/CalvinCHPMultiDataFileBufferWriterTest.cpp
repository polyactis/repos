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

//
#include "CHPMultiDataFileReader.h"
#include "CalvinCHPMultiDataFileBufferWriter.h"
#include "CalvinCHPMultiDataFileBufferWriterTest.h"
#include "CalvinCHPMultiDataFileWriter.h"
#include "StringUtils.h"

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_data;

CPPUNIT_TEST_SUITE_REGISTRATION( CalvinCHPMultiDataFileBufferWriterTest );

#pragma warning(disable: 4996) // ignore deprecated functions warning

#define TEST_FILE1 "multi_data_file1"
#define TEST_FILE2a "multi_data_file2a"
#define TEST_FILE2b "multi_data_file2b"
#define TEST_FILE3a "multi_data_file3a"
#define TEST_FILE3b "multi_data_file3b"

static string IntToString(int i)
{
   char buf[64];
   sprintf(buf, "%d", i);
   return buf;
}

void CalvinCHPMultiDataFileBufferWriterTest::setUp()
{
}

void CalvinCHPMultiDataFileBufferWriterTest::tearDown()
{
}

void CalvinCHPMultiDataFileBufferWriterTest::CreateReferenceFile1()
{
    int ng = 10000;
    int nx = 5000;
	CHPMultiDataData data(TEST_FILE1);
	data.SetEntryCount(GenotypeMultiDataType, ng, 10);
	data.SetEntryCount(ExpressionMultiDataType, nx, 10);
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	ProbeSetMultiDataGenotypeData gn;
	ProbeSetMultiDataExpressionData ex;

    writer->SeekToDataSet(GenotypeMultiDataType);
    gn.call = 0;
    gn.confidence = 0.0f;
    for (int i=0; i<ng; i++)
    {
        gn.name = IntToString(i);
	    writer->WriteEntry(gn);
    }
	writer->SeekToDataSet(ExpressionMultiDataType);
    ex.quantification = 0;
    for (int i=0; i<nx; i++)
    {
        ex.name = IntToString(i);
	    writer->WriteEntry(ex);
    }

	delete writer;
}

void CalvinCHPMultiDataFileBufferWriterTest::testBuffer1()
{
	CreateReferenceFile1();
    UpdateFile1();

    ProbeSetMultiDataExpressionData ex;
	ProbeSetMultiDataGenotypeData gn;
	CHPMultiDataData data;
	CHPMultiDataFileReader reader;
	reader.SetFilename(TEST_FILE1);
	reader.Read(data);

    int ng = 10000;
    int nx = 5000;

    CPPUNIT_ASSERT(data.GetEntryCount(ExpressionMultiDataType) == nx);
	CPPUNIT_ASSERT(data.GetEntryCount(GenotypeMultiDataType) == ng);

    for (int i=0; i<nx; i++)
    {
        data.GetExpressionEntry(ExpressionMultiDataType, i, ex);
    	CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.quantification, i, 0.0001f);
	    CPPUNIT_ASSERT(ex.name == IntToString(i));
    }
    for (int i=0; i<ng; i++)
    {
        data.GetGenotypeEntry(GenotypeMultiDataType, i, gn);
    	CPPUNIT_ASSERT(gn.call == i%4);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.confidence, i, 0.0001f);
	    CPPUNIT_ASSERT(gn.name == IntToString(i));
    }
}

void CalvinCHPMultiDataFileBufferWriterTest::UpdateFile1()
{
    CHPMultiDataFileBufferWriter writer;
    vector<string> chps;
    vector<MultiDataType> dataTypes;
    map<MultiDataType, int> maxNameLengths;
    dataTypes.push_back(GenotypeMultiDataType);
    dataTypes.push_back(ExpressionMultiDataType);
    chps.push_back(TEST_FILE1);
    maxNameLengths[GenotypeMultiDataType] = 10;
    maxNameLengths[ExpressionMultiDataType] = 10;
    writer.Initialize(&chps, dataTypes, maxNameLengths);
    writer.SetMaxBufferSize(102400);
    ProbeSetMultiDataExpressionData ex;
	ProbeSetMultiDataGenotypeData gn;

    int ng = 10000;
    int nx = 5000;

    for (int i=0; i<nx; i++)
    {
        ex.quantification = (float)i;
        ex.name = IntToString(i);
        writer.WriteMultiDataExpressionEntry(ExpressionMultiDataType, 0, ex);
    }
    for (int i=0; i<ng; i++)
    {
        gn.name = IntToString(i);
        gn.call = i%4;
        gn.confidence = (float) i;
        writer.WriteMultiDataGenotypeEntry(GenotypeMultiDataType, 0, gn);
    }
    writer.FlushBuffer();

}

void CalvinCHPMultiDataFileBufferWriterTest::CreateReferenceFile2(const char *fileName, int offset)
{
    vector<ColumnInfo> cols;
    ParameterNameValueType nv;

    IntColumn icol(L"int");
    cols.push_back(icol);

    FloatColumn fcol(L"float");
    cols.push_back(fcol);

    int ng = 10000;
    int nx = 5000;
	CHPMultiDataData data(fileName);
	data.SetEntryCount(ExpressionMultiDataType, nx, 10, cols);
	data.SetEntryCount(GenotypeMultiDataType, ng, 10, cols);
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	ProbeSetMultiDataGenotypeData gn;
	ProbeSetMultiDataExpressionData ex;

    nv.SetName(L"int");
    nv.SetValueInt32(0);
    ex.metrics.push_back(nv);
    gn.metrics.push_back(nv);
    nv.SetName(L"float");
    nv.SetValueFloat(0);
    ex.metrics.push_back(nv);
    gn.metrics.push_back(nv);

    writer->SeekToDataSet(GenotypeMultiDataType);
    gn.call = 0;
    gn.confidence = 0.0f;
    for (int i=0; i<ng; i++)
    {
        gn.name = IntToString(i+offset);
	    writer->WriteEntry(gn);
    }

	writer->SeekToDataSet(ExpressionMultiDataType);
    ex.quantification = 0;
    for (int i=0; i<nx; i++)
    {
        ex.name = IntToString(i+offset);
        writer->WriteEntry(ex);
    }


	delete writer;

}

void CalvinCHPMultiDataFileBufferWriterTest::testBuffer2()
{
	CreateReferenceFile2(TEST_FILE2a, 0);
	CreateReferenceFile2(TEST_FILE2b, 1001);
    vector<string> chps;
    vector<int> offset;
    chps.push_back(TEST_FILE2a);
    offset.push_back(0);
    chps.push_back(TEST_FILE2b);
    offset.push_back(1001);
    UpdateFile2(chps, offset);
    ChecFile2(TEST_FILE2a, 0);
    ChecFile2(TEST_FILE2b, 1001);
}

void CalvinCHPMultiDataFileBufferWriterTest::ChecFile2(const char *fileName, int offset)
{
    ProbeSetMultiDataExpressionData ex;
	ProbeSetMultiDataGenotypeData gn;
	CHPMultiDataData data;
	CHPMultiDataFileReader reader;
	reader.SetFilename(fileName);
	reader.Read(data);

    int ng = 10000;
    int nx = 5000;

    CPPUNIT_ASSERT(data.GetEntryCount(ExpressionMultiDataType) == nx);
	CPPUNIT_ASSERT(data.GetEntryCount(GenotypeMultiDataType) == ng);

    for (int i=0; i<nx; i++)
    {
        data.GetExpressionEntry(ExpressionMultiDataType, i, ex);
    	CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.quantification, i+offset, 0.0001f);
	    CPPUNIT_ASSERT(ex.name == IntToString(i+offset));

        CPPUNIT_ASSERT(ex.metrics[0].GetValueInt32() == i+offset);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.metrics[1].GetValueFloat(), 2*i+offset, 0.00001f);
    }
    for (int i=0; i<ng; i++)
    {
        data.GetGenotypeEntry(GenotypeMultiDataType, i, gn);
    	CPPUNIT_ASSERT(gn.call == (i+offset)%4);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.confidence, i+offset, 0.0001f);
	    CPPUNIT_ASSERT(gn.name == IntToString(i+offset));

        CPPUNIT_ASSERT(gn.metrics[0].GetValueInt32() == i+offset);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.metrics[1].GetValueFloat(), 2*i+offset, 0.00001f);
    }
}

void CalvinCHPMultiDataFileBufferWriterTest::UpdateFile2(std::vector<std::string> &fileNames, vector<int> offset)
{
    ParameterNameValueType nv;
    CHPMultiDataFileBufferWriter writer;
    vector<MultiDataType> dataTypes;
    dataTypes.push_back(ExpressionMultiDataType);
    dataTypes.push_back(GenotypeMultiDataType);
    map<MultiDataType, int> maxNameLengths;
    maxNameLengths[ExpressionMultiDataType]=10;
    maxNameLengths[GenotypeMultiDataType]=10;
    writer.Initialize(&fileNames, dataTypes, maxNameLengths);
    writer.SetMaxBufferSize(102400);

    int ng = 10000;
    int nx = 5000;

    ProbeSetMultiDataExpressionData ex;
    ProbeSetMultiDataGenotypeData gn;

    nv.SetName(L"int");
    nv.SetValueInt32(0);
    ex.metrics.push_back(nv);
    gn.metrics.push_back(nv);
    nv.SetName(L"float");
    nv.SetValueFloat(0);
    ex.metrics.push_back(nv);
    gn.metrics.push_back(nv);
    for (int i=0; i<nx; i++)
    {
        for (int ifile=0; ifile<(int)fileNames.size(); ifile++)
        {
            ex.name = IntToString(i+offset[ifile]);
            ex.quantification = (float)(i+offset[ifile]);
            ex.metrics[0].SetValueInt32(i+offset[ifile]);
            ex.metrics[1].SetValueFloat(2.0f*i+offset[ifile]);
            writer.WriteMultiDataExpressionEntry(ExpressionMultiDataType, ifile, ex);
        }
    }

    for (int i=0; i<ng; i++)
    {
        for (int ifile=0; ifile<(int)fileNames.size(); ifile++)
        {
            gn.name = IntToString(i+offset[ifile]);
            gn.call = (i+offset[ifile])%4;
            gn.confidence = (float)(i+offset[ifile]);
            gn.metrics[0].SetValueInt32(i+offset[ifile]);
            gn.metrics[1].SetValueFloat(2.0f*i+offset[ifile]);
            writer.WriteMultiDataGenotypeEntry(GenotypeMultiDataType, ifile, gn);
        }
    }
    writer.FlushBuffer();
}

void CalvinCHPMultiDataFileBufferWriterTest::CreateReferenceFile3(const char *fileName, int offset)
{
    vector<ColumnInfo> cols;
    ParameterNameValueType nv;

    IntColumn icol(L"int");
    cols.push_back(icol);

    FloatColumn fcol(L"float");
    cols.push_back(fcol);

    int ng = 10000;
    int nc = 10000;
    int nx = 5000;
	CHPMultiDataData data(fileName);
	data.SetEntryCount(ExpressionMultiDataType, nx, 10, cols);
	data.SetEntryCount(GenotypeMultiDataType, ng, 10, cols);
	data.SetEntryCount(CopyNumberMultiDataType, nc, 10, cols);
	data.SetEntryCount(CytoMultiDataType, nc, 10);
   
    DataSetHeader *dsh = data.GetDataSetHeader(CopyNumberMultiDataType);
    for (int i=0; i<5; ++i)
    {
        ParameterNameValueType param;
        param.SetName(StringUtils::ConvertMBSToWCS(IntToString(i)));
        param.SetValueAscii(IntToString(i));
        dsh->AddNameValParam(param);
    }

	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	ProbeSetMultiDataGenotypeData gn;
	ProbeSetMultiDataCopyNumberData cn;
	ProbeSetMultiDataCytoRegionData cy;
    ProbeSetMultiDataExpressionData ex;

    nv.SetName(L"int");
    nv.SetValueInt32(0);
    ex.metrics.push_back(nv);
    gn.metrics.push_back(nv);
    nv.SetName(L"float");
    nv.SetValueFloat(0);
    ex.metrics.push_back(nv);
    gn.metrics.push_back(nv);

    writer->SeekToDataSet(GenotypeMultiDataType);
    gn.call = 0;
    gn.confidence = 0.0f;
    //for (int i=0; i<ng; i++)
    //{
    //    gn.name = IntToString(i+offset);
	   // writer->WriteEntry(gn);
    //}

    writer->SeekToDataSet(CopyNumberMultiDataType);
    cn.chr = 0;
    cn.position = 0;
    /*for (int i=0; i<nc; i++)
    {
        cn.name = IntToString(i+offset);
	    writer->WriteEntry(cn);
    }*/

    writer->SeekToDataSet(CytoMultiDataType);
    cy.call = 0;
    cy.confidenceScore = 0;
    cy.chr = 0;
    cy.startPosition = 0;
    cy.stopPosition = 0;
    /*for (int i=0; i<nc; i++)
    {
        cy.name = IntToString(i+offset);
	    writer->WriteEntry(cy);
    }*/

	writer->SeekToDataSet(ExpressionMultiDataType);
    ex.quantification = 0;
    for (int i=0; i<nx; i++)
    {
        ex.name = IntToString(i+offset);
        writer->WriteEntry(ex);
    }

	delete writer;

}

void CalvinCHPMultiDataFileBufferWriterTest::testBuffer3()
{
	CreateReferenceFile3(TEST_FILE3a, 0);
	CreateReferenceFile3(TEST_FILE3b, 1001);
    vector<string> chps;
    vector<int> offset;
    chps.push_back(TEST_FILE3a);
    offset.push_back(0);
    chps.push_back(TEST_FILE3b);
    offset.push_back(1001);
    UpdateFile3(chps, offset);
    ChecFile3(TEST_FILE3a, 0);
    ChecFile3(TEST_FILE3b, 1001);
}

void CalvinCHPMultiDataFileBufferWriterTest::ChecFile3(const char *fileName, int offset)
{
    ProbeSetMultiDataExpressionData ex;
	ProbeSetMultiDataGenotypeData gn;
	ProbeSetMultiDataCopyNumberData cn;
	ProbeSetMultiDataCytoRegionData cy;
	CHPMultiDataData data;
	CHPMultiDataFileReader reader;
	reader.SetFilename(fileName);
	reader.Read(data);

    int ng = 10000;
    int nc = 10000;
    int nx = 5000;

    CPPUNIT_ASSERT(data.GetEntryCount(ExpressionMultiDataType) == nx);
	CPPUNIT_ASSERT(data.GetEntryCount(GenotypeMultiDataType) == ng);
	CPPUNIT_ASSERT(data.GetEntryCount(CopyNumberMultiDataType) == nc);
	CPPUNIT_ASSERT(data.GetEntryCount(CytoMultiDataType) == nc);

    for (int i=0; i<nx; i++)
    {
        data.GetExpressionEntry(ExpressionMultiDataType, i, ex);
    	CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.quantification, i+offset, 0.0001f);
	    CPPUNIT_ASSERT(ex.name == IntToString(i+offset));

        CPPUNIT_ASSERT(ex.metrics[0].GetValueInt32() == i+offset);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.metrics[1].GetValueFloat(), 2*i+offset, 0.00001f);
    }
    for (int i=0; i<ng; i++)
    {
        data.GetGenotypeEntry(GenotypeMultiDataType, i, gn);
    	CPPUNIT_ASSERT(gn.call == (i+offset)%4);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.confidence, i+offset, 0.0001f);
	    CPPUNIT_ASSERT(gn.name == IntToString(i+offset));

        CPPUNIT_ASSERT(gn.metrics[0].GetValueInt32() == i+offset);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.metrics[1].GetValueFloat(), 2*i+offset, 0.00001f);
    }
    for (int i=0; i<nc; i++)
    {
        data.GetCopyNumberEntry(CopyNumberMultiDataType, i, cn);
        CPPUNIT_ASSERT(cn.chr == (i+offset)%4);
        CPPUNIT_ASSERT(cn.position == i+offset);
	    CPPUNIT_ASSERT(cn.name == IntToString(i+offset));

        CPPUNIT_ASSERT(cn.metrics[0].GetValueInt32() == i+offset);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(cn.metrics[1].GetValueFloat(), 2*i+offset, 0.00001f);
    }
    for (int i=0; i<nc; i++)
    {
        data.GetCytoEntry(CytoMultiDataType, i, cy);
        CPPUNIT_ASSERT(cy.call == (i+offset)%4);
        CPPUNIT_ASSERT(cy.chr == (i+offset)%10);
        CPPUNIT_ASSERT(cy.startPosition == i);
        CPPUNIT_ASSERT(cy.stopPosition == i+1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(cy.confidenceScore, i+offset, 0.00001f);
	    CPPUNIT_ASSERT(cy.name == IntToString(i+offset));
    }
}

void CalvinCHPMultiDataFileBufferWriterTest::UpdateFile3(std::vector<std::string> &fileNames, vector<int> offset)
{
    ParameterNameValueType nv;
    CHPMultiDataFileBufferWriter writer;
    vector<MultiDataType> dataTypes;
    dataTypes.push_back(ExpressionMultiDataType);
    dataTypes.push_back(GenotypeMultiDataType);
    dataTypes.push_back(CopyNumberMultiDataType);
    dataTypes.push_back(CytoMultiDataType);
    map<MultiDataType, int> maxNameLengths;
    maxNameLengths[ExpressionMultiDataType]=10;
    maxNameLengths[GenotypeMultiDataType]=10;
    maxNameLengths[CopyNumberMultiDataType]=10;
    maxNameLengths[CytoMultiDataType]=10;
    writer.Initialize(&fileNames, dataTypes, maxNameLengths);
    writer.SetMaxBufferSize(102400);

    int ng = 10000;
    int nc = 10000;
    int nx = 5000;

    ProbeSetMultiDataExpressionData ex;
    ProbeSetMultiDataGenotypeData gn;
    ProbeSetMultiDataCopyNumberData cn;
    ProbeSetMultiDataCytoRegionData cy;

    nv.SetName(L"int");
    nv.SetValueInt32(0);
    ex.metrics.push_back(nv);
    gn.metrics.push_back(nv);
    cn.metrics.push_back(nv);
    nv.SetName(L"float");
    nv.SetValueFloat(0);
    ex.metrics.push_back(nv);
    gn.metrics.push_back(nv);
    cn.metrics.push_back(nv);
    for (int i=0; i<nx; i++)
    {
        for (int ifile=0; ifile<(int)fileNames.size(); ifile++)
        {
            ex.name = IntToString(i+offset[ifile]);
            ex.quantification = (float)(i+offset[ifile]);
            ex.metrics[0].SetValueInt32(i+offset[ifile]);
            ex.metrics[1].SetValueFloat(2.0f*i+offset[ifile]);
            writer.WriteMultiDataExpressionEntry(ExpressionMultiDataType, ifile, ex);
        }
    }

    for (int i=0; i<ng; i++)
    {
        for (int ifile=0; ifile<(int)fileNames.size(); ifile++)
        {
            gn.name = IntToString(i+offset[ifile]);
            gn.call = (i+offset[ifile])%4;
            gn.confidence = (float)(i+offset[ifile]);
            gn.metrics[0].SetValueInt32(i+offset[ifile]);
            gn.metrics[1].SetValueFloat(2.0f*i+offset[ifile]);
            writer.WriteMultiDataGenotypeEntry(GenotypeMultiDataType, ifile, gn);
        }
    }
    for (int i=0; i<nc; i++)
    {
        for (int ifile=0; ifile<(int)fileNames.size(); ifile++)
        {
            cn.name = IntToString(i+offset[ifile]);
            cn.chr = ((i+offset[ifile])%4);
            cn.position = (i+offset[ifile]);
            cn.metrics[0].SetValueInt32(i+offset[ifile]);
            cn.metrics[1].SetValueFloat(2.0f*i+offset[ifile]);
            writer.WriteMultiDataCopyNumberEntry(CopyNumberMultiDataType, ifile, cn);
        }
    }
    for (int i=0; i<nc; i++)
    {
        for (int ifile=0; ifile<(int)fileNames.size(); ifile++)
        {
            cy.name = IntToString(i+offset[ifile]);
            cy.call = ((i+offset[ifile])%4);
            cy.chr = ((i+offset[ifile])%10);
            cy.startPosition = i;
            cy.stopPosition = i+1;
            cy.confidenceScore = (float)(i+offset[ifile]);
            writer.WriteMultiDataCytoRegionEntry(CytoMultiDataType, ifile, cy);
        }
    }
    writer.FlushBuffer();
}

