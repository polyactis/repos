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
#include "CalvinCHPMultiDataFileWriter.h"
#include "CalvinCHPMultiDataFileWriterTest.h"
#include "GenericFileReader.h"
#include "ProbeSetMultiDataData.h"

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;

CPPUNIT_TEST_SUITE_REGISTRATION( CHPMultiDataFileWriterTest );

void CHPMultiDataFileWriterTest::setUp() {}

void CHPMultiDataFileWriterTest::tearDown(){}

void CHPMultiDataFileWriterTest::testCreation()
{
	CHPMultiDataData fHdr("CHP_MultiData_file_empty");
	CHPMultiDataFileWriter* w = new CHPMultiDataFileWriter(fHdr);
	CPPUNIT_ASSERT(1);
	delete w;
}

void CHPMultiDataFileWriterTest::WriteTestGeno()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPMultiDataData data("CHP_MultiData_file");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(GenotypeMultiDataType, 2, 10);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);


	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::ProbeSetMultiDataGenotypeData e;

	writer->SeekToDataSet(GenotypeMultiDataType);
	e.name = "abc";
	e.confidence = 10.0f;
	writer->WriteEntry(e);
	e.name = "xyz";
	e.confidence = 20.0f;
	writer->WriteEntry(e);

	delete writer;

	CPPUNIT_ASSERT(1);

    /////////////////////////////

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_file");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 0);

	ParameterNameValueTypeList p = data2.GetAlgParams();
	ParameterNameValueTypeList::iterator it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"an1");
	CPPUNIT_ASSERT(param.GetValueText() == L"av1");

	p = data2.GetSummaryParams();
	it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"sn1");
	CPPUNIT_ASSERT(param.GetValueText() == L"sv1");


    data2.GetGenotypeEntry(GenotypeMultiDataType, 0, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidence, 10.0f, 0.0001f);
	CPPUNIT_ASSERT(e.name == "abc");
	data2.GetGenotypeEntry(GenotypeMultiDataType, 1, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidence, 20.0f, 0.0001f);
	CPPUNIT_ASSERT(e.name == "xyz");

}

void CHPMultiDataFileWriterTest::WriteTestCN()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPMultiDataData data("CHP_MultiData_file_cn");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(CopyNumberMultiDataType, 2, 10);
    data.SetEntryCount(CytoMultiDataType, 2, 10);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);


	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::ProbeSetMultiDataCopyNumberData e;
	affymetrix_calvin_data::ProbeSetMultiDataCytoRegionData c;

	writer->SeekToDataSet(CopyNumberMultiDataType);
	e.name = "abc";
    e.chr = 10;
    e.position = 11;
	writer->WriteEntry(e);
	e.name = "xyz";
    e.chr = 20;
    e.position = 21;
	writer->WriteEntry(e);

	writer->SeekToDataSet(CytoMultiDataType);
	c.name = "abc";
    c.call = 10;
    c.confidenceScore = 11.0f;
    c.chr = 10;
    c.startPosition = 10;
    c.stopPosition = 11;
	writer->WriteEntry(c);
	c.name = "xyz";
    c.call = 20;
    c.confidenceScore = 21.0f;
    c.chr = 20;
    c.startPosition = 20;
    c.stopPosition = 21;
	writer->WriteEntry(c);

	delete writer;

	CPPUNIT_ASSERT(1);

    /////////////////////////////

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_file_cn");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(CopyNumberMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 0);
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 0);

	ParameterNameValueTypeList p = data2.GetAlgParams();
	ParameterNameValueTypeList::iterator it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"an1");
	CPPUNIT_ASSERT(param.GetValueText() == L"av1");

	p = data2.GetSummaryParams();
	it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"sn1");
	CPPUNIT_ASSERT(param.GetValueText() == L"sv1");

    data2.GetCopyNumberEntry(CopyNumberMultiDataType, 0, e);
	CPPUNIT_ASSERT(e.chr == 10);
	CPPUNIT_ASSERT(e.position == 11);
	CPPUNIT_ASSERT(e.name == "abc");
	data2.GetCopyNumberEntry(CopyNumberMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.chr == 20);
	CPPUNIT_ASSERT(e.position == 21);
	CPPUNIT_ASSERT(e.name == "xyz");

    data2.GetCytoEntry(CytoMultiDataType, 0, c);
	CPPUNIT_ASSERT(c.call == 10);
	CPPUNIT_ASSERT(c.chr == 10);
	CPPUNIT_ASSERT(c.startPosition == 10);
	CPPUNIT_ASSERT(c.stopPosition == 11);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(c.confidenceScore, 11.0f, 0.0001f);
	CPPUNIT_ASSERT(c.name == "abc");
	data2.GetCytoEntry(CytoMultiDataType, 1, c);
	CPPUNIT_ASSERT(c.call == 20);
	CPPUNIT_ASSERT(c.chr == 20);
	CPPUNIT_ASSERT(c.startPosition == 20);
	CPPUNIT_ASSERT(c.stopPosition == 21);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(c.confidenceScore, 21.0f, 0.0001f);
	CPPUNIT_ASSERT(c.name == "xyz");
}

void CHPMultiDataFileWriterTest::WriteTestExp()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPMultiDataData data("CHP_MultiData_file");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(ExpressionMultiDataType, 2, 10);

	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::ProbeSetMultiDataExpressionData e;

	writer->SeekToDataSet(ExpressionMultiDataType);
	e.name = "abc";
    e.quantification = 10.0f;
	writer->WriteEntry(e);
	e.name = "xyz";
	e.quantification = 20.0f;
	writer->WriteEntry(e);

	delete writer;

	CPPUNIT_ASSERT(1);

    /////////////////////////////

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_file");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 0);

    data2.GetExpressionEntry(ExpressionMultiDataType, 0, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 10.0f, 0.0001f);
	CPPUNIT_ASSERT(e.name == "abc");
	data2.GetExpressionEntry(ExpressionMultiDataType, 1, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 20.0f, 0.0001f);
	CPPUNIT_ASSERT(e.name == "xyz");

}

void CHPMultiDataFileWriterTest::WriteTestAll()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPMultiDataData data("CHP_MultiData_file");

    vector<ColumnInfo> cols;
    IntColumn icol(L"int");
    cols.push_back(icol);
    FloatColumn fcol(L"float");
    cols.push_back(fcol);

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(ExpressionMultiDataType, 1, 10);
	data.SetEntryCount(GenotypeMultiDataType, 2, 10, cols);

	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::ProbeSetMultiDataExpressionData ex;
    affymetrix_calvin_data::ProbeSetMultiDataGenotypeData gn;

    // expression
	writer->SeekToDataSet(ExpressionMultiDataType);
	ex.name = "ex1";
    ex.quantification = 10.0f;
	writer->WriteEntry(ex);

    // genotype
    writer->SeekToDataSet(GenotypeMultiDataType);
	gn.name = "gn1";
	gn.call = 1;
    gn.confidence = 11.0f;
    gn.metrics.clear();
    param.SetValueInt32(1);
    gn.metrics.push_back(param);
    param.SetValueFloat(2.0f);
    gn.metrics.push_back(param);
	writer->WriteEntry(gn);

	gn.name = "gn2";
	gn.call = 2;
    gn.confidence = 22.0f;;
    gn.metrics.clear();
    param.SetValueInt32(2);
    gn.metrics.push_back(param);
    param.SetValueFloat(3.0f);
    gn.metrics.push_back(param);
	writer->WriteEntry(gn);

	delete writer;

	CPPUNIT_ASSERT(1);

    /////////////////////////////

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_file");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 1);
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 2);

    // expression
    CPPUNIT_ASSERT(data2.GetNumMetricColumns(ExpressionMultiDataType) == 0);
    data2.GetExpressionEntry(ExpressionMultiDataType, 0, ex);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.quantification, 10.0f, 0.0001f);
	CPPUNIT_ASSERT(ex.name == "ex1");

    // genotype
    CPPUNIT_ASSERT(data2.GetNumMetricColumns(GenotypeMultiDataType) == 2);
    CPPUNIT_ASSERT(data2.GetMetricColumnName(GenotypeMultiDataType, 0) == L"int");
    CPPUNIT_ASSERT(data2.GetMetricColumnName(GenotypeMultiDataType, 1) == L"float");

    data2.GetGenotypeEntry(GenotypeMultiDataType, 0, gn);
	CPPUNIT_ASSERT(gn.call == 1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.confidence, 11.0f, 0.0001f);
	CPPUNIT_ASSERT(gn.name == "gn1");
    CPPUNIT_ASSERT(gn.metrics.size() == 2);
    param = gn.metrics[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 1);
    param = gn.metrics[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 2.0f, 0.00001f);

    data2.GetGenotypeEntry(GenotypeMultiDataType, 1, gn);
	CPPUNIT_ASSERT(gn.call == 2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.confidence, 22.0f, 0.0001f);
	CPPUNIT_ASSERT(gn.name == "gn2");
    CPPUNIT_ASSERT(gn.metrics.size() == 2);
    param = gn.metrics[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 2);
    param = gn.metrics[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 3.0f, 0.00001f);
}
