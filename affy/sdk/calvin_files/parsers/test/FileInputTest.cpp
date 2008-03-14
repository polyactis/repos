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
#include <istream>
#include <math.h>
#include <string>
//
#include "FileInput.h"
#include "FileInputTest.h"

using namespace std;
using namespace affymetrix_calvin_io;

#ifdef _MSC_VER
#pragma warning(disable: 4996)
#endif


CPPUNIT_TEST_SUITE_REGISTRATION( FileInputTest );

static bool CompareFloats(float f1, float f2)
{
	const float EPS = 0.0000001f;
	return (fabs(f1-f2) < EPS);
}

#define TEST_DATA_FILE "../data/test.file"

void FileInputTest::setUp()
{
}

void FileInputTest::tearDown()
{
}

void FileInputTest::testdefine_SIZES()
{
	CPPUNIT_ASSERT( sizeof(int8_t) == 1 );
	CPPUNIT_ASSERT( sizeof(int16_t) == 2 );
	CPPUNIT_ASSERT( sizeof(int32_t) == 4 );
	CPPUNIT_ASSERT( sizeof(int64_t) == 8 );
	CPPUNIT_ASSERT( sizeof(u_int8_t) == 1 );
	CPPUNIT_ASSERT( sizeof(u_int16_t) == 2 );
	CPPUNIT_ASSERT( sizeof(u_int32_t) == 4 );
	CPPUNIT_ASSERT( sizeof(u_int64_t) == 8 );
}

void FileInputTest::testmethod_ReadInt8()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".char");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	int8_t val=0;
	val = FileInput::ReadInt8(instr);
	CPPUNIT_ASSERT(val == 11 );
	val = FileInput::ReadInt8(instr);
	CPPUNIT_ASSERT(val == 22 );
	val = FileInput::ReadInt8(instr);
	CPPUNIT_ASSERT(val == 33 );
	instr.close();
}

void FileInputTest::testmethod_ReadInt16()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".short");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	int16_t val=0;
	val = FileInput::ReadInt16(instr);
	CPPUNIT_ASSERT(val == -21 );
	val = FileInput::ReadInt16(instr);
	CPPUNIT_ASSERT(val == 21 );
	val = FileInput::ReadInt16(instr);
	CPPUNIT_ASSERT(val == -31 );
	instr.close();
}

void FileInputTest::testmethod_ReadInt32()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".int");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	int32_t val=0;
	val = FileInput::ReadInt32(instr);
	CPPUNIT_ASSERT(val == 10 );
	val = FileInput::ReadInt32(instr);
	CPPUNIT_ASSERT(val == 20 );
	val = FileInput::ReadInt32(instr);
	CPPUNIT_ASSERT(val == 30 );
	instr.close();
}

void FileInputTest::testmethod_ReadUInt8()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".uchar");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	u_int8_t val=0;
	val = FileInput::ReadUInt8(instr);
	CPPUNIT_ASSERT(val == 11 );
	val = FileInput::ReadUInt8(instr);
	CPPUNIT_ASSERT(val == 22 );
	val = FileInput::ReadUInt8(instr);
	CPPUNIT_ASSERT(val == 33 );
	instr.close();
}

void FileInputTest::testmethod_ReadUInt16()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".ushort");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	u_int16_t val=0;
	val = FileInput::ReadUInt16(instr);
	CPPUNIT_ASSERT(val == 45 );
	val = FileInput::ReadUInt16(instr);
	CPPUNIT_ASSERT(val == 54 );
	val = FileInput::ReadUInt16(instr);
	CPPUNIT_ASSERT(val == 66 );
	instr.close();
}

void FileInputTest::testmethod_ReadUInt32()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".uint");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	u_int32_t val=0;
	val = FileInput::ReadUInt32(instr);
	CPPUNIT_ASSERT(val == 101 );
	val = FileInput::ReadUInt32(instr);
	CPPUNIT_ASSERT(val == 202 );
	val = FileInput::ReadUInt32(instr);
	CPPUNIT_ASSERT(val == 303 );
	instr.close();
}

void FileInputTest::testmethod_ReadInt8_stream()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".char");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	int8_t strm[3];
	int8_t val=0;
	instr.read((char *)strm, sizeof(strm[0])*3);
	char *p = (char *)strm;
	val = FileInput::ReadInt8(p);
	CPPUNIT_ASSERT(val == 11 );
	val = FileInput::ReadInt8(p);
	CPPUNIT_ASSERT(val == 22 );
	val = FileInput::ReadInt8(p);
	CPPUNIT_ASSERT(val == 33 );
	instr.close();
}

void FileInputTest::testmethod_ReadInt16_stream()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".short");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	int16_t strm[3];
	int16_t val=0;
	instr.read((char *)strm, sizeof(strm[0])*3);
	char *p = (char *)strm;
	val = FileInput::ReadInt16(p);
	CPPUNIT_ASSERT(val == -21 );
	val = FileInput::ReadInt16(p);
	CPPUNIT_ASSERT(val == 21 );
	val = FileInput::ReadInt16(p);
	CPPUNIT_ASSERT(val == -31 );
	instr.close();
}

void FileInputTest::testmethod_ReadInt32_stream()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".int");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	int32_t strm[3];
	int32_t val=0;
	instr.read((char *)strm, sizeof(strm[0])*3);
	char *p = (char *)strm;
	val = FileInput::ReadInt32(p);
	CPPUNIT_ASSERT(val == 10 );
	val = FileInput::ReadInt32(p);
	CPPUNIT_ASSERT(val == 20 );
	val = FileInput::ReadInt32(p);
	CPPUNIT_ASSERT(val == 30 );
	instr.close();
}

void FileInputTest::testmethod_ReadUInt8_stream()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".uchar");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	u_int8_t strm[3];
	u_int8_t val=0;
	instr.read((char *)strm, sizeof(strm[0])*3);
	char *p = (char *)strm;
	val = FileInput::ReadUInt8(p);
	CPPUNIT_ASSERT(val == 11 );
	val = FileInput::ReadUInt8(p);
	CPPUNIT_ASSERT(val == 22 );
	val = FileInput::ReadUInt8(p);
	CPPUNIT_ASSERT(val == 33 );
	instr.close();
}

void FileInputTest::testmethod_ReadUInt16_stream()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".ushort");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	u_int16_t strm[3];
	u_int16_t val=0;
	instr.read((char *)strm, sizeof(strm[0])*3);
	char *p = (char *)strm;
	val = FileInput::ReadUInt16(p);
	CPPUNIT_ASSERT(val == 45 );
	val = FileInput::ReadUInt16(p);
	CPPUNIT_ASSERT(val == 54 );
	val = FileInput::ReadUInt16(p);
	CPPUNIT_ASSERT(val == 66 );
	instr.close();
}

void FileInputTest::testmethod_ReadUInt32_stream()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".uint");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	u_int32_t strm[3];
	u_int32_t val=0;
	instr.read((char *)strm, sizeof(strm[0])*3);
	char *p = (char *)strm;
	val = FileInput::ReadUInt32(p);
	CPPUNIT_ASSERT(val == 101 );
	val = FileInput::ReadUInt32(p);
	CPPUNIT_ASSERT(val == 202 );
	val = FileInput::ReadUInt32(p);
	CPPUNIT_ASSERT(val == 303 );
	instr.close();
}

void FileInputTest::testmethod_ReadFloat()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".float");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	float val=0.0f;
	val = FileInput::ReadFloat(instr);
	CPPUNIT_ASSERT( CompareFloats(val, 321.123f) );
	val = FileInput::ReadFloat(instr);
	CPPUNIT_ASSERT( CompareFloats(val, 123.123f) );
	val = FileInput::ReadFloat(instr);
	CPPUNIT_ASSERT( CompareFloats(val, 543.56f) );
	instr.close();
}

void FileInputTest::testmethod_ReadFloat_stream()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".float");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	float strm[3];
	float val=0.0f;
	instr.read((char *)strm, sizeof(strm[0])*3);
	char *p = (char *)strm;
	val = FileInput::ReadFloat(p);
	CPPUNIT_ASSERT( CompareFloats(val, 321.123f) );
	val = FileInput::ReadFloat(p);
	CPPUNIT_ASSERT( CompareFloats(val, 123.123f) );
	val = FileInput::ReadFloat(p);
	CPPUNIT_ASSERT( CompareFloats(val, 543.56f) );
	instr.close();
}

void FileInputTest::testmethod_ReadString8()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".string8");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	string val;
	val = FileInput::ReadString8(instr);
	CPPUNIT_ASSERT( val == string("string") );
	val = FileInput::ReadString8(instr);
	CPPUNIT_ASSERT( val == string("test") );
	val = FileInput::ReadString8(instr);
	CPPUNIT_ASSERT( val == string("case") );
	instr.close();
}

void FileInputTest::testmethod_ReadString8_stream()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".string8");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	instr.seekg(0, ios::end);
	int flen = instr.tellg();
	instr.seekg(0, ios::beg);
	char *strm = new char[flen+1];
	char *s = strm;
	instr.read(strm, flen);
	string sval;
	sval = FileInput::ReadString8(s);
	CPPUNIT_ASSERT( sval == string("string") );
	sval = FileInput::ReadString8(s);
	CPPUNIT_ASSERT( sval == string("test") );
	sval = FileInput::ReadString8(s);
	CPPUNIT_ASSERT( sval == string("case") );
	instr.close();
	delete[] strm;
}

void FileInputTest::testmethod_ReadString8_fixedlen()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".string8");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	string val;
	int32_t ival;
	ival = FileInput::ReadInt32(instr);
	val = FileInput::ReadString8(instr, ival);
	CPPUNIT_ASSERT( val == string("string") );
	ival = FileInput::ReadInt32(instr);
	val = FileInput::ReadString8(instr, ival);
	CPPUNIT_ASSERT( val == string("test") );
	ival = FileInput::ReadInt32(instr);
	val = FileInput::ReadString8(instr, ival);
	CPPUNIT_ASSERT( val == string("case") );
	instr.close();
}

void FileInputTest::testmethod_ReadString8_fixedlen_stream()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".string8");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	instr.seekg(0, ios::end);
	int flen = instr.tellg();
	instr.seekg(0, ios::beg);
	char *strm = new char[flen+1];
	char *s = strm;
	instr.read(strm, flen);
	string sval;
	int32_t ival;
	ival = FileInput::ReadInt32(s);
	sval = FileInput::ReadString8(s, ival);
	CPPUNIT_ASSERT( sval == string("string") );
	ival = FileInput::ReadInt32(s);
	sval = FileInput::ReadString8(s, ival);
	CPPUNIT_ASSERT( sval == string("test") );
	ival = FileInput::ReadInt32(s);
	sval = FileInput::ReadString8(s, ival);
	CPPUNIT_ASSERT( sval == string("case") );
	instr.close();
	delete[] strm;
}

void FileInputTest::testmethod_ReadString16()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".string16");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	wstring val;
	val = FileInput::ReadString16(instr);
	CPPUNIT_ASSERT( val == L"string" );
	val = FileInput::ReadString16(instr);
	CPPUNIT_ASSERT( val == L"test" );
	val = FileInput::ReadString16(instr);
	CPPUNIT_ASSERT( val == L"case" );
	instr.close();
}

void FileInputTest::testmethod_ReadString16_stream()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".string16");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	instr.seekg(0, ios::end);
	int flen = instr.tellg();
	instr.seekg(0, ios::beg);
	char *strm = new char[flen+1];
	char *s = strm;
	instr.read(strm, flen);
	wstring sval;
	sval = FileInput::ReadString16(s);
	CPPUNIT_ASSERT( sval == L"string" );
	sval = FileInput::ReadString16(s);
	CPPUNIT_ASSERT( sval == L"test" );
	sval = FileInput::ReadString16(s);
	CPPUNIT_ASSERT( sval == L"case" );
	instr.close();
	delete[] strm;
}

void FileInputTest::testmethod_ReadString16_fixedlen()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".string16");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	wstring val;
	int32_t ival;
	ival = FileInput::ReadInt32(instr);
	val = FileInput::ReadString16(instr, ival);
	CPPUNIT_ASSERT( val == L"string" );
	ival = FileInput::ReadInt32(instr);
	val = FileInput::ReadString16(instr, ival);
	CPPUNIT_ASSERT( val == L"test" );
	ival = FileInput::ReadInt32(instr);
	val = FileInput::ReadString16(instr, ival);
	CPPUNIT_ASSERT( val == L"case" );
	instr.close();
}

void FileInputTest::testmethod_ReadString16_fixedlen_stream()
{
	ifstream instr;
	string file = TEST_DATA_FILE + string(".string16");
	instr.open(file.c_str(), ios::in | ios::binary);
	CPPUNIT_ASSERT( instr );
	instr.seekg(0, ios::end);
	int flen = instr.tellg();
	instr.seekg(0, ios::beg);
	char *strm = new char[flen+1];
	char *s = strm;
	instr.read(strm, flen);
	wstring sval;
	int32_t ival;
	ival = FileInput::ReadInt32(s);
	sval = FileInput::ReadString16(s, ival);
	CPPUNIT_ASSERT( sval == L"string" );
	ival = FileInput::ReadInt32(s);
	sval = FileInput::ReadString16(s, ival);
	CPPUNIT_ASSERT( sval == L"test" );
	ival = FileInput::ReadInt32(s);
	sval = FileInput::ReadString16(s, ival);
	CPPUNIT_ASSERT( sval == L"case" );
	instr.close();
	delete[] strm;
}

