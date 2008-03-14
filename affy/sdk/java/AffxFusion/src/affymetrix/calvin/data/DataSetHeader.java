/////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 2.1 of the License,
// or (at your option) any later version.
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
/////////////////////////////////////////////////////////////////

package affymetrix.calvin.data;

import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.portability.UInt;

/** The DataSetHeader information container class. */
public class DataSetHeader {

	/** total rows in the dataGroup */
	private int rowCount;

	/** data dataGroup name */
	private String name;

	/** name/value pairs */
	private List<ParameterNameValue> nameValParams = new ArrayList<ParameterNameValue>();

	/** column information */
	private List<ColumnInfo> columnTypes = new ArrayList<ColumnInfo>();

	/** file position of the start of the dataSet header */
	private UInt headerStartFilePos;

	/** file position of the start of the data */
	private UInt dataStartFilePos;

	/** file position of the next dataSet header */
	private UInt nextSetFilePos;

	/** Constructs a new DataSetHeader class. */
	public DataSetHeader() {
		clear();
	}

	public void clear() {
		rowCount = 0;
		name = "";
		nameValParams.clear();
		columnTypes.clear();
		headerStartFilePos = UInt.ZERO;
		dataStartFilePos = UInt.ZERO;
		nextSetFilePos = UInt.ZERO;
	}

	/**  */
	public int getDataSize() {
		return getRowSize() * rowCount;
	}

	/**  */
	public int getRowSize() {
		int result = 0;
		int sz = getColumnCnt();
		for (int i = 0; i < sz; i++) {
			result += getColumnInfo(i).getSize();
		}
		return result;
	}

	/**  */
	public void setName(String p) {
		name = p;
	}

	/**  */
	public String getName() {
		return name;
	}

	/**  */
	public int getNameValParamCnt() {
		return nameValParams.size();
	}

	/**  */
	public void addNameValParam(ParameterNameValue p) {
		nameValParams.add(p);
	}

	/**  */
	public List<ParameterNameValue> getNameValParameters() {
		return nameValParams;
	}

	/**
	 * Finds a ParameterNameValueType by name in the nameValPairs collection
	 * 
	 * @param name
	 *          The name of the NameValPair to find
	 * @return Reference to a ParameterNameValue if found, null otherwise.
	 */
	public ParameterNameValue findNameValParam(String name) {
		int n = nameValParams.size();
		for (int i = 0; i < n; i++) {
			ParameterNameValue param = nameValParams.get(i);
			if (name.equals(param.getName())) {
				return param;
			}
		}
		return null;
	}

	/**  */
	public void addColumn(ColumnInfo colInfo) {
		columnTypes.add(colInfo);
	}

	/**  */
	public void addIntColumn(String name) {
		columnTypes.add(new IntColumn(name));
	}

	/**  */
	public void addUIntColumn(String name) {
		columnTypes.add(new UIntColumn(name));
	}

	/**  */
	public void addShortColumn(String name) {
		columnTypes.add(new ShortColumn(name));
	}

	/**  */
	public void addUShortColumn(String name) {
		columnTypes.add(new UShortColumn(name));
	}

	/**  */
	public void addByteColumn(String name) {
		columnTypes.add(new ByteColumn(name));
	}

	/**  */
	public void addUByteColumn(String name) {
		columnTypes.add(new UByteColumn(name));
	}

	/**  */
	public void addFloatColumn(String name) {
		columnTypes.add(new FloatColumn(name));
	}

	/**
	 * @param len
	 *          Maximum number of char in string
	 */
	public void addAsciiColumn(String name, int len) {
		columnTypes.add(new ASCIIColumn(name, len));
	}

	/**
	 * @param len
	 *          Maximum number of wchar_t in string
	 */
	public void addUnicodeColumn(String name, int len) {
		columnTypes.add(new UnicodeColumn(name, len));
	}

	/**  */
	public ColumnInfo getColumnInfo(int index) {
		if (columnTypes.size() > 0) {
			return columnTypes.get(index);
		}
		return null;
	}

	/**  */
	public int getRowCnt() {
		return rowCount;
	}

	/**  */
	public void setRowCnt(int p) {
		rowCount = p;
	}

	/**  */
	public int getColumnCnt() {
		return columnTypes.size();
	}

	/**
	 * Set the file position of the start of the DataSet header. The value set here is not written to the file.
	 */
	public void setHeaderStartFilePos(UInt pos) {
		headerStartFilePos = pos;
	}

	/** Get the file position of the start of the DataSet header. */
	public UInt getHeaderStartFilePos() {
		return headerStartFilePos;
	}

	/**
	 * Set the file position of the start of the DataSet data. The value set here is not written to the file.
	 */
	public void setDataStartFilePos(UInt pos) {
		dataStartFilePos = pos;
	}

	/** Get the file position of the start of the DataSet data. */
	public UInt getDataStartFilePos() {
		return dataStartFilePos;
	}

	/** Set the file position of the next DataSet header. */
	public void setNextSetFilePos(UInt pos) {
		nextSetFilePos = pos;
	}

	/** Get the file position of the next DataSet header. */
	public UInt getNextSetFilePos() {
		return nextSetFilePos;
	}

}
