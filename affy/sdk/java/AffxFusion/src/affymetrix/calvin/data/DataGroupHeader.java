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

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.utils.IOUtils;
import affymetrix.portability.UInt;

/**  */
public class DataGroupHeader {

	public DataGroupHeader() {
		clear();
	}

	public DataGroupHeader(String n) {
		clear();
		name = n;
	}

	/** data dataGroup name */
	private String name;

	/** file position of the 1st data dataSet */
	private UInt dataSetPos = new UInt();

	/** file position of the next dataGroup */
	private UInt nextGrpPos = new UInt();

	/** data dataSets in this dataGroup */
	private List<DataSetHeader> dataSetHdrs = new ArrayList<DataSetHeader>();

	/**  */
	public void clear() {
		name = IOUtils.EMPTY;
		dataSetHdrs.clear();
	}

	/**  */
	public void setName(String p) {
		name = p;
	}

	/**  */
	public String getName() {
		return name;
	}

	/** Get the data set count */
	public int getDataSetCnt() {
		if (dataSetHdrs != null) {
			return dataSetHdrs.size();
		}
		return 0;
	}

	/**  */
	public void addDataSetHdr(DataSetHeader p) {
		dataSetHdrs.add(p);
	}

	/**  */
	public DataSetHeader getDataSet(int index) {
		return dataSetHdrs.get(index);
	}

	/**  */
	public List<DataSetHeader> getDataSets() {
		return dataSetHdrs;
	}

	/**
	 * Set the file position of the DataSet header. The value set here is not necessarily the value written to the file.
	 */
	public void setDataSetPos(UInt pos) throws UnsignedOutOfLimitsException {
		dataSetPos.set(pos.toLong());
	}

	/** Get the file position of the DataSet header. */
	public UInt getDataSetPos() {
		return dataSetPos;
	}

	/** Set the file position of the next DataGroup header. */
	public void setNextGroupPos(UInt pos) throws UnsignedOutOfLimitsException {
		nextGrpPos.set(pos.toLong());
	}

	/** Get the file position of the next DataGroup header. */
	public UInt getNextGroupPos() {
		return nextGrpPos;
	}

	/**
	 */
	public DataSetHeader findDataSetHeader(String dataSetName) {
		int n = getDataSetCnt();
		for (int i = 0; i < n; i++) {
			DataSetHeader dph = getDataSet(i);
			if (dataSetName.equals(dph.getName())) {
				return dph;
			}
		}
		return null;
	}
}
