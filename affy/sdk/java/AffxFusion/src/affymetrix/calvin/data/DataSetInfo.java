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

package affymetrix.calvin.data;

import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.data.CHPMultiDataData.MultiDataType;

/** Holds data set information for a multi-data chp file data set. */
public class DataSetInfo {

	/** The data type. */
	private MultiDataType dataType;

	/** chp data sets */
	private DataSet entries;

	/** The maximum length of a probe set name. */
	private int maxProbeSetName;

	/** The data set index. */
	private int dataSetIndex;

	/** An array of extra metric columns. */
	private List<ColumnInfo> metricColumns = new ArrayList<ColumnInfo>();

	/** constructor */
	public DataSetInfo() {
		entries = null;
		maxProbeSetName = -1;
		dataSetIndex = -1;
	}

	public int getDataSetIndex() {
		return dataSetIndex;
	}

	public void setDataSetIndex(int dataSetIndex) {
		this.dataSetIndex = dataSetIndex;
	}

	public MultiDataType getDataType() {
		return dataType;
	}

	public void setDataType(MultiDataType dataType) {
		this.dataType = dataType;
	}

	public DataSet getEntries() {
		return entries;
	}

	public void setEntries(DataSet entries) {
		this.entries = entries;
	}

	public int getMaxProbeSetName() {
		return maxProbeSetName;
	}

	public void setMaxProbeSetName(int maxProbeSetName) {
		this.maxProbeSetName = maxProbeSetName;
	}

	public void addMetricColumn(ColumnInfo mColumn) {
		this.metricColumns.add(mColumn);
	}

	public void setMetricColumns(List<ColumnInfo> mColumns) {
		this.metricColumns.addAll(mColumns);
	}

	public List<ColumnInfo> getMetricColumns() {
		return metricColumns;
	}

	public void clearMetricColumns() {
		this.metricColumns.clear();
	}
}
