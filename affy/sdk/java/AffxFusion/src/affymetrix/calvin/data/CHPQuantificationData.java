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

import java.io.IOException;
import java.util.List;

import affymetrix.calvin.data.ColumnInfo.DataSetColumnTypes;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.AffymetrixParameterConsts;
import affymetrix.calvin.parameter.ParameterNameValue;

/** Holds data associated with quantification array CHP files. */
public class CHPQuantificationData extends ChpDataBase {

	/** The identifier to identify a quantification CHP file. */
	public static final String CHP_QUANTIFICATION_TYPE = "affymetrix-quantification-analysis";

	/** The group and data set and column name to store the quantification data. */
	public static final String QUANTIFICATION_NAME = "Quantification";

	/** The column name for the probe set name. */
	public static final String QUANTIFICATION_PROBE_SET_NAME = "ProbesetName";

	/** The column name for the probe set id. */
	public static final String QUANTIFICATION_PROBE_SET_ID = "ProbeSetId";

	// Flag indicating if the probe set names were stored in wide character format.
	private DataSetColumnTypes firstColumnType;

	// chp data sets
	private DataSet entries;

	// The maximum length of a probe set name.
	private int maxProbeSetName;

	/** Constructor */
	public CHPQuantificationData() {
		entries = null;
		maxProbeSetName = -1;
		firstColumnType = DataSetColumnTypes.UnicodeCharColType;
		clear();
	}

	/**
	 * ructor with file name.
	 * 
	 * @param filename
	 *          The name of the CHP file.
	 */
	public CHPQuantificationData(String filename) {
		this();
		setFilename(filename);
		DataGroupHeader dcHdr = new DataGroupHeader(QUANTIFICATION_NAME);
		genericData.getHeader().addDataGroupHdr(dcHdr);
		genericData.getHeader().getGenericDataHdr().setFileTypeId(CHP_QUANTIFICATION_TYPE);
	}

	/** The maximum length of a probe set name. */
	public int getMaxProbeSetName() {
		return maxProbeSetName;
	}

	/** Clears the members. */
	public void clear() {
		entries = null;
		genericData.getHeader().clear();
	}

	/**
	 * sets the file name.
	 * 
	 * @param p
	 *          The name of the CHP file
	 */
	public void setFilename(String p) {
		genericData.getHeader().setFilename(p);
	}

	/**
	 * gets the file name.
	 * 
	 * @return The file name.
	 */
	public String getFilename() {
		return genericData.getHeader().getFilename();
	}

	/** sets the array type */
	public String getArrayType() {
		return getStringFromGenericHdr(AffymetrixParameterConsts.ARRAY_TYPE_PARAM_NAME);
	}

	public void setArrayType(String value) {
		setStringToGenericHdr(AffymetrixParameterConsts.ARRAY_TYPE_PARAM_NAME, value,
				AffymetrixParameterConsts.ARRAY_TYPE_MAX_LEN);
	}

	/** gets the number of entries (probe sets) */
	public int getEntryCount() {
		DataGroupHeader dcHdr = genericData.getHeader().getDataGroup(0);
		DataSetHeader dpHdr = dcHdr.getDataSet(0);
		return dpHdr.getRowCnt();
	}

	public void setEntryCount(int ln, int maxln) {
		firstColumnType = DataSetColumnTypes.ASCIICharColType;
		maxProbeSetName = maxln;
		DataGroupHeader grpHdr = genericData.getHeader().getDataGroup(0);
		DataSetHeader setHdr = new DataSetHeader();
		setHdr.setRowCnt(ln);
		setHdr.setName(QUANTIFICATION_NAME);
		addColumns(setHdr, false);
		grpHdr.addDataSetHdr(setHdr);
	}

	public void setEntryCount(int ln) {
		firstColumnType = DataSetColumnTypes.IntColType;
		maxProbeSetName = -1;
		DataGroupHeader grpHdr = genericData.getHeader().getDataGroup(0);
		DataSetHeader setHdr = new DataSetHeader();
		setHdr.setRowCnt(ln);
		setHdr.setName(QUANTIFICATION_NAME);
		addColumns(setHdr, true);
		grpHdr.addDataSetHdr(setHdr);
	}

	/**
	 * gets the name of the algorithm.
	 * 
	 * @return The algorithm name.
	 */
	public String getAlgName() {
		return getStringFromGenericHdr(AffymetrixParameterConsts.ALGORITHM_NAME_PARAM_NAME);
	}

	public void setAlgName(String value) {
		setStringToGenericHdr(AffymetrixParameterConsts.ALGORITHM_NAME_PARAM_NAME, value);
	}

	/**
	 * gets the algorithm version.
	 * 
	 * @return The version.
	 */
	public String getAlgVersion() {
		return getStringFromGenericHdr(AffymetrixParameterConsts.ALG_VERSION_PARAM_NAME);
	}

	public void setAlgVersion(String value) {
		setStringToGenericHdr(AffymetrixParameterConsts.ALG_VERSION_PARAM_NAME, value);
	}

	/**
	 * gets the algorithm parameters
	 * 
	 * @return The algoirhtm parameters.
	 */
	public List<ParameterNameValue> getAlgParams() {
		return super.getParams(AffymetrixParameterConsts.ALGORITHM_PARAM_NAME_PREFIX);
	}

	/**
	 * gets the summary parameters
	 * 
	 * @return The summary parameters.
	 */
	public List<ParameterNameValue> getSummaryParams() {
		return getParams(AffymetrixParameterConsts.CHIP_SUMMARY_PARAM_NAME_PREFIX);
	}

	public void addSummaryParams(List<ParameterNameValue> params) {
		super.addParams(AffymetrixParameterConsts.CHIP_SUMMARY_PARAM_NAME_PREFIX, params);
	}

	public void addAlgParams(List<ParameterNameValue> params) {
		super.addParams(AffymetrixParameterConsts.ALGORITHM_PARAM_NAME_PREFIX, params);
	}

	/**
	 * gets the file header.
	 * 
	 * @return The file header.
	 */
	public FileHeader getFileHeader() {
		return genericData.getHeader();
	}

	/**
	 * gets the generic data object.
	 * 
	 * @return The data object.
	 */
	public GenericData getGenericData() {
		return genericData;
	}

	/**
	 * gets the sequence data.
	 * 
	 * @param index
	 *          The row index.
	 * @return The quantification value.
	 */
	public ProbeSetQuantificationData getQuantificationEntry(int index) throws IOException, UnsignedOutOfLimitsException {
		ProbeSetQuantificationData entry = null;
		openQuantificationDataSet();
		if ((entries != null) && entries.isOpen()) {
			entry = new ProbeSetQuantificationData();
			if (firstColumnType == DataSetColumnTypes.ASCIICharColType) {
				entry.setName(entries.getDataString8(index, 0));
			}
			else if (firstColumnType == DataSetColumnTypes.UnicodeCharColType) {
				entry.setName(entries.getDataString16(index, 0));
			}
			else if (firstColumnType == DataSetColumnTypes.IntColType) {
				entry.setId(entries.getDataInt(index, 0));
			}
			entry.setQuantification(entries.getDataFloat(index, 1));
		}
		return entry;
	}

	/** Opens a group for reading. */
	public void openQuantificationDataSet() {
		if (entries == null) {
			try {
				entries = genericData.getDataSet(0, 0);
				if (entries != null) {
					entries.open();
					firstColumnType = entries.getHeader().getColumnInfo(0).getColumnType();
				}
			}
			catch (Throwable t) {
				entries = null;
			}
		}
	}

	private void addColumns(DataSetHeader hdr, boolean keyIsID) {
		if (keyIsID) {
			hdr.addIntColumn(QUANTIFICATION_PROBE_SET_ID);
		}
		else {
			hdr.addAsciiColumn(QUANTIFICATION_PROBE_SET_NAME, maxProbeSetName);
		}
		hdr.addFloatColumn(QUANTIFICATION_NAME);
	}
}
