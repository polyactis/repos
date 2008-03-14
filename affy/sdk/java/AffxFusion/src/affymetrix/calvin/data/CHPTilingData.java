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
import java.util.Iterator;
import java.util.List;

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.ParameterNameValue;

/** Holds data associated with tiling array CHP files. */
public class CHPTilingData extends ChpDataBase {

	/** The identifier to identify a tiling CHP file. */
	public static final String CHP_TILING_TYPE = "affymetrix-tiling-analysis";

	/** The identifier for the type of data stored in the file. */
	public static final String TILING_DATA_TYPE = "file_type";

	/** The identifier for the scale of the data. */
	public static final String TILING_PARAM_PREFIX = "Param-";

	/** The identifier for the algorithms tail type. */
	public static final String TILING_NUM_SEQS = "NumberSequences";

	/** The identifier for the algorithm name. */
	public static final String TILING_ALG_NAME = "Algorithm-Name";

	/** The identifier for the algorithm version. */
	public static final String TILING_ALG_VERSION = "Algorithm-Version";

	/** The name of the tiling data group. */
	public static final String CHP_TILING_GROUP = "Tiling Results";

	/** The id for the sequence name. */
	public static final String TILING_SEQ_NAME = "Name";

	/** The id for the sequence group name. */
	public static final String TILING_SEQ_GROUP_NAME = "GroupName";

	/** The id for the sequence version. */
	public static final String TILING_SEQ_VERSION = "Version";

	/** The value to indicate signal values are stored in the CHP file. */
	public static final String TILING_SIGNAL_VALUES = "Signal";

	/** The value to indicate p-values are stored in the CHP file. */
	public static final String TILING_PVALUE_VALUES = "p-value";

	/** column name for position */
	private static final String GENOMIC_POS_COL_NAME = "Genomic Position";

	/** column name for p-value */
	private static final String VAL_COL_NAME = "Result";

	/** keep number of sequences from being read from the header all the time */
	private int cachedNumSequences;

	/** chp data sets */
	private DataSet entries;

	/** Constructor */
	public CHPTilingData() {
		entries = null;
		clear();
	}

	/**
	 * Constructor with file name.
	 * 
	 * @param filename
	 *          The name of the CHP file.
	 */
	public CHPTilingData(String filename) {
		entries = null;
		clear();
		setFilename(filename);
		DataGroupHeader dcHdr = new DataGroupHeader(CHP_TILING_GROUP);
		genericData.getHeader().addDataGroupHdr(dcHdr);
		genericData.getHeader().getGenericDataHdr().setFileTypeId(CHP_TILING_TYPE);
	}

	/** Clears the members. */
	public void clear() {
		entries = null;
		genericData.getHeader().clear();
		cachedNumSequences = -1;
	}

	/**
	 * Sets the file name.
	 * 
	 * @param p
	 *          The name of the CHP file
	 */
	public void setFilename(String p) {
		genericData.getHeader().setFilename(p);
	}

	/**
	 * Gets the file name.
	 * 
	 * @return The file name.
	 */
	public String getFilename() {
		return genericData.getHeader().getFilename();
	}

	/**
	 * Gets the number of sequences.
	 * 
	 * @return The number of sequences.
	 */
	public int getNumberSequences() {
		if (cachedNumSequences == -1) {
			cachedNumSequences = getIntFromGenericHdr(TILING_NUM_SEQS);
		}
		return cachedNumSequences;
	}

	public void setNumberSequences(int value) {
		setIntToGenericHdr(TILING_NUM_SEQS, value);
		cachedNumSequences = value;
	}

	/**
	 * Gets the name of the algorithm.
	 * 
	 * @return The algorithm name.
	 */
	public String getAlgName() {
		return getStringFromGenericHdr(TILING_ALG_NAME);
	}

	public void setAlgName(String value) {
		setStringToGenericHdr(TILING_ALG_NAME, value);
	}

	/**
	 * Gets the algorithm version.
	 * 
	 * @return The version.
	 */
	public String getAlgVersion() {
		return getStringFromGenericHdr(TILING_ALG_VERSION);
	}

	public void setAlgVersion(String value) {
		setStringToGenericHdr(TILING_ALG_VERSION, value);
	}

	/**
	 * gets the algorithm parameters
	 * 
	 * @return The algoirhtm parameters.
	 */
	public List<ParameterNameValue> getAlgParams() {
		return getParams(TILING_PARAM_PREFIX);
	}

	public void addAlgParams(List<ParameterNameValue> params) {
		super.addParams(TILING_PARAM_PREFIX, params);
	}

	/**
	 * Gets the file header.
	 * 
	 * @return The file header.
	 */
	public FileHeader getFileHeader() {
		return genericData.getHeader();
	}

	/**
	 * Gets the generic data object.
	 * 
	 * @return The data object.
	 */
	public GenericData getGenericData() {
		return genericData;
	}

	/**
	 * Gets the sequence data.
	 * 
	 * @return The data associated with the sequence.
	 */
	public TilingSequenceData getTilingSequenceData() {
		TilingSequenceData data = new TilingSequenceData();
		if ((entries != null) && (entries.isOpen() == true)) {
			DataSetHeader hdr = entries.getHeader();
			List<ParameterNameValue> params = hdr.getNameValParameters();
			List<ParameterNameValue> dataParams = new ArrayList<ParameterNameValue>();
			for (int i = 0; i < params.size(); i++) {
				ParameterNameValue param = params.get(i);
				if (param.getName().equals(TILING_SEQ_NAME)) {
					data.setName(param.getValueText());
				}
				else if (param.getName().equals(TILING_SEQ_GROUP_NAME)) {
					data.setGroupName(param.getValueText());
				}
				else if (param.getName().equals(TILING_SEQ_VERSION)) {
					data.setVersion(param.getValueText());
				}
				else {
					dataParams.add(param);
				}
			}
			data.addParameters(dataParams);
		}
		return data;
	}

	public void addTilingSequenceData(int numEntries, TilingSequenceData data) {
		DataGroupHeader dcHdr = genericData.getHeader().getDataGroup(0);
		DataSetHeader dpHdr = new DataSetHeader();
		dpHdr.setRowCnt(numEntries);

		// int nSets = dcHdr.getDataSetCnt();
		// wchar_t name[65];
		// FormatString1(name, 65, L"%d", nSets);
		// dpHdr.setName(name);

		int nSets = dcHdr.getDataSetCnt();
		String s = String.valueOf(nSets);
		if (s.length() < 65) {
			char[] buf = new char[65];

			// copy s into buf
			char[] sBuf = s.toCharArray();
			System.arraycopy(sBuf, 0, buf, 0, sBuf.length);

			dpHdr.setName(new String(buf));
		}
		else if (s.length() == 65) {
			dpHdr.setName(s);
		}
		else {

			throw new RuntimeException("DataSet count is greater than 65 characters in length.");
		}
		ParameterNameValue param1 = new ParameterNameValue();
		param1.setName(TILING_SEQ_NAME);
		param1.setValueText(data.getName());
		dpHdr.addNameValParam(param1);

		ParameterNameValue param2 = new ParameterNameValue();
		param2.setName(TILING_SEQ_GROUP_NAME);
		param2.setValueText(data.getGroupName());
		dpHdr.addNameValParam(param2);

		ParameterNameValue param3 = new ParameterNameValue();
		param3.setName(TILING_SEQ_VERSION);
		param3.setValueText(data.getVersion());
		dpHdr.addNameValParam(param3);

		Iterator<ParameterNameValue> it = data.getParameters().iterator();
		while (it.hasNext()) {
			dpHdr.addNameValParam(it.next());
		}
		addColumns(dpHdr);
		dcHdr.addDataSetHdr(dpHdr);
	}

	/**
	 * Gets the number of entries in a tiling sequence.
	 * 
	 * @param index
	 *          The sequence index.
	 * @return The number of entries in the sequence.
	 */
	public int getTilingSequenceEntryCount(int index) {
		DataGroupHeader dcHdr = genericData.getHeader().getDataGroup(0);
		DataSetHeader dpHdr = dcHdr.getDataSet(index);
		return dpHdr.getRowCnt();
	}

	/**
	 * Opens a group for reading.
	 * 
	 * @param index
	 *          The index to the sequence.
	 */
	public void openTilingSequenceDataSet(int index) {
		if (entries != null) {
			entries.delete();
		}
		try {
			entries = genericData.getDataSet(0, index);
			if (entries != null) {
				entries.open();
			}
		}
		catch (Throwable t) {
		}
	}

	/**
	 * Returns the entry for the given row. The data set must be open.
	 * 
	 * @param row
	 *          The row index.
	 * @return The entry.
	 */
	public CHPTilingEntry getTilingSequenceEntry(int row) throws UnsignedOutOfLimitsException {
		if ((entries != null) && (entries.isOpen() == true)) {
			CHPTilingEntry e = new CHPTilingEntry();
			e.setPosition(entries.getDataInt(row, 0));
			e.setValue(entries.getDataFloat(row, 1));
			return e;
		}
		return null;
	}

	private void addColumns(DataSetHeader hdr) {
		hdr.addUIntColumn(GENOMIC_POS_COL_NAME); // genomic position - unsigned int 32
		hdr.addFloatColumn(VAL_COL_NAME); // value - float
	}
}
