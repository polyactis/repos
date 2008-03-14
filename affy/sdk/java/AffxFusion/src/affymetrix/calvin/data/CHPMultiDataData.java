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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.AffymetrixParameterConsts;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.portability.UByte;
import affymetrix.portability.UInt;
import affymetrix.portability.UShort;

/** Holds data associated with a multi-data CHP files. */
public class CHPMultiDataData extends ChpDataBase {
	/** The identifier to identify a multi-data CHP file. */
	public static final String CHP_MULTI_DATA_TYPE = "affymetrix-multi-data-type-analysis";

	/** The group and data set to store the generic genotype data. */
	private static final String MULTI_DATA_NAME = "MultiData";

	/** The data set names for a multi-data chp file */
	public static final String[] DATA_SET_NAMES = { "Expression", "ExpressionControl", "Genotype", "GenotypeControl",
			"CopyNumber", "Cyto" };

	// An enumerant to store the types of data stored in the file.
	public enum MultiDataType {
		ExpressionMultiDataType, ExpressionControlMultiDataType, GenotypeMultiDataType, GenotypeControlMultiDataType, CopyNumberMultiDataType, CytoMultiDataType
	};

	/** The column name for the probe set name. */
	private static final String PROBE_SET_NAME = "ProbeSetName";

	/** The column name for the call. */
	private static final String GENOTYPE_CALL_NAME = "Call";

	/** The column name for the confidence value. */
	private static final String GENOTYPE_CONFIDENCE_NAME = "Confidence";

	/** The column name for the quantification value. */
	private static final String EXPRESSION_QUANTIFICATION_NAME = "Quantification";

	/** The column name for the chromosome name. */
	private static final String COPY_NUMBER_CHR_NAME = "Chromosome";

	/* ! The column name for the physical position of the SNP. */
	private static final String COPY_NUMBER_POSITION_NAME = "Position";

	/* ! The column name for the call value. */
	private static final String CYTO_CALL_NAME = "Call";

	/* ! The column name for the confidence value. */
	private static final String CYTO_CONFIDENCE_NAME = "Confidence";

	/* ! The column name for the chromosome name. */
	private static final String CYTO_CHR_NAME = "Chromosome";

	/* ! The column name for the physical position of the cyto region. */
	private static final String CYTO_START_POSITION_NAME = "StartPosition";

	/* ! The column name for the physical position of the cyto region. */
	private static final String CYTO_STOP_POSITION_NAME = "StopPosition";

	/* ! The column name for the cyto region. */
	private static final String CYTO_REGION_NAME = "Region";

	/** a map of chp data sets. */
	public Map<MultiDataType, DataSetInfo> dataSetInfo = new HashMap<MultiDataType, DataSetInfo>();

	/** Constructor */
	public CHPMultiDataData() {
		clear();
	}

	/**
	 * Constructor with file name.
	 * 
	 * @param filename
	 *          The name of the CHP file.
	 */
	public CHPMultiDataData(String filename) {
		this();
		setFilename(filename);
		DataGroupHeader dcHdr = new DataGroupHeader(MULTI_DATA_NAME);
		genericData.getHeader().addDataGroupHdr(dcHdr);
		genericData.getHeader().getGenericDataHdr().setFileTypeId(CHP_MULTI_DATA_TYPE);
	}

	/** Get the data set information map */
	public Map<MultiDataType, DataSetInfo> getDataSetInfo() {
		return dataSetInfo;
	}

	/**
	 * The maximum length of a probe set name.
	 * 
	 * @param dataType
	 *          The data type
	 */
	public int getMaxProbeSetName(MultiDataType dataType) {
		DataSetInfo info = openMultiDataDataSet(dataType);
		if (info != null) {
			return info.getMaxProbeSetName();
		}
		return 0;
	}

	/** Clears the members. */
	public void clear() {
		dataSetInfo.clear();
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

	/**
	 * gets the number of entries (probe sets)
	 * 
	 * @param dataType
	 *          The data type
	 */
	public int getEntryCount(MultiDataType dataType) {
		DataSetHeader h = getDataSetHeader(dataType);
		return (h == null ? 0 : h.getRowCnt());
	}

	public void setEntryCount(MultiDataType dataType, int ln, int maxln) {
		List<ColumnInfo> empty = new ArrayList<ColumnInfo>();
		setEntryCount(dataType, ln, maxln, empty);
	}

	public void setEntryCount(MultiDataType dataType, int ln, int maxln, List<ColumnInfo> columns) {
		DataSetInfo info = new DataSetInfo();
		info.setMaxProbeSetName(maxln);
		info.setMetricColumns(columns);
		info.setEntries(null);
		info.setDataType(dataType);
		info.setDataSetIndex(dataSetInfo.size());
		dataSetInfo.put(dataType, info);

		DataSetHeader dsHdr = new DataSetHeader();
		dsHdr.setRowCnt(ln);
		dsHdr.setName(DATA_SET_NAMES[dataType.ordinal()]);
		addColumns(info, dsHdr);

		DataGroupHeader dgHdr = genericData.getHeader().getDataGroup(0);
		dgHdr.addDataSetHdr(dsHdr);
	}

	private void addColumns(DataSetInfo info, DataSetHeader hdr) {
		switch (info.getDataType()) {
		case ExpressionMultiDataType:
		case ExpressionControlMultiDataType:
			hdr.addAsciiColumn(PROBE_SET_NAME, info.getMaxProbeSetName());
			hdr.addFloatColumn(EXPRESSION_QUANTIFICATION_NAME);
			break;

		case GenotypeMultiDataType:
		case GenotypeControlMultiDataType:
			hdr.addAsciiColumn(PROBE_SET_NAME, info.getMaxProbeSetName());
			hdr.addUByteColumn(GENOTYPE_CALL_NAME);
			hdr.addFloatColumn(GENOTYPE_CONFIDENCE_NAME);
			break;

		case CopyNumberMultiDataType:
			hdr.addAsciiColumn(PROBE_SET_NAME, info.getMaxProbeSetName());
			hdr.addUByteColumn(COPY_NUMBER_CHR_NAME);
			hdr.addUIntColumn(COPY_NUMBER_POSITION_NAME);
			break;

		case CytoMultiDataType:
			hdr.addAsciiColumn(CYTO_REGION_NAME, info.getMaxProbeSetName());
			hdr.addUByteColumn(CYTO_CHR_NAME);
			hdr.addUIntColumn(CYTO_START_POSITION_NAME);
			hdr.addUIntColumn(CYTO_STOP_POSITION_NAME);
			hdr.addUByteColumn(CYTO_CALL_NAME);
			hdr.addFloatColumn(CYTO_CONFIDENCE_NAME);
			break;

		default:
			break;
		}
		Iterator<ColumnInfo> it = info.getMetricColumns().iterator();
		while (it.hasNext()) {
			hdr.addColumn(it.next());
		}
	}

	private DataSetHeader getDataSetHeader(MultiDataType dataType) {
		DataGroupHeader dcHdr = genericData.getHeader().getDataGroup(0);
		int n = dcHdr.getDataSetCnt();
		for (int i = 0; i < n; i++) {
			DataSetHeader dpHdr = dcHdr.getDataSet(i);
			if (dpHdr.getName().equals(DATA_SET_NAMES[dataType.ordinal()])) {
				return dpHdr;
			}
		}
		return null;
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
		return getParams(AffymetrixParameterConsts.ALGORITHM_PARAM_NAME_PREFIX);
	}

	public void addSummaryParams(List<ParameterNameValue> params) {
		super.addParams(AffymetrixParameterConsts.CHIP_SUMMARY_PARAM_NAME_PREFIX, params);
	}

	public void addAlgParams(List<ParameterNameValue> params) {
		super.addParams(AffymetrixParameterConsts.ALGORITHM_PARAM_NAME_PREFIX, params);
	}

	/**
	 * gets the summary parameters
	 * 
	 * @return The summary parameters.
	 */
	public List<ParameterNameValue> getSummaryParams() {
		return getParams(AffymetrixParameterConsts.CHIP_SUMMARY_PARAM_NAME_PREFIX);
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
	 * gets the probe set data.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @return The genotype results.
	 */
	public ProbeSetMultiDataGenotypeData getGenotypeEntry(MultiDataType dataType, int index) throws IOException,
			UnsignedOutOfLimitsException {
		return getGenericGenotypeEntry(dataType, index);
	}

	/**
	 * gets the probe set data.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @return The results.
	 */
	public ProbeSetMultiDataExpressionData getExpressionEntry(MultiDataType dataType, int index) throws IOException,
			UnsignedOutOfLimitsException {
		return getGenericExpressionEntry(dataType, index);
	}

	/**
	 * Gets the copy number data.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @return The results.
	 */
	public ProbeSetMultiDataCopyNumberData getCopyNumberEntry(MultiDataType dataType, int index) throws IOException,
			UnsignedOutOfLimitsException {
		return getGenericCopyNumberEntry(dataType, index);
	}

	public ProbeSetMultiDataCytoRegionData getCytoEntry(MultiDataType dataType, int index) throws IOException,
			UnsignedOutOfLimitsException {
		return getGenericCytoRegionEntry(dataType, index);
	}

	private ProbeSetMultiDataCytoRegionData getGenericCytoRegionEntry(MultiDataType dataType, int index)
			throws IOException, UnsignedOutOfLimitsException {
		DataSetInfo ds = openMultiDataDataSet(dataType);
		ProbeSetMultiDataCytoRegionData entry = null;
		if ((ds != null) && (ds.getEntries() != null) && ds.getEntries().isOpen()) {
			int colIndex = 0;
			entry = new ProbeSetMultiDataCytoRegionData();
			entry.setName(ds.getEntries().getDataString8(index, colIndex++));
			entry.setChr(ds.getEntries().getDataByte(index, colIndex++));
			entry.setStartPosition(ds.getEntries().getDataUInt(index, colIndex++));
			entry.setStopPosition(ds.getEntries().getDataUInt(index, colIndex++));
			entry.setCall(ds.getEntries().getDataByte(index, colIndex++));
			entry.setConfidence(ds.getEntries().getDataFloat(index, colIndex++));
			List<ParameterNameValue> extraMetrics = getExtraMetricEntries(ds, index, colIndex);
			if (extraMetrics != null) {
				entry.addMetrics(extraMetrics);
			}
		}
		return entry;
	}

	/**
	 * Gets the probe set data.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @return The genotype results.
	 */
	private ProbeSetMultiDataGenotypeData getGenericGenotypeEntry(MultiDataType dataType, int index) throws IOException,
			UnsignedOutOfLimitsException {
		DataSetInfo ds = openMultiDataDataSet(dataType);
		ProbeSetMultiDataGenotypeData entry = null;
		if ((ds != null) && (ds.getEntries() != null) && ds.getEntries().isOpen()) {
			int colIndex = 0;
			entry = new ProbeSetMultiDataGenotypeData();
			entry.setName(ds.getEntries().getDataString8(index, colIndex++));
			entry.setCall(ds.getEntries().getDataByte(index, colIndex++));
			entry.setConfidence(ds.getEntries().getDataFloat(index, colIndex++));
			List<ParameterNameValue> extraMetrics = getExtraMetricEntries(ds, index, colIndex);
			if (extraMetrics != null) {
				entry.addMetrics(extraMetrics);
			}
		}
		return entry;
	}

	/**
	 * Gets the probe set data.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @param entry
	 *          The expression results.
	 */
	private ProbeSetMultiDataExpressionData getGenericExpressionEntry(MultiDataType dataType, int index)
			throws IOException, UnsignedOutOfLimitsException {
		DataSetInfo ds = openMultiDataDataSet(dataType);
		ProbeSetMultiDataExpressionData entry = null;
		if ((ds != null) && (ds.getEntries() != null) && ds.getEntries().isOpen()) {
			int colIndex = 0;
			entry = new ProbeSetMultiDataExpressionData();
			entry.setName(ds.getEntries().getDataString8(index, colIndex++));
			entry.setQuantification(ds.getEntries().getDataFloat(index, colIndex++));
			List<ParameterNameValue> metrics = getExtraMetricEntries(ds, index, colIndex);
			if (metrics != null) {
				entry.addMetrics(metrics);
			}
		}
		return entry;
	}

	private ProbeSetMultiDataCopyNumberData getGenericCopyNumberEntry(MultiDataType dataType, int index)
			throws IOException, UnsignedOutOfLimitsException {
		DataSetInfo ds = openMultiDataDataSet(dataType);
		ProbeSetMultiDataCopyNumberData entry = null;
		if ((ds != null) && (ds.getEntries() != null) && ds.getEntries().isOpen()) {
			int colIndex = 0;
			entry = new ProbeSetMultiDataCopyNumberData();
			entry.setName(ds.getEntries().getDataString8(index, colIndex++));
			entry.setChr(ds.getEntries().getDataByte(index, colIndex++));
			entry.setPosition(ds.getEntries().getDataInt(index, colIndex++));
			List<ParameterNameValue> extraMetrics = getExtraMetricEntries(ds, index, colIndex);
			if (extraMetrics != null) {
				entry.addMetrics(extraMetrics);
			}
		}
		return entry;
	}

	/**
	 * Get the extra metric columns.
	 * 
	 * @param ds
	 *          The data set info.
	 * @param rowIndex
	 *          The row index.
	 * @param colIndex
	 *          The column index
	 * @param metrics
	 *          The results.
	 */
	private List<ParameterNameValue> getExtraMetricEntries(DataSetInfo ds, int rowIndex, int colIndex)
			throws IOException, UnsignedOutOfLimitsException {
		int ncols = ds.getMetricColumns().size();
		if (ncols == 0) {
			return null;
		}
		List<ParameterNameValue> metrics = new ArrayList<ParameterNameValue>();
		for (int icol = 0; icol < ncols; icol++) {
			ParameterNameValue nv = new ParameterNameValue();
			ColumnInfo cinfo = ds.getMetricColumns().get(icol);
			nv.setName(cinfo.getName());
			switch (cinfo.getColumnType()) {
			case ByteColType: {
				byte val = ds.getEntries().getDataByte(rowIndex, colIndex++);
				nv.setValueInt8(val);
			}
				break;

			case UByteColType: {
				UByte val = new UByte((short)(0xFF & (int)ds.getEntries().getDataByte(rowIndex, colIndex++)));
				nv.setValueUInt8(val);
			}
				break;

			case ShortColType: {
				short val = ds.getEntries().getDataShort(rowIndex, colIndex++);
				nv.setValueInt16(val);
			}
				break;

			case UShortColType: {
				UShort val = new UShort(0xFFFF & (int)ds.getEntries().getDataShort(rowIndex, colIndex++));
				nv.setValueUInt16(val);
			}
				break;

			case IntColType: {
				int val = ds.getEntries().getDataInt(rowIndex, colIndex++);
				nv.setValueInt32(val);
			}
				break;

			case UIntColType: {
				UInt val = new UInt(0xFFFFFFFFL & (long)ds.getEntries().getDataInt(rowIndex, colIndex++));
				nv.setValueUInt32(val);
			}
				break;

			case FloatColType: {
				float val = ds.getEntries().getDataFloat(rowIndex, colIndex++);
				nv.setValueFloat(val);
			}
				break;

			case ASCIICharColType: {
				String val = ds.getEntries().getDataString8(rowIndex, colIndex++);
				nv.setValueAscii(val);
			}
				break;

			case UnicodeCharColType: {
				String val = ds.getEntries().getDataString16(rowIndex, colIndex++);
				nv.setValueText(val);
			}
				break;
			}
			metrics.add(nv);
		}
		return metrics;
	}

	/**
	 * gets the call of the probe set.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @return The call.
	 */
	public byte getGenoCall(MultiDataType dataType, int index) throws UnsignedOutOfLimitsException {
		byte call = 0;
		DataSetInfo ds = openMultiDataDataSet(dataType);
		if ((ds != null) && (ds.getEntries() != null) && ds.getEntries().isOpen()) {
			call = ds.getEntries().getDataByte(index, 1);
		}
		return call;
	}

	/**
	 * gets the confidence in the call of the probe set.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @return The confidence.
	 */
	public float getGenoConfidence(MultiDataType dataType, int index) throws UnsignedOutOfLimitsException {
		float conf = 0.0f;
		DataSetInfo ds = openMultiDataDataSet(dataType);
		if ((ds != null) && (ds.getEntries() != null) && ds.getEntries().isOpen()) {
			conf = ds.getEntries().getDataFloat(index, 2);
		}
		return conf;
	}

	/**
	 * Gets the quantification of the probe set.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @return The quantification.
	 */
	public float getExpressionQuantification(MultiDataType dataType, int index) throws UnsignedOutOfLimitsException {
		float quant = 0.0f;
		DataSetInfo ds = openMultiDataDataSet(dataType);
		if ((ds != null) && (ds.getEntries() != null) && ds.getEntries().isOpen()) {
			quant = ds.getEntries().getDataFloat(index, 1);
		}
		return quant;
	}

	/**
	 * get the probe set name.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @return The probe set name.
	 */
	public String getProbeSetName(MultiDataType dataType, int index) throws IOException, UnsignedOutOfLimitsException {
		String name = null;
		DataSetInfo ds = openMultiDataDataSet(dataType);
		if ((ds != null) && (ds.getEntries() != null) && ds.getEntries().isOpen()) {
			name = ds.getEntries().getDataString8(index, 0);
		}
		return name;
	}

	/**
	 * Opens a group for reading.
	 * 
	 * @param dataType
	 *          The data type
	 */
	public DataSetInfo openMultiDataDataSet(MultiDataType dataType) {
		if (dataSetInfo.containsKey(dataType)) {
			return dataSetInfo.get(dataType);
		}
		else {
			try {
				DataSetInfo ds = new DataSetInfo();
				ds.setEntries(genericData.getDataSet(MULTI_DATA_NAME, DATA_SET_NAMES[dataType.ordinal()]));
				if (ds.getEntries() != null) {
					ds.getEntries().open();
					int ncols = ds.getEntries().getHeader().getColumnCnt();
					ds.clearMetricColumns();
					int startCol = 0;
					if ((dataType == MultiDataType.ExpressionMultiDataType)
							|| (dataType == MultiDataType.ExpressionControlMultiDataType)) {
						startCol = 2;
					}
					else if ((dataType == MultiDataType.GenotypeMultiDataType)
							|| (dataType == MultiDataType.GenotypeControlMultiDataType)) {
						startCol = 3;
					}
					else if (dataType == MultiDataType.CopyNumberMultiDataType) {
						startCol = 3;
					}
					else if (dataType == MultiDataType.CytoMultiDataType) {
						startCol = 6;
					}

					for (int icol = startCol; icol < ncols; icol++) {
						ds.addMetricColumn((ds.getEntries().getHeader().getColumnInfo(icol)));
					}
					dataSetInfo.put(dataType, ds);
					return ds;
				}
			}
			catch (Throwable t) {
			}
		}
		return null;
	}

	/**
	 * get the length of the metric columns.
	 * 
	 * @param col
	 *          The column index (of the metric columns)
	 * @return The length.
	 */
	public int getMetricColumnLength(MultiDataType dataType, int col) {
		DataSetInfo ds = openMultiDataDataSet(dataType);
		if ((ds != null) && (ds.getEntries() != null) && ds.getEntries().isOpen()) {
			ColumnInfo info = ds.getMetricColumns().get(col);
			return info.getLength();
		}
		return 0;
	}

	/**
	 * Get the length of the metric columns.
	 * 
	 * @param dataType
	 *          The data type
	 * @return The number of columns.
	 */
	public int getNumMetricColumns(MultiDataType dataType) {
		DataSetInfo ds = openMultiDataDataSet(dataType);
		if ((ds != null) && (ds.getEntries() != null) && ds.getEntries().isOpen()) {
			return ds.getMetricColumns().size();
		}
		return 0;
	}

	/**
	 * Get the metric column name.
	 * 
	 * @param dataType
	 *          The data type
	 * @param colIndex
	 *          the metric column index
	 * @return The column name
	 */
	public String getMetricColumnName(MultiDataType dataType, int colIndex) {
		DataSetInfo ds = openMultiDataDataSet(dataType);
		if ((ds != null) && (ds.getEntries() != null) && ds.getEntries().isOpen()) {
			ColumnInfo col = ds.getMetricColumns().get(colIndex);
			return col.getName();
		}
		return null;
	}
};
