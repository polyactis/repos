////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
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

package affymetrix.fusion.chp;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import affymetrix.calvin.data.CHPMultiDataData;
import affymetrix.calvin.data.DataSetInfo;
import affymetrix.calvin.data.GenericData;
import affymetrix.calvin.data.ProbeSetMultiDataCopyNumberData;
import affymetrix.calvin.data.ProbeSetMultiDataExpressionData;
import affymetrix.calvin.data.ProbeSetMultiDataGenotypeData;
import affymetrix.calvin.data.CHPMultiDataData.MultiDataType;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.parsers.CHPMultiDataFileReader;
import affymetrix.calvin.utils.AffymetrixGuidType;

/** This class provides storage and reading capabilities for multi-data CHP files */
public class FusionCHPMultiDataData extends FusionCHPData {

	/** Creates a new instance of FusionCHPMultiDataData */
	public FusionCHPMultiDataData() {
		chpData = null;
	}

	/** The chip object. */
	private CHPMultiDataData chpData;

	/**
	 * Get the id of the file (only valid for Command Console "calvin" files)
	 * 
	 * @return The unique file id.
	 */
	@Override
	public AffymetrixGuidType getFileId() {
		return chpData.getGenericData().getFileIdentifier();
	}

	/** Returns the GenericData object associated with a Calvin file, NULL for GCOS files. */
	public GenericData getGenericData() {
		return chpData.getGenericData();
	}

	/** The data set information */
	public Map<MultiDataType, DataSetInfo> getDataSetInfo() {
		return chpData.getDataSetInfo();
	}

	/**
	 * The maximum length of a probe set name.
	 * 
	 * @param dataType
	 *          The data type
	 * @return The maximum probe set name length
	 */
	public int getMaxProbeSetName(MultiDataType dataType) {
		return chpData.getMaxProbeSetName(dataType);
	}

	/** Clears the members. */
	public void clear() {
		chpData.clear();
	}

	/**
	 * Gets the name of the algorithm.
	 * 
	 * @return The algorithm name.
	 */
	public String getAlgName() {
		return chpData.getAlgName();
	}

	/**
	 * Gets the algorithm version.
	 * 
	 * @return The version.
	 */
	public String getAlgVersion() {
		return chpData.getAlgVersion();
	}

	/** Sets the array type */
	public String getArrayType() {
		return chpData.getArrayType();
	}

	/**
	 * Gets the algorithm parameters
	 * 
	 * @return The algoirhtm parameters.
	 */
	public List<ParameterNameValue> getAlgParams() {
		return chpData.getAlgParams();
	}

	/**
	 * Gets the summary parameters
	 * 
	 * @return The summary parameters.
	 */
	public List<ParameterNameValue> getSummaryParams() {
		return chpData.getSummaryParams();
	}

	/**
	 * Gets the number of entries (probe sets)
	 * 
	 * @param dataType
	 *          The data type
	 */
	public int getEntryCount(MultiDataType dataType) {
		return chpData.getEntryCount(dataType);
	}

	/**
	 * Gets the data for the given row.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @return The entry.
	 */
	public ProbeSetMultiDataGenotypeData getGenotypeEntry(MultiDataType dataType, int index) throws IOException,
			UnsignedOutOfLimitsException {
		return chpData.getGenotypeEntry(dataType, index);
	}

	/**
	 * Gets the probe set data.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @param entry
	 *          The copy number results.
	 */
	public ProbeSetMultiDataCopyNumberData getCopyNumberEntry(MultiDataType dataType, int index) throws IOException,
			UnsignedOutOfLimitsException {
		return chpData.getCopyNumberEntry(dataType, index);
	}

	/**
	 * Gets the probe set data.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @return The expression results.
	 */
	public ProbeSetMultiDataExpressionData getExpressionEntry(MultiDataType dataType, int index) throws IOException,
			UnsignedOutOfLimitsException {
		return chpData.getExpressionEntry(dataType, index);
	}

	/**
	 * Gets the call of the probe set.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @return The call.
	 */
	public byte getGenoCall(MultiDataType dataType, int index) throws IOException, UnsignedOutOfLimitsException {
		return chpData.getGenoCall(dataType, index);
	}

	/**
	 * Gets the confidence in the call of the probe set.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @return The confidence.
	 */
	public float getGenoConfidence(MultiDataType dataType, int index) throws IOException, UnsignedOutOfLimitsException {
		return chpData.getGenoConfidence(dataType, index);
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
	public float getExpressionQuantification(MultiDataType dataType, int index) throws IOException,
			UnsignedOutOfLimitsException {
		return chpData.getExpressionQuantification(dataType, index);
	}

	/**
	 * Get the probe set name.
	 * 
	 * @param dataType
	 *          The data type
	 * @param index
	 *          The row index.
	 * @return The probe set name.
	 */
	public String getProbeSetName(MultiDataType dataType, int index) throws IOException, UnsignedOutOfLimitsException {
		return chpData.getProbeSetName(dataType, index);
	}

	/**
	 * Get the length of the metric columns.
	 * 
	 * @param dataType
	 *          The data type
	 * @param col
	 *          The column index (of the metric columns)
	 * @return The length.
	 */
	public int getMetricColumnLength(MultiDataType dataType, int col) {
		return chpData.getMetricColumnLength(dataType, col);
	}

	/**
	 * Get the length of the metric columns.
	 * 
	 * @param dataType
	 *          The data type
	 * @return The number of columns.
	 */
	public int getNumMetricColumns(MultiDataType dataType) {
		return chpData.getNumMetricColumns(dataType);
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
		return chpData.getMetricColumnName(dataType, colIndex);
	}

	/**
	 * Reads the CHP file.
	 * 
	 * @return True if successful.
	 */
	@Override
	protected boolean read() {
		CHPMultiDataFileReader reader = new CHPMultiDataFileReader();
		chpData = new CHPMultiDataData();
		reader.setFilename(filename);
		try {
			reader.read(chpData);
			return true;
		}
		catch (Throwable t) {
			return false;
		}
	}

	/**
	 * Reads the header of the CHP file
	 * 
	 * @return True if successful
	 */
	@Override
	protected boolean readHeader() {
		return read();
	}

	/** A class to register the multi-data CHP reader. */
	private static class Reg extends FusionCHPDataReg {
		/** Constructor - register the multi-data file type. */
		public Reg() {
			super();
			List<AffymetrixGuidType> ids = new ArrayList<AffymetrixGuidType>();
			AffymetrixGuidType guid = new AffymetrixGuidType();
			guid.setGuid(CHPMultiDataData.CHP_MULTI_DATA_TYPE);
			ids.add(guid);
			setFileTypeIds(ids);
		}

		/**
		 * Creates a multi-data CHP object.
		 * 
		 * @return The multi-data CHP object.
		 */
		@Override
		public FusionCHPData makeObject() {
			return new FusionCHPMultiDataData();
		}
	};

	/** The one and only registration object. This registers the class as a CHP reader. */
	private static Reg reg = null;

	/** Register the reader with the system. */
	public static void registerReader() {
		if (FusionCHPMultiDataData.reg == null) {
			FusionCHPMultiDataData.reg = new Reg();
		}
	}

	/**
	 * Converts the type to the multi-data CHP type.
	 * 
	 * @param chip
	 *          The pointer to the CHP data object.
	 * @return The multi-data CHP data type or NULL if not compatible.
	 */
	public static FusionCHPMultiDataData fromBase(FusionCHPData chip) {
		if (chip == null) {
			return null;
		}
		String chipName = chip.getClass().getName();
		String genName = FusionCHPMultiDataData.class.getName();
		if (chipName.compareTo(genName) == 0) {
			return (FusionCHPMultiDataData)chip;
		}
		return null;
	}

}
