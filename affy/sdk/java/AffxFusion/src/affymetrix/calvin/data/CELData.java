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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import affymetrix.calvin.data.ColumnInfo.DataSetColumnTypes;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.AffymetrixParameterConsts;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.portability.UByte;

/** This is the container class for CEL data. */
public class CELData {

	/** Data Group name */
	public static final String CelDataGroupName = "Default Group";

	/** Name of the intensity data set */
	public static final String CelIntensityLabel = "Intensity";

	/** Name of the standard deviation data set */
	public static final String CelStdDevLabel = "StdDev";

	/** Name of the pixel data set */
	public static final String CelPixelLabel = "Pixel";

	/** Name of the outlier data set */
	public static final String CelOutlierLabel = "Outlier";

	/** Name of the mask data set */
	public static final String CelMaskLabel = "Mask";

	/** Cel file version */
	public static final byte CurrentCelFileVersion = 1;

	/** The row parameter label. */
	private static String CelRowsLabel = "affymetrix-cel-rows";

	/** The col parameter label. */
	private static String CelColsLabel = "affymetrix-cel-cols";

	/** The attribute name of the cel file version number */
	private static String CelFileVersionNumberName = "affymetrix-file-version";

	/** The guid for the CEL file. */
	private static String INTENSITY_DATA_TYPE = "affymetrix-calvin-intensity";

	/** A set of outlier cell coordinates. */
	private Map<Integer, Boolean> outliers = new HashMap<Integer, Boolean>();

	/** Generic layer object */
	private GenericData genericData = new GenericData();

	private DataSetColumnTypes intensityColumnType;

	// DataSet cache - initialized on first use and Delete in destructor
	/** Intensity DataSet */
	private DataSet dpInten = null;

	/** Stdev DataSet */
	private DataSet dpStdev = null;

	/** NumPixels DataSet */
	private DataSet dpPixels = null;

	/** Indicates whether an attempt to read the outlier data set has been made. */
	private boolean outlierPlaneRead = false;

	/** Indicates whether an attempt to read the mask data set has been made. */
	private boolean maskPlaneRead = false;

	/** A set of masked cell coordinates. */
	private Map<Integer, Boolean> masked = new HashMap<Integer, Boolean>();

	/** keep rows from being read from the header all the time */
	private int cachedRows = -1;

	/** keep cols from being read from the header all the time */
	private int cachedCols = -1;

	/** Default constructor */
	public CELData() {
	}

	/**
	 * Constructor
	 * 
	 * @param filename
	 *          Name of the cel file.
	 */
	public CELData(String filename) {
		setFilename(filename);
		genericData.getHeader().getGenericDataHdr().setFileTypeId(INTENSITY_DATA_TYPE);
		DataGroupHeader dcHdr = new DataGroupHeader(CelDataGroupName);
		genericData.getHeader().addDataGroupHdr(dcHdr);
	}

	/**
	 * Clear the object members
	 */
	public void clear() {
		genericData.getHeader().clear();
		dpInten = null;
		dpStdev = null;
		dpPixels = null;
		outlierPlaneRead = false;
		outliers.clear();
		maskPlaneRead = false;
		masked.clear();
		cachedRows = -1;
		cachedCols = -1;
	}

	/**
	 * set the file name
	 * 
	 * @param p
	 *          file name
	 */
	public void setFilename(String p) {
		genericData.getHeader().setFilename(p);
	}

	/**
	 * get the file name
	 * 
	 * @return file name
	 */
	public String getFilename() {
		return genericData.getHeader().getFilename();
	}

	public FileHeader getFileHeader() {
		return genericData.getHeader();
	}

	/**
	 * get the version of the CEL file
	 * 
	 * @return CEL file version
	 */
	public UByte getVersion() throws UnsignedOutOfLimitsException {
		GenericDataHeader hdr = genericData.getHeader().getGenericDataHdr();
		ParameterNameValue paramType = hdr.findNameValParam(CelFileVersionNumberName);
		if (paramType != null) {
			return paramType.getValueUInt8();
		}
		return UByte.ZERO;
	}

	/**
	 * get the array type
	 * 
	 * @return array type name
	 */
	public String getArrayType() {
		return getWStringFromGenericHdr(AffymetrixParameterConsts.ARRAY_TYPE_PARAM_NAME);
	}

	/**
	 * Get the master file
	 * 
	 * @return The name of the master file.
	 */
	public String getMasterFileName() {
		return getWStringFromGenericHdr(AffymetrixParameterConsts.MASTER_FILE_PARAM_NAME);
	}

	/**
	 * Get the library package
	 * 
	 * @return The name of the library package.
	 */
	public String getLibraryPackageName() {
		return getWStringFromGenericHdr(AffymetrixParameterConsts.LIBRARY_PACKAGE_PARAM_NAME);
	}

	/**
	 * get the name of the algorithm used to generate the results.
	 * 
	 * @return Algorithm name
	 */
	public String getAlgorithmName() {
		return getWStringFromGenericHdr(AffymetrixParameterConsts.ALGORITHM_NAME_PARAM_NAME);
	}

	/**
	 * get the number of rows of cells on the array.
	 * 
	 * @return Number of rows of cells.
	 */
	public int getRows() {
		if (cachedRows == -1) {
			cachedRows = getInt32FromGenericHdrParameterList(CelRowsLabel);
		}
		return cachedRows;
	}

	/**
	 * get the number of columns of cells on the array.
	 * 
	 * @return Number of columns of cells.
	 */
	public int getCols() {
		if (cachedCols == -1) {
			cachedCols = getInt32FromGenericHdrParameterList(CelColsLabel);
		}
		return cachedCols;
	}

	/**
	 * Return the number of cells on the array. This is the number of the intensity data elements and is == getRows() *
	 * getCols() Stdev, NumPixels, Outlier and Masked data is optional, but if present they will have getNumCells
	 * elements.
	 * 
	 * @return Number of cells in the array.
	 */
	public int getNumCells() {
		int rows = 0;
		try {
			prepareIntensityPlane();
			if (dpInten != null) {
				rows = dpInten.getRows();
			}
		}
		catch (Throwable t) {
		}
		return rows;
	}

	/**
	 * Return the algorithm parameters. The algorithm parameter prefix is removed from the name.
	 * 
	 * @return Vector with algorithm parameters.
	 */
	public List<ParameterNameValue> getAlgorithmParameters() {
		List<ParameterNameValue> algParams = new ArrayList<ParameterNameValue>();
		List<ParameterNameValue> allParams = genericData.getHeader().getGenericDataHdr().getNameValParams();
		for (int i = 0; i < allParams.size(); i++) {
			ParameterNameValue param = allParams.get(i);
			String name = param.getName();
			if (name.startsWith(AffymetrixParameterConsts.ALGORITHM_PARAM_NAME_PREFIX)) {
				ParameterNameValue algParam = new ParameterNameValue(param);
				algParam.setName(name.substring(AffymetrixParameterConsts.ALGORITHM_PARAM_NAME_PREFIX.length(), name.length()));
				algParams.add(algParam);
			}
		}
		return algParams;
	}

	/**
	 * Return an algorithm parameter given a name.
	 * 
	 * @param name
	 *          Name of the parameter to find.
	 * @return The found parameter.
	 */
	public ParameterNameValue findAlgorithmParameter(String name) {
		String paramName = AffymetrixParameterConsts.ALGORITHM_PARAM_NAME_PREFIX + name;
		GenericDataHeader hdr = genericData.getHeader().getGenericDataHdr();
		ParameterNameValue param = hdr.findNameValParam(paramName);
		if (param != null) {
			ParameterNameValue foundParam = new ParameterNameValue(param);
			int len = AffymetrixParameterConsts.ALGORITHM_PARAM_NAME_PREFIX.length();
			foundParam.setName(param.getName().substring(len));
			return foundParam;
		}
		return null;
	}

	/**
	 * get the cell intensity.
	 * 
	 * @param cellIdx
	 *          Index of the cell.
	 * @return The intensity of the cell.
	 */
	public float getIntensity(int cellIdx) throws UnsignedOutOfLimitsException {
		prepareIntensityPlane();
		if (dpInten != null) {
			if (intensityColumnType == DataSetColumnTypes.FloatColType) {
				return dpInten.getDataFloat(cellIdx, 0);
			}
			else {
				return dpInten.getDataShort(cellIdx, 0);
			}
		}
		return 0.0f;
	}

	/**
	 * get the cell stdv.
	 * 
	 * @param cellIdx
	 *          Index of the cell.
	 * @return The stdv of the cell.
	 */
	public float getStdv(int cellIdx) throws UnsignedOutOfLimitsException {
		prepareStdevPlane();
		if (dpStdev != null) {
			return dpStdev.getDataFloat(cellIdx, 0);
		}
		return 0.0f;
	}

	/**
	 * get the cell stdv.
	 * 
	 * @param cellIdx
	 *          Index of the cell.
	 * @return The stdv of the cell.
	 */
	public short getPixels(int cellIdx) throws UnsignedOutOfLimitsException {
		prepareNumPixelPlane();
		if (dpPixels != null) {
			return dpPixels.getDataShort(cellIdx, 0);
		}
		return 0;
	}

	/**
	 * Indicates whether there are standard deviation values.
	 * 
	 * @return True if there are standard deviation values.
	 */
	public boolean hasStdev() {
		DataSetHeader dph = findDataSetHeader(CelStdDevLabel);
		if (dph != null) {
			return (dph.getRowCnt() > 0);
		}
		return false;
	}

	/**
	 * Indicates whether there are number of pixel values.
	 * 
	 * @return True if there are number of pixels values.
	 */
	public boolean hasNumPixels() {
		DataSetHeader dph = findDataSetHeader(CelPixelLabel);
		if (dph != null) {
			return (dph.getRowCnt() > 0);
		}
		return false;
	}

	/**
	 * Return a reference to the generic layer object
	 * 
	 * @return Generic layer object
	 */
	public GenericData getGenericData() {
		return genericData;
	} // should be a friend method only

	/** Prepare to read intensity data */
	private void prepareIntensityPlane() {
		if (dpInten == null) {
			try {
				dpInten = genericData.getDataSet(CelDataGroupName, CelIntensityLabel);
				if (dpInten != null) {
					dpInten.open();
					intensityColumnType = dpInten.getHeader().getColumnInfo(0).getColumnType();
				}
			}
			catch (Throwable t) {
			}
		}
	}

	/** Prepare to read the standard deviation data */
	private void prepareStdevPlane() {
		if (dpStdev == null) {
			try {
				dpStdev = genericData.getDataSet(CelDataGroupName, CelStdDevLabel);
				if (dpStdev != null) {
					dpStdev.open();
				}
			}
			catch (Throwable t) {
			}
		}
	}

	/** Prepare to read the number of pixel data */
	private void prepareNumPixelPlane() {
		if (dpPixels == null) {
			try {
				dpPixels = genericData.getDataSet(CelDataGroupName, CelPixelLabel);
				if (dpPixels != null) {
					dpPixels.open();
				}
			}
			catch (Throwable t) {
			}
		}
	}

	/** Prepare to read the outlier data */
	private void prepareOutlierPlane() {
		if (outlierPlaneRead) {
			return;
		}
		outlierPlaneRead = true; // Read attempted

		try {
			DataSet dpOutlier = genericData.getDataSet(CelDataGroupName, CelOutlierLabel);
			if (dpOutlier != null) {
				if (dpOutlier.open()) {
					outliers.clear();
					int rows = dpOutlier.getRows();
					for (int row = 0; row < rows; ++row) {
						short x = dpOutlier.getDataShort(row, 0);
						short y = dpOutlier.getDataShort(row, 1);
						outliers.put(new Integer(xyToIndex(x, y)), new Boolean(true));
					}
				}
				dpOutlier.delete();
				dpOutlier = null;
			}
		}
		catch (Throwable t) {
		}
	}

	/** Prepare to read the mask data */
	private void prepareMaskedPlane() {
		if (maskPlaneRead) {
			return;
		}

		maskPlaneRead = true; // Read attempted

		try {
			DataSet dpMasked = genericData.getDataSet(CelDataGroupName, CelMaskLabel);
			if (dpMasked != null) {
				if (dpMasked.open()) {
					masked.clear();
					int rows = dpMasked.getRows();
					for (int row = 0; row < rows; ++row) {
						short x = dpMasked.getDataShort(row, 0);
						short y = dpMasked.getDataShort(row, 1);
						masked.put(new Integer(xyToIndex(x, y)), new Boolean(true));
					}
				}
				dpMasked = null;
			}
		}
		catch (Throwable t) {
		}
	}

	/** Prepare to read all data */
	// private void prepareAllPlanes() {
	// prepareIntensityPlane();
	// prepareStdevPlane();
	// prepareNumPixelPlane();
	// prepareOutlierPlane();
	// prepareMaskedPlane();
	// }
	/**
	 * Read an int value from the GenericDataHeader parameter list.
	 * 
	 * @param name
	 *          Name of the parameter to read
	 * @return int value of the named parameter
	 */
	private int getInt32FromGenericHdrParameterList(String name) {
		int result = 0;
		GenericDataHeader hdr = genericData.getHeader().getGenericDataHdr();
		ParameterNameValue paramType = hdr.findNameValParam(name);
		if (paramType != null) {
			result = paramType.getValueInt32();
		}
		return result;
	}

	/**
	 * Read an wstring value from the GenericDataHeader parameter list.
	 * 
	 * @param name
	 *          Name of the parameter to read
	 * @return wstring value of the named parameter
	 */
	private String getWStringFromGenericHdr(String name) {
		String result = null;
		GenericDataHeader hdr = genericData.getHeader().getGenericDataHdr();
		ParameterNameValue paramType = hdr.findNameValParam(name);
		if (paramType != null) {
			result = paramType.getValueText();
		}
		return result;
	}

	/** Get the number of outliers */
	public int getNumOutliers() {
		prepareOutlierPlane();
		return outliers.size();
	}

	/**
	 * Check if the cell is an outlier (outlier flag is true)
	 * 
	 * @param cellIdx
	 *          Cell index
	 * @return True if the cell outlier flag is true
	 */
	public boolean isOutlier(int cellIdx) {
		prepareOutlierPlane();
		Boolean outlier = outliers.get(new Integer(cellIdx));
		if (outlier != null) {
			return outlier.booleanValue();
		}
		return false;
	}

	/** Get the number of outliers */
	public int getNumMasked() {
		prepareMaskedPlane();
		return masked.size();
	}

	/**
	 * Check if the cell is masked (mask flag is true)
	 * 
	 * @param cellIdx
	 *          cell index
	 * @return True if the cell mask flag is true
	 */
	public boolean isMasked(int cellIdx) {
		prepareMaskedPlane();
		Boolean m = masked.get(new Integer(cellIdx));
		if (m != null) {
			return m.booleanValue();
		}
		return false;
	}

	/**
	 * Gets the X coordinates from index.
	 * 
	 * @param index
	 *          The 0 based index to the entry array.
	 * @return X coordinate
	 */
	public int indexToX(int index) {
		return index % getCols();
	}

	/**
	 * Gets the Y coordinates from index.
	 * 
	 * @param index
	 *          The 0 based index to the entry array.
	 * @return Y coordinate
	 */
	public int indexToY(int index) {
		return index / getCols();
	}

	/**
	 * Maps X/Y coordinates to CEL file index.
	 * 
	 * @param x
	 *          The x coordinate
	 * @param y
	 *          The y coordinate.
	 * @return The index to the entry array.
	 */
	public int xyToIndex(int x, int y) {
		return xyToIndex(x, y, getRows(), getCols());
	}

	/**
	 * Maps X/Y coordinates to CEL file index.
	 * 
	 * @param x
	 *          The x coordinate.
	 * @param y
	 *          The y coordinate.
	 * @param r
	 *          The number of rows.
	 * @param c
	 *          The number of columns.
	 * @return The index to the intensity arrays.
	 */
	public static int xyToIndex(int x, int y, int r, int c) {
		return ((y * c) + x);
	}

	/**
	 * Determine the xy coordinate given a cell index.
	 * 
	 * @param cellIdx
	 *          Cell index
	 * @return Cell coordinate
	 */
	// private XYCoord computeXY(int cellIdx) {
	// XYCoord coord = new XYCoord();
	// coord.setY((short) (cellIdx / getCols()));
	// coord.setX((short) (cellIdx - getCols() * coord.getY()));
	// return coord;
	// }
	/**
	 * Find a DataSetHeader by name.
	 * 
	 * @param name
	 *          DataSetHeader name
	 * @return Pointer to a DataSetHeader with name parameter, otherwise 0
	 */
	private DataSetHeader findDataSetHeader(String name) {
		DataGroupHeader dch = genericData.findDataGroupHeader(CelDataGroupName);
		if (dch != null) {
			DataSetHeader dph = GenericData.findDataSetHeader(dch, name);
			if (dph != null) {
				return dph;
			}
		}
		return null;
	}

}
