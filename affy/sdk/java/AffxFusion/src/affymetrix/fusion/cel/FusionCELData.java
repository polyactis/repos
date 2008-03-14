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

package affymetrix.fusion.cel;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.data.CELData;
import affymetrix.calvin.data.GenericData;
import affymetrix.calvin.data.GenericDataHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.AffymetrixParameterConsts;
import affymetrix.calvin.parameter.CELAlgorithmParameterNames;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.parameter.ParameterNameValue.ParameterType;
import affymetrix.calvin.parsers.CELFileReader;
import affymetrix.calvin.utils.AffymetrixGuidType;
import affymetrix.calvin.utils.FGridCoords;
import affymetrix.calvin.utils.FPoint;
import affymetrix.calvin.utils.IOUtils;
import affymetrix.fusion.FusionTagValuePair;
import affymetrix.gcos.GridCoordinates;
import affymetrix.gcos.TagValuePair;
import affymetrix.gcos.cel.CELFileData;
import affymetrix.portability.UShort;

/** Provides an interface to read CEL files in either GCOS or Calvin format. */
public class FusionCELData {

	/** Flag to read all of the data in a CEL file. */
	public static final int CEL_ALL = 1;

	/** Flag to read only the data from a CEL file. */
	public static final int CEL_DATA = 2;

	/** Flag to read the outlier and data sections from a CEL file. */
	public static final int CEL_OUTLIER = 4;

	/** Flag to read the mask and data sections from a CEL file. */
	public static final int CEL_MASK = 8;

	/** The GCOS file object. */
	private CELFileData gcosFile;

	/** The Calvin file object. */
	private CELData calvinFile;

	/** Error string */
	private String strError;

	/**
	 * Gets the error.
	 * 
	 * @return The last error message.
	 */
	public String getError() {
		return strError;
	}

	/** The file name. */
	private String fileName;

	/**
	 * Gets the file name.
	 * 
	 * @return The file name.
	 */
	public String getFileName() {
		return fileName;
	}

	/**
	 * Sets the file name.
	 * 
	 * @param value
	 *          The name of the CEL file to read.
	 */
	public void setFileName(String value) {
		fileName = value;
	}

	/**
	 * Get the id of the file (only valid for Command Console "calvin" files)
	 * 
	 * @return The unique file id.
	 */
	public AffymetrixGuidType getFileId() {
		if (calvinFile != null) {
			return calvinFile.getFileHeader().getGenericDataHdr().getFileId();
		}
		return null;
	}

	public GenericData getGenericData() {
		if (calvinFile != null) {
			return calvinFile.getGenericData();
		}
		return null;
	}

	/** Number of columns in array */
	public int getCols() {
		if (gcosFile != null) {
			return gcosFile.getHeader().getCols();
		}
		else if (calvinFile != null) {
			return calvinFile.getCols();
		}
		return 0;
	}

	/** Number of rows in array */
	public int getRows() {
		if (gcosFile != null) {
			return gcosFile.getHeader().getRows();
		}
		else if (calvinFile != null) {
			return calvinFile.getRows();
		}
		return 0;
	}

	/** Number of cells in array */
	public int getCells() {
		if (gcosFile != null) {
			return gcosFile.getHeader().getCells();
		}
		else if (calvinFile != null) {
			return calvinFile.getNumCells();
		}
		return 0;
	}

	/** Header information concatenated in a string */
	public String getHeader() {
		if (gcosFile != null) {
			return gcosFile.getHeader().getHeader();
		}
		else if (calvinFile != null) {
			return "";
		}
		return null;
	}

	/** Algorithm name */
	public String getAlg() {
		if (gcosFile != null) {
			return gcosFile.getHeader().getAlg();
		}
		else if (calvinFile != null) {
			return calvinFile.getAlgorithmName();
		}
		return null;
	}

	/**
	 * Retrieve algorithm parameters
	 * 
	 * @return Algorithm parameters
	 */
	public String getParams() {
		if (gcosFile != null) {
			String strparam = IOUtils.EMPTY;
			int n = gcosFile.getHeader().getParameters().size();
			for (int i = 0; i < n; i++) {
				if (i > 0) {
					strparam += ";";
				}
				TagValuePair param = gcosFile.getHeader().getParameters().get(i);
				strparam += param.getTag();
				strparam += ":";
				if (param.getValue() != null) {
					strparam += param.getValue();
				}
			}
			return strparam;
		}
		else if (calvinFile != null) {
			String strparam = "";
			List<ParameterNameValue> params = calvinFile.getAlgorithmParameters();
			for (int i = 0; i < params.size(); i++) {
				if (i > 0) {
					strparam += ";";
				}
				ParameterNameValue param = params.get(i);
				strparam += param.getName();
				strparam += ":";
				if (param.toString() != null) {
					strparam += param.toString();
				}
			}
		}
		return null;
	}

	/**
	 * Retrieve algorithm parameter of specified tag
	 * 
	 * @param tag
	 *          Algorithm parameter tag
	 * @return Algorithm parameter value
	 */
	public String getAlgorithmParameter(String tag) {
		if (gcosFile != null) {
			int n = gcosFile.getHeader().getParameters().size();
			for (int i = 0; i < n; i++) {
				TagValuePair param = gcosFile.getHeader().getParameters().get(i);
				if (tag.compareTo(param.getTag()) == 0) {
					return param.getValue();
				}
			}
		}
		else if (calvinFile != null) {
			List<ParameterNameValue> params = calvinFile.getAlgorithmParameters();
			for (int i = 0; i < params.size(); i++) {
				ParameterNameValue param = params.get(i);
				if (tag.compareTo(param.getName()) == 0) {
					return param.toString();
				}
			}
		}
		return null;
	}

	/**
	 * Retrieves the algorithm parameter name (tag) for a given index position.
	 * 
	 * @param index
	 *          The zero based index to the parameter array (0 to the number of alg parameters - 1).
	 * @return The parameter name (tag).
	 */
	public String getAlgorithmParameterTag(int index) {
		if (gcosFile != null) {
			TagValuePair param = gcosFile.getHeader().getParameters().get(index);
			return param.getTag();
		}
		else if (calvinFile != null) {
			ParameterNameValue param = calvinFile.getAlgorithmParameters().get(index);
			return param.getName();
		}
		return null;
	}

	/**
	 * Retrieves the number of algorithm parameters.
	 * 
	 * @return The number of algorithm parameters.
	 */
	public int getNumberAlgorithmParameters() {
		if (gcosFile != null) {
			return gcosFile.getHeader().getParameters().size();
		}
		else if (calvinFile != null) {
			return calvinFile.getAlgorithmParameters().size();
		}
		return 0;
	}

	/**
	 * Retrieve algorithm parameters
	 * 
	 * @return Algorithm parameters
	 */
	public String getAlgorithmParameters() {
		return getParams();
	}

	/** Algorithm parameters as a vector of FusionTagValuePair objects */
	public List<FusionTagValuePair> getParameters() {
		if (gcosFile != null) {
			List<TagValuePair> gcosParams = gcosFile.getHeader().getParameters();
			int n = gcosParams.size();
			List<FusionTagValuePair> params = new ArrayList<FusionTagValuePair>();
			for (int i = 0; i < n; i++) {
				FusionTagValuePair param = new FusionTagValuePair();
				TagValuePair gcosParam = gcosParams.get(i);
				param.setTag(gcosParam.getTag());
				param.setValue(gcosParam.getValue());
				params.add(param);
			}
			return params;
		}
		else if (calvinFile != null) {
			List<ParameterNameValue> calvinParams = calvinFile.getAlgorithmParameters();
			int n = calvinParams.size();
			List<FusionTagValuePair> params = new ArrayList<FusionTagValuePair>();
			for (int i = 0; i < n; i++) {
				FusionTagValuePair param = new FusionTagValuePair();
				ParameterNameValue calvinParam = calvinParams.get(i);
				param.setTag(calvinParam.getName());
				param.setValue(calvinParam.toString());
				param.setDetailed(calvinParam);
				params.add(param);
			}
			return params;
		}
		return null;
	}

	/** Chip type of array */
	public String getChipType() {
		if (gcosFile != null) {
			return gcosFile.getHeader().getChipType();
		}
		else if (calvinFile != null) {
			return calvinFile.getArrayType();
		}
		return null;
	}

	/**
	 * Get the master file
	 * 
	 * @return The name of the master file (null if GCOS format CEL file)
	 */
	public String getMasterFileName() {
		if (calvinFile != null) {
			return calvinFile.getMasterFileName();
		}
		return null;
	}

	/**
	 * Get the library package
	 * 
	 * @return The name of the library package (null if GCOS format CEL file)
	 */
	public String getLibraryPackageName() {
		if (calvinFile != null) {
			return calvinFile.getLibraryPackageName();
		}
		return null;
	}

	/** DAT header string */
	public String getDatHeader() throws UnsignedOutOfLimitsException {
		if (gcosFile != null) {
			return gcosFile.getHeader().getDatHeader();
		}
		else if (calvinFile != null) {
			GenericDataHeader gdh = calvinFile.getFileHeader().getGenericDataHdr().findParent(
					AffymetrixParameterConsts.SCAN_ACQUISITION_DATA_TYPE);
			if (gdh != null) {
				// found the right header, now look for the parameter
				ParameterNameValue nvt = gdh.findNameValParam(AffymetrixParameterConsts.DAT_HEADER_PARAM_NAME);
				if (nvt != null) {
					if (nvt.getParameterType() == ParameterType.TextType) {
						return nvt.getValueText();
					}
				}
				nvt = gdh.findNameValParam(AffymetrixParameterConsts.PARTIAL_DAT_HEADER_PARAM_NAME);
				if (nvt != null) {
					if (nvt.getParameterType() == ParameterType.TextType) {
						String partialDatHeader = nvt.getValueText();

						UShort min = UShort.ZERO;
						UShort max = UShort.ZERO;
						// Find the max and min parameters and append to the string.
						nvt = gdh.findNameValParam(AffymetrixParameterConsts.MAX_PIXEL_INTENSITY_PARAM_NAME);
						if (nvt != null) {
							if (nvt.getParameterType() == ParameterType.UInt16Type) {
								max = nvt.getValueUInt16();
							}
						}
						nvt = gdh.findNameValParam(AffymetrixParameterConsts.MIN_PIXEL_INTENSITY_PARAM_NAME);
						if (nvt != null) {
							if (nvt.getParameterType() == ParameterType.UInt16Type) {
								min = nvt.getValueUInt16();
							}
						}
						return "[" + Integer.toString(min.toInt()) + ".." + Integer.toString(max.toInt()) + "]" + partialDatHeader;
					}
				}
			}
		}
		return null;
	}

	/** Cell margin */
	public int getCellMargin() throws UnsignedOutOfLimitsException {
		if (gcosFile != null) {
			return gcosFile.getHeader().getMargin();
		}
		else if (calvinFile != null) {
			ParameterNameValue nvt = calvinFile.findAlgorithmParameter(CELAlgorithmParameterNames.CELLMARGIN_PARAM_NAME);
			if (nvt != null) {
				switch (nvt.getParameterType()) {
				case Int32Type:
					return nvt.getValueInt32();
				case Int16Type:
					return nvt.getValueInt16();
				case Int8Type:
					return nvt.getValueInt8();
				case UInt32Type:
					return (int)nvt.getValueUInt32().toLong();
				case UInt16Type:
					return nvt.getValueUInt16().toInt();
				case UInt8Type:
					return nvt.getValueUInt8().toShort();
				case AsciiType:
					return Integer.parseInt(nvt.toString());
				default:
					return 0;
				}
			}
		}
		return 0;
	}

	/** Number of outliers */
	public int getNumOutliers() {
		if (gcosFile != null) {
			return gcosFile.getHeader().getOutliers();
		}
		else if (calvinFile != null) {
			return calvinFile.getNumOutliers();
		}
		return 0;
	}

	/** Number of masked cells */
	public int getNumMasked() {
		if (gcosFile != null) {
			return gcosFile.getHeader().getMasked();
		}
		else if (calvinFile != null) {
			return calvinFile.getNumMasked();
		}
		return 0;
	}

	/**
	 * Get the grid coordinates.
	 * 
	 * @return Returns the grid coordinates.
	 */
	public FGridCoords getGridCorners() {
		if (gcosFile != null) {
			GridCoordinates grid = gcosFile.getHeader().getGrid();
			FGridCoords fgrid = new FGridCoords();
			FPoint fpt = new FPoint();
			fpt.setX(grid.getUpperLeft().getX());
			fpt.setY(grid.getUpperLeft().getY());
			fgrid.setUpperLeft(fpt);
			fpt = new FPoint();
			fpt.setX(grid.getUpperRight().getX());
			fpt.setY(grid.getUpperRight().getY());
			fgrid.setUpperRight(fpt);
			fpt = new FPoint();
			fpt.setX(grid.getLowerRight().getX());
			fpt.setY(grid.getLowerRight().getY());
			fgrid.setLowerRight(fpt);
			fpt = new FPoint();
			fpt.setX(grid.getLowerLeft().getX());
			fpt.setY(grid.getLowerLeft().getY());
			fgrid.setLowerLeft(fpt);
			return fgrid;
		}
		else if (calvinFile != null) {
			FGridCoords grid = new FGridCoords();
			ParameterNameValue nvt;

			nvt = calvinFile.findAlgorithmParameter(CELAlgorithmParameterNames.GRIDULX_PARAM_NAME);
			grid.getUpperLeft().setX(nvt.getValueFloat());

			nvt = calvinFile.findAlgorithmParameter(CELAlgorithmParameterNames.GRIDULY_PARAM_NAME);
			grid.getUpperLeft().setY(nvt.getValueFloat());

			nvt = calvinFile.findAlgorithmParameter(CELAlgorithmParameterNames.GRIDURX_PARAM_NAME);
			grid.getUpperRight().setX(nvt.getValueFloat());

			nvt = calvinFile.findAlgorithmParameter(CELAlgorithmParameterNames.GRIDURY_PARAM_NAME);
			grid.getUpperRight().setY(nvt.getValueFloat());

			nvt = calvinFile.findAlgorithmParameter(CELAlgorithmParameterNames.GRIDLRX_PARAM_NAME);
			grid.getLowerRight().setX(nvt.getValueFloat());

			nvt = calvinFile.findAlgorithmParameter(CELAlgorithmParameterNames.GRIDLRY_PARAM_NAME);
			grid.getLowerRight().setY(nvt.getValueFloat());

			nvt = calvinFile.findAlgorithmParameter(CELAlgorithmParameterNames.GRIDLLX_PARAM_NAME);
			grid.getLowerLeft().setX(nvt.getValueFloat());

			nvt = calvinFile.findAlgorithmParameter(CELAlgorithmParameterNames.GRIDLLY_PARAM_NAME);
			grid.getLowerLeft().setY(nvt.getValueFloat());

			return grid;
		}
		return null;
	}

	/**
	 * Gets the X coordinates from index.
	 * 
	 * @param index
	 *          The 0 based index to the entry array.
	 * @return X coordinate
	 */
	public int indexToX(int index) {
		if (gcosFile != null) {
			return gcosFile.indexToX(index);
		}
		else if (calvinFile != null) {
			return calvinFile.indexToX(index);
		}
		return 0;
	}

	/**
	 * Gets the Y coordinates from index.
	 * 
	 * @param index
	 *          The 0 based index to the entry array.
	 * @return Y coordinate
	 */
	public int indexToY(int index) {
		if (gcosFile != null) {
			return gcosFile.indexToY(index);
		}
		else if (calvinFile != null) {
			return calvinFile.indexToY(index);
		}
		return 0;
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
		if (gcosFile != null) {
			return gcosFile.xyToIndex(x, y);
		}
		else if (calvinFile != null) {
			return calvinFile.xyToIndex(x, y);
		}
		return 0;
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
		return CELFileData.xyToIndex(x, y, r, c);
	}

	/**
	 * Retrieves a CEL file entry.
	 * 
	 * @param index
	 *          The index to the CEL file entries.
	 * @param entry
	 *          The CEL file entry.
	 */
	public void getEntry(int index, FusionCELFileEntryType entry) throws UnsignedOutOfLimitsException {
		entry.setIntensity(getIntensity(index));
		entry.setStdv(getStdv(index));
		entry.setPixels(getPixels(index));
	}

	/**
	 * Retrieves a CEL file entry.
	 * 
	 * @param x
	 *          The X coordinate.
	 * @param y
	 *          The Y coordinate.
	 * @param entry
	 *          The CEL file entry.
	 */
	public void getEntry(int x, int y, FusionCELFileEntryType entry) throws UnsignedOutOfLimitsException {
		getEntry(xyToIndex(x, y), entry);
	}

	/**
	 * Retrieves a CEL file intensity.
	 * 
	 * @param index
	 *          The index to the CEL file entries.
	 * @return The CEL file intensity.
	 */
	public float getIntensity(int index) throws UnsignedOutOfLimitsException {
		if (gcosFile != null) {
			return gcosFile.getIntensity(index);
		}
		else if (calvinFile != null) {
			return calvinFile.getIntensity(index);
		}
		return 0.0f;
	}

	/**
	 * Retrieves a CEL file intensity.
	 * 
	 * @param x
	 *          The X coordinate.
	 * @param y
	 *          The Y coordinate.
	 * @return The CEL file intensity.
	 */
	public float getIntensity(int x, int y) throws UnsignedOutOfLimitsException {
		return getIntensity(xyToIndex(x, y));
	}

	/**
	 * Retrieves a CEL file stdv value.
	 * 
	 * @param index
	 *          The index to the CEL file entries.
	 * @return The CEL file stdv value.
	 */
	public float getStdv(int index) throws UnsignedOutOfLimitsException {
		if (gcosFile != null) {
			return gcosFile.getStdv(index);
		}
		else if (calvinFile != null) {
			return calvinFile.getStdv(index);
		}
		return 0.0f;
	}

	/**
	 * Retrieves a CEL file stdv value.
	 * 
	 * @param x
	 *          The X coordinate.
	 * @param y
	 *          The Y coordinate.
	 * @return The CEL file stdv value.
	 */
	public float getStdv(int x, int y) throws UnsignedOutOfLimitsException {
		return getStdv(xyToIndex(x, y));
	}

	/**
	 * Retrieves a CEL file pixel count.
	 * 
	 * @param index
	 *          The index to the CEL file entries.
	 * @return The CEL file pixel count.
	 */
	public short getPixels(int index) throws UnsignedOutOfLimitsException {
		if (gcosFile != null) {
			return gcosFile.getPixels(index);
		}
		else if (calvinFile != null) {
			return calvinFile.getPixels(index);
		}
		return 0;
	}

	/**
	 * Retrieves a CEL file pixel count.
	 * 
	 * @param x
	 *          The X coordinate.
	 * @param y
	 *          The Y coordinate.
	 * @return The CEL file pixel count.
	 */
	public short getPixels(int x, int y) throws UnsignedOutOfLimitsException {
		return getPixels(xyToIndex(x, y));
	}

	/**
	 * Retrieves a CEL file mask flag.
	 * 
	 * @param x
	 *          The X coordinate.
	 * @param y
	 *          The Y coordinate.
	 * @return True if the feature is masked.
	 */
	public boolean isMasked(int x, int y) {
		return isMasked(xyToIndex(x, y));
	}

	/**
	 * Retrieves a CEL file mask flag.
	 * 
	 * @param index
	 *          The index to the CEL file entries.
	 * @return True if the feature is masked.
	 */
	public boolean isMasked(int index) {
		if (gcosFile != null) {
			return gcosFile.isMasked(index);
		}
		else if (calvinFile != null) {
			return calvinFile.isMasked(index);
		}
		return false;
	}

	/**
	 * Retrieves a CEL file outlier flag.
	 * 
	 * @param x
	 *          The X coordinate.
	 * @param y
	 *          The Y coordinate.
	 * @return True if the feature is an outlier.
	 */
	public boolean isOutlier(int x, int y) {
		return isOutlier(xyToIndex(x, y));
	}

	/**
	 * Retrieves a CEL file outlier flag.
	 * 
	 * @param index
	 *          The index to the CEL file entries.
	 * @return True if the feature is an outlier.
	 */
	public boolean isOutlier(int index) {
		if (gcosFile != null) {
			return gcosFile.isOutlier(index);
		}
		else if (calvinFile != null) {
			return calvinFile.isOutlier(index);
		}
		return false;
	}

	/**
	 * Checks if the file exists.
	 * 
	 * @return True if the file exists.
	 */
	public boolean exists() {
		return new File(fileName).exists();
	}

	/** Creates the GCOS or Calvin parser object as needed. */
	private void createFileParser() {
		clear();
		CELFileData gcosCel = new CELFileData();
		gcosCel.setFileName(fileName);
		if (gcosCel.isVersion3CompatibleFile() || gcosCel.isXDACompatibleFile()) {
			gcosFile = gcosCel;
			return;
		}
		gcosCel = null;

		calvinFile = new CELData();
		CELFileReader reader = new CELFileReader();
		reader.setFilename(fileName);
		try {
			reader.read(calvinFile);
		}
		catch (Throwable t) {
			calvinFile = null;
		}
	}

	/**
	 * Reads the header of the CEL file.
	 * 
	 * @return True if successful.
	 */
	public boolean readHeader() {
		createFileParser();
		if (gcosFile != null) {
			return gcosFile.readHeader();
		}
		else if (calvinFile != null) {
			return true;
		}
		return false;
	}

	/**
	 * Reads the CEL file.
	 * 
	 * @return True if successful.
	 */
	public boolean read() {
		createFileParser();
		if (gcosFile != null) {
			return gcosFile.read();
		}
		else if (calvinFile != null) {
			return true;
		}
		return false;
	}

	/**
	 * Reads the CEL file.
	 * 
	 * @param bIncludeMaskAndOutliers
	 *          Flag to indicate if the mask and outlier sections should also be read.
	 * @return True if successful.
	 */
	public boolean read(boolean bIncludeMaskAndOutliers) {
		createFileParser();
		if (gcosFile != null) {
			return gcosFile.read(bIncludeMaskAndOutliers);
		}
		else if (calvinFile != null) {
			return true;
		}
		return false;
	}

	/** Clears the members. */
	public void clear() {
		gcosFile = null;
		calvinFile = null;
	}

	/** Close the members. */
	public void close() {
		if (gcosFile != null) {
			gcosFile.clear();
		}
		gcosFile = null;
		if (calvinFile != null) {
			calvinFile.clear();
		}
		calvinFile = null;
	}

	/** Creates a new instance of CELFileData */
	public FusionCELData() {
		clear();
	}
}
