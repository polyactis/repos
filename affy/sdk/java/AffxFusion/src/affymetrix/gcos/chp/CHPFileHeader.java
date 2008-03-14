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

package affymetrix.gcos.chp;

import java.util.List;

import affymetrix.calvin.utils.IOUtils;
import affymetrix.gcos.TagValuePair;

/** Stores the contents of the header of a CHP file. */
public class CHPFileHeader {

	/** Expression assay type. */
	public static final int EXPRESSION_ASSAY = 0;

	/** Genotyping assay type. */
	public static final int GENOTYPING_ASSAY = 1;

	/** Resequencing assay type. */
	public static final int RESEQUENCING_ASSAY = 2;

	/** Universal assay type. */
	public static final int UNIVERSAL_ASSAY = 3;

	/** Unknown assay type. */
	public static final int UNKNOWN_ASSAY = 4;

	/** The magic number in the file */
	private int magic;

	/** The version number in the file */
	private int version;

	/** The number of feature rows in the array */
	private int rows;

	/** The number of feature columns in the array */
	private int cols;

	/** The number of probe set results */
	private int numProbeSets;

	/** The type of results stored in the CHP file */
	private int assayType;

	/** The chip type or probe array type of the CHP file */
	private String chipType;

	/** The version number of the algorithm used to create the CHP file */
	private String algorithmVersion;

	/** The name of the algorithm used to create the CHP file */
	private String algorithmName;

	/** The programmatic identifier of the algorithm used to create the CHP file */
	private String progID;

	/** The vector of algorithm parameters */
	private List<TagValuePair> algorithmParameters;

	/** The background's for each of the zones (calculated by the expression algorithm) */
	private BackgroundZoneInfo backgroundZoneInfo;

	/**
	 * Gets the magic number in the file.
	 * 
	 * @return The magic number of the file.
	 */
	public int getMagic() {
		return magic;
	}

	/**
	 * Sets the magic number in the file.
	 * 
	 * @param m
	 *          The magic number of the file.
	 */
	public void setMagic(int m) {
		magic = m;
	}

	/**
	 * Gets the version number in the file.
	 * 
	 * @return The version number of the file.
	 */
	public int getVersion() {
		return version;
	}

	/**
	 * Sets the version number in the file.
	 * 
	 * @param v
	 *          The version number of the file.
	 */
	public void setVersion(int v) {
		version = v;
	}

	/**
	 * Gets the number of columns of features.
	 * 
	 * @return The number of columns of features.
	 */
	public int getCols() {
		return cols;
	}

	/**
	 * Sets the number of columns of features.
	 * 
	 * @param c
	 *          The number of columns of features.
	 */
	public void setCols(int c) {
		cols = c;
	}

	/**
	 * Gets the number of rows of features.
	 * 
	 * @return The number of rows of features.
	 */
	public int getRows() {
		return rows;
	}

	/**
	 * Sets the number of rows of features.
	 * 
	 * @param r
	 *          The number of rows of features.
	 */
	public void setRows(int r) {
		rows = r;
	}

	/**
	 * Gets the number of probe sets.
	 * 
	 * @return The number of probe sets.
	 */
	public int getNumProbeSets() {
		return numProbeSets;
	}

	/**
	 * Sets the number of probe sets.
	 * 
	 * @param n
	 *          The number of probe sets.
	 */
	public void setNumProbeSets(int n) {
		numProbeSets = n;
	}

	/**
	 * Gets the assay type.
	 * 
	 * @return The assay type.
	 */
	public int getAssayType() {
		return assayType;
	}

	/**
	 * Sets the assay type.
	 * 
	 * @param t
	 *          The assay type.
	 */
	public void setAssayType(int t) {
		assayType = t;
	}

	/**
	 * Gets the chip type.
	 * 
	 * @return The chip type.
	 */
	public String getChipType() {
		return chipType;
	}

	/**
	 * Sets the chip type.
	 * 
	 * @param str
	 *          The chip type.
	 */
	public void setChipType(String str) {
		chipType = str;
	}

	/**
	 * Gets the algorithm name.
	 * 
	 * @return The algorithm name.
	 */
	public String getAlgName() {
		return algorithmName;
	}

	/**
	 * Sets the algorithm name.
	 * 
	 * @param str
	 *          The algorithm name.
	 */
	public void setAlgName(String str) {
		algorithmName = str;
	}

	/**
	 * Gets the algorithm version.
	 * 
	 * @return The algorithm version.
	 */
	public String getAlgVersion() {
		return algorithmVersion;
	}

	/**
	 * Sets the algorithm version.
	 * 
	 * @param str
	 *          The algorithm version.
	 */
	public void setAlgVersion(String str) {
		algorithmVersion = str;
	}

	/** The name of the CEL file used in the creation of the CHP file */
	private String parentCellFile;

	/**
	 * Gets the name of the CEL file used to create the CHP file.
	 * 
	 * @return The CEL file name.
	 */
	public String getParentCellFile() {
		return parentCellFile;
	}

	/**
	 * Sets the name of the CEL file used to create the CHP file.
	 * 
	 * @param str
	 *          The CEL file name.
	 */
	public void setParentCellFile(String str) {
		parentCellFile = str;
	}

	/**
	 * Gets the algorithm prog ID (COM components only).
	 * 
	 * @return The id.
	 */
	public String getProgID() {
		return progID;
	}

	/**
	 * Sets the prog ID.
	 * 
	 * @param str
	 *          The prog ID.
	 */
	public void setProgID(String str) {
		progID = str;
	}

	/**
	 * Gets the algorithm parameters.
	 * 
	 * @return The parameters.
	 */
	public List<TagValuePair> getAlgorithmParameters() {
		return algorithmParameters;
	}

	/**
	 * Sets the algorithm parameters.
	 * 
	 * @param params
	 *          The parameters.
	 */
	public void setAlgorithmParameters(List<TagValuePair> params) {
		algorithmParameters = params;
	}

	/** A vector of summary parameters generated by the CHP file generating algorithm */
	private List<TagValuePair> summaryParameters;

	/**
	 * Gets the summary parameters.
	 * 
	 * @return The parameters.
	 */
	public List<TagValuePair> getSummaryParameters() {
		return summaryParameters;
	}

	/**
	 * Sets the summary parameters.
	 * 
	 * @param params
	 *          The parameters.
	 */
	public void setSummaryParameters(List<TagValuePair> params) {
		summaryParameters = params;
	}

	/**
	 * Gets the background zone information.
	 * 
	 * @return The background information.
	 */
	public BackgroundZoneInfo getBackgroundZoneInfo() {
		return backgroundZoneInfo;
	}

	/**
	 * Sets the background zone information.
	 * 
	 * @param zones
	 *          The background information.
	 */
	public void setBackgroundZoneInfo(BackgroundZoneInfo zones) {
		backgroundZoneInfo = zones;
	}

	/**
	 * Gets a specific algorithm parameter given a name/tag
	 * 
	 * @return The specific algorithm parameter given a name/tag
	 */
	public String getAlgorithmParameter(String tag) {
		if (algorithmParameters != null) {
			int n = algorithmParameters.size();
			for (int i = 0; i < n; i++) {
				TagValuePair param = algorithmParameters.get(i);
				if (param.getTag().compareTo(tag) == 0) {
					return param.getValue();
				}
			}
		}
		return null;
	}

	/**
	 * Gets a specific summary parameter given a name/tag
	 * 
	 * @return The specific summary parameter given a name/tag
	 */
	public String getSummaryParameter(String tag) {
		if (summaryParameters != null) {
			int n = summaryParameters.size();
			for (int i = 0; i < n; i++) {
				TagValuePair param = summaryParameters.get(i);
				if (param.getTag().compareTo(tag) == 0) {
					return param.getValue();
				}
			}
		}
		return null;
	}

	/**
	 * Gets the background value for a given center coordinate
	 * 
	 * @return The background value for a given center coordinate
	 */
	public BackgroundZoneType getBackgroundZone(int x, int y) {
		int n = backgroundZoneInfo.getNumberZones();
		for (int i = 0; i < n; i++) {
			BackgroundZoneType zone = backgroundZoneInfo.getZone(i);
			if (((int)zone.getCenterX() == x) && ((int)zone.getCenterY() == y)) {
				return zone;
			}
		}
		return null;
	}

	/** Creates a new instance of CHPFileHeader */
	public CHPFileHeader() {
		magic = 0;
		version = 0;
		cols = 0;
		rows = 0;
		numProbeSets = 0;
		assayType = UNKNOWN_ASSAY;
		chipType = IOUtils.EMPTY;
		algorithmName = IOUtils.EMPTY;
		algorithmVersion = IOUtils.EMPTY;
		parentCellFile = IOUtils.EMPTY;
		progID = IOUtils.EMPTY;
		algorithmParameters = null;
		summaryParameters = null;
		backgroundZoneInfo = null;
	}
}
