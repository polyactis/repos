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

package affymetrix.fusion.chp;

import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.data.CHPBackgroundZone;
import affymetrix.calvin.data.CHPData;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.fusion.FusionTagValuePair;
import affymetrix.gcos.TagValuePair;
import affymetrix.gcos.chp.BackgroundZoneInfo;
import affymetrix.gcos.chp.BackgroundZoneType;
import affymetrix.gcos.chp.CHPFileHeader;

/** Stores the contents of the header of a CHP file. */
public class FusionCHPHeader {

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

	/** The GCOS header object. */
	private CHPFileHeader gcosHeader;

	/**
	 * Set the gcos header object.
	 * 
	 * @param header
	 *          The GCOS header object.
	 */
	public void setGCOSHeader(CHPFileHeader header) {
		gcosHeader = header;
	}

	/** The calvin object. */
	private CHPData calvinData;

	/**
	 * Set the calvin object.
	 * 
	 * @param data
	 *          The calvin object.
	 */
	public void setCalvinObject(CHPData data) {
		calvinData = data;
	}

	/**
	 * Gets the number of columns of features.
	 * 
	 * @return The number of columns of features.
	 */
	public int getCols() {
		if (gcosHeader != null) {
			return gcosHeader.getCols();
		}
		else if (calvinData != null) {
			return calvinData.getCols();
		}
		return 0;
	}

	/**
	 * Gets the number of rows of features.
	 * 
	 * @return The number of rows of features.
	 */
	public int getRows() {
		if (gcosHeader != null) {
			return gcosHeader.getRows();
		}
		else if (calvinData != null) {
			return calvinData.getRows();
		}
		return 0;
	}

	/**
	 * Gets the number of probe sets.
	 * 
	 * @return The number of probe sets.
	 */
	public int getNumProbeSets() {
		if (gcosHeader != null) {
			return gcosHeader.getNumProbeSets();
		}
		else if (calvinData != null) {
			return calvinData.getEntryCount();
		}
		return 0;
	}

	/**
	 * Gets the assay type.
	 * 
	 * @return The assay type.
	 */
	public int getAssayType() {
		if (gcosHeader != null) {
			return gcosHeader.getAssayType();
		}
		else if (calvinData != null) {
			String str = calvinData.getAssayType();
			if (str.equals(CHPData.CHP_EXPRESSION_ASSAY_TYPE)) {
				return EXPRESSION_ASSAY;
			}
			if (str.equals(CHPData.CHP_GENOTYPING_ASSAY_TYPE)) {
				return GENOTYPING_ASSAY;
			}
			if (str.equals(CHPData.CHP_UNIVERSAL_ASSAY_TYPE)) {
				return UNIVERSAL_ASSAY;
			}
			if (str.equals(CHPData.CHP_RESEQUENCING_ASSAY_TYPE)) {
				return RESEQUENCING_ASSAY;
			}
		}
		return UNKNOWN_ASSAY;
	}

	/**
	 * Gets the chip type.
	 * 
	 * @return The chip type.
	 */
	public String getChipType() {
		if (gcosHeader != null) {
			return gcosHeader.getChipType();
		}
		else if (calvinData != null) {
			return calvinData.getArrayType();
		}
		return null;
	}

	/**
	 * Gets the algorithm name.
	 * 
	 * @return The algorithm name.
	 */
	public String getAlgName() {
		if (gcosHeader != null) {
			return gcosHeader.getAlgName();
		}
		else if (calvinData != null) {
			return calvinData.getAlgName();
		}
		return null;
	}

	/**
	 * Gets the algorithm version.
	 * 
	 * @return The algorithm version.
	 */
	public String getAlgVersion() {
		if (gcosHeader != null) {
			return gcosHeader.getAlgVersion();
		}
		else if (calvinData != null) {
			return calvinData.getAlgVersion();
		}
		return null;
	}

	/**
	 * Gets the name of the CEL file used to create the CHP file.
	 * 
	 * @return The CEL file name.
	 */
	public String getParentCellFile() {
		if (gcosHeader != null) {
			return gcosHeader.getParentCellFile();
		}
		else if (calvinData != null) {
			return calvinData.getParentCell();
		}
		return null;
	}

	/**
	 * Gets the algorithm prog ID (COM components only).
	 * 
	 * @return The id.
	 */
	public String getProgID() {
		if (gcosHeader != null) {
			return gcosHeader.getProgID();
		}
		else if (calvinData != null) {
			return calvinData.getProgId();
		}
		return null;
	}

	/**
	 * Gets the algorithm parameters.
	 * 
	 * @return The parameters.
	 */
	public List<FusionTagValuePair> getAlgorithmParameters() {
		if (gcosHeader != null) {
			List<TagValuePair> gParams = gcosHeader.getAlgorithmParameters();
			if (gParams == null) {
				return null;
			}
			int n = gParams.size();
			List<FusionTagValuePair> params = new ArrayList<FusionTagValuePair>();
			for (int i = 0; i < n; i++) {
				TagValuePair gparam = gParams.get(i);
				FusionTagValuePair param = new FusionTagValuePair();
				param.setTag(gparam.getTag());
				param.setValue(gparam.getValue());
				params.add(param);
			}
			return params;
		}
		else if (calvinData != null) {
			List<ParameterNameValue> cParams = calvinData.getAlgParams();
			if (cParams == null) {
				return null;
			}
			int n = cParams.size();
			List<FusionTagValuePair> params = new ArrayList<FusionTagValuePair>();
			for (int i = 0; i < n; i++) {
				ParameterNameValue cparam = cParams.get(i);
				FusionTagValuePair param = new FusionTagValuePair();
				param.setTag(cparam.getName());
				param.setValue(cparam.toString());
				param.setDetailed(cparam);
				params.add(param);
			}
			return params;
		}
		return null;
	}

	/**
	 * Gets the summary parameters.
	 * 
	 * @return The parameters.
	 */
	public List<FusionTagValuePair> getSummaryParameters() {
		if (gcosHeader != null) {
			List<TagValuePair> gParams = gcosHeader.getSummaryParameters();
			if (gParams == null) {
				return null;
			}
			int n = gParams.size();
			List<FusionTagValuePair> params = new ArrayList<FusionTagValuePair>();
			for (int i = 0; i < n; i++) {
				TagValuePair gparam = gParams.get(i);
				FusionTagValuePair param = new FusionTagValuePair();
				param.setTag(gparam.getTag());
				param.setValue(gparam.getValue());
				params.add(param);
			}
			return params;
		}
		else if (calvinData != null) {
			List<ParameterNameValue> cParams = calvinData.getChipSums();
			if (cParams == null) {
				return null;
			}
			int n = cParams.size();
			List<FusionTagValuePair> params = new ArrayList<FusionTagValuePair>();
			for (int i = 0; i < n; i++) {
				ParameterNameValue cparam = cParams.get(i);
				FusionTagValuePair param = new FusionTagValuePair();
				param.setTag(cparam.getName());
				param.setValue(cparam.toString());
				param.setDetailed(cparam);
				params.add(param);
			}
			return params;
		}
		return null;
	}

	/**
	 * Gets the background zone information.
	 * 
	 * @return The background information.
	 */
	public BackgroundZoneInfo getBackgroundZoneInfo() throws UnsignedOutOfLimitsException {
		if (gcosHeader != null) {
			return gcosHeader.getBackgroundZoneInfo();
		}
		else if (calvinData != null) {
			CHPBackgroundZone calvinZone = new CHPBackgroundZone();
			BackgroundZoneInfo info = new BackgroundZoneInfo();
			int n = calvinData.getBackgroundZoneCnt();
			for (int i = 0; i < n; i++) {
				calvinData.getBackgroundZone(i, calvinZone);
				info.setSmoothFactor(calvinZone.getSmoothFactor());
				BackgroundZoneType zone = new BackgroundZoneType();
				zone.setBackground(calvinZone.getBackground());
				zone.setCenterX(calvinZone.getCenterX());
				zone.setCenterY(calvinZone.getCenterY());
				info.addZone(zone);
			}
			return info;
		}
		return null;
	}

	/**
	 * Gets a specific algorithm parameter given a name/tag
	 * 
	 * @return The specific algorithm parameter given a name/tag
	 */
	public String getAlgorithmParameter(String tag) {
		if (gcosHeader != null) {
			return gcosHeader.getAlgorithmParameter(tag);
		}
		else if (calvinData != null) {
			ParameterNameValue nvt = calvinData.getAlgParam(tag);
			if (nvt != null) {
				return nvt.toString();
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
		if (gcosHeader != null) {
			return gcosHeader.getSummaryParameter(tag);
		}
		else if (calvinData != null) {
			ParameterNameValue nvt = calvinData.getChipSum(tag);
			if (nvt != null) {
				return nvt.toString();
			}
		}
		return null;
	}

	/**
	 * Gets the background value for a given center coordinate
	 * 
	 * @return The background value for a given center coordinate
	 */
	public BackgroundZoneType getBackgroundZone(int x, int y) throws UnsignedOutOfLimitsException {
		if (gcosHeader != null) {
			return gcosHeader.getBackgroundZone(x, y);
		}
		else if (calvinData != null) {
			CHPBackgroundZone zn = new CHPBackgroundZone();
			int n = calvinData.getBackgroundZoneCnt();
			for (int i = 0; i < n; i++) {
				calvinData.getBackgroundZone(i, zn);
				if ((Float.compare(zn.getCenterX(), x) == 0) && (Float.compare(zn.getCenterY(), y) == 0)) {
					BackgroundZoneType bg = new BackgroundZoneType();
					bg.setBackground(zn.getBackground());
					bg.setCenterX((int)zn.getCenterX());
					bg.setCenterY((int)zn.getCenterY());
					return bg;
				}
			}
		}
		return null;
	}

	/** Creates a new instance of FusionCHPHeader */
	public FusionCHPHeader() {
		gcosHeader = null;
		calvinData = null;
	}
}
