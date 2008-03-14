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

import affymetrix.calvin.parameter.AffymetrixParameterConsts;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.utils.AffymetrixGuidType;
import affymetrix.calvin.utils.IOUtils;

/**
 * 
 * @author ljevon
 */
public class GenericDataHeader {

	/** Creates a new instance of GenericDataHeader */
	public GenericDataHeader() {
		clear();
		locale = AffymetrixParameterConsts.US_ENGLISH_LOCALE;
	}

	/**  */
	private String fileTypeId;

	/**  */
	private AffymetrixGuidType fileId = new AffymetrixGuidType();

	/**  */
	private String fileCreationTime;

	/**  */
	private String locale;

	/**  */
	private List<ParameterNameValue> nameValParams = new ArrayList<ParameterNameValue>();

	/**  */
	private List<GenericDataHeader> genericDataHdrs = new ArrayList<GenericDataHeader>();

	/**  */
	public void clear() {
		fileTypeId = IOUtils.EMPTY;
		fileId.clear();
		fileCreationTime = IOUtils.EMPTY;
		nameValParams.clear();
		genericDataHdrs.clear();
	}

	/**  */
	public void setFileTypeId(String p) {
		fileTypeId = p;
	}

	/**  */
	public String getFileTypeId() {
		return fileTypeId;
	}

	/**  */
	public void setFileId(AffymetrixGuidType p) {
		fileId = p;
	}

	public void setFileId(byte[] affyGuid) {
		fileId.setGuid(affyGuid);
	}

	/**  */
	public AffymetrixGuidType getFileId() {
		return fileId;
	}

	/**  */
	public void setFileCreationTime(String f) {
		fileCreationTime = f;
	}

	/**  */
	public String getFileCreationTime() {
		return fileCreationTime;
	}

	/**  */
	public void setLocale(String p) {
		locale = p;
	}

	/**  */
	public String getLocale() {
		return locale;
	}

	/**  */
	public void addNameValParam(ParameterNameValue p) {
		nameValParams.add(p);
	}

	/**  */
	public ParameterNameValue getNameValParam(int index) {
		if (nameValParams.size() > 0) {
			return nameValParams.get(index);
		}
		return null;
	}

	/**  */
	public int getNameValParamCnt() {
		return nameValParams.size();
	}

	/**  */
	public List<ParameterNameValue> getNameValParams() {
		return nameValParams;
	}

	/**  */
	public int getParentCnt() {
		return genericDataHdrs.size();
	}

	/**  */
	public void addParent(GenericDataHeader hdr) {
		genericDataHdrs.add(hdr);
	}

	/**  */
	public GenericDataHeader getParent(int index) {
		if (genericDataHdrs.size() > 0) {
			return genericDataHdrs.get(index);
		}
		return null;
	}

	/**  */
	public List<GenericDataHeader> getParents() {
		return genericDataHdrs;
	}

	/**
	 * Finds a ParameterNameValue by name in the nameValPairs collection
	 * 
	 * @param name
	 *          The name of the NameValPair to find
	 * @return Reference to a ParameterNameValue with the found ParameterNameValue, null if not found.
	 */
	public ParameterNameValue findNameValParam(String name) {
		int n = getNameValParamCnt();
		for (int i = 0; i < n; i++) {
			ParameterNameValue param = getNameValParam(i);
			if (name.equals(param.getName())) {
				return param;
			}
		}
		return null;
	}

	/**
	 * Find an immediate parent GenericDataHeader based on file type id. Does not search grand-parents or above.
	 * 
	 * @param fileTypeId
	 *          The fileTypeId of the parent header to find.
	 * @return Returns a pointer to the parent GenericDataHeader if found, otherwise returns 0.
	 */
	public GenericDataHeader findParent(String fileTypeId) {
		List<GenericDataHeader> parents = getParents();
		for (int i = 0; i < parents.size(); i++) {
			GenericDataHeader p = parents.get(i);
			if (fileTypeId.equals(p.getFileTypeId())) {
				return p;
			}
		}
		return null;
	}

	/**
	 * Gets all ParameterNameValue where the name starts with a given string.
	 * 
	 * @param beginsWith
	 *          The string that the beginning of the ParameterNameValue name needs to match.
	 * @return A result vector of ParameterNameValues where the name begins with the beginsWith argument.
	 */
	public List<ParameterNameValue> getNameValParamsBeginsWith(String beginsWith) {

		List<ParameterNameValue> p = null;
		int n = getNameValParamCnt();
		for (int i = 0; i < n; i++) {
			ParameterNameValue param = getNameValParam(i);
			if (param.getName().startsWith(beginsWith)) {
				if (p == null) {
					p = new ArrayList<ParameterNameValue>();
				}
				p.add(param);
			}
		}
		return p;
	}

}
