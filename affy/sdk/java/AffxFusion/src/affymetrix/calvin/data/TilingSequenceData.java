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

import affymetrix.calvin.parameter.ParameterNameValue;

/** Stores data for a tiling array sequence. */
public class TilingSequenceData {

	/** The name of the sequence. */
	private String name;

	/** The version associated to the sequence. */
	private String version;

	/** The group name for the sequence. */
	private String groupName;

	/** The parameter name/value array. */
	private List<ParameterNameValue> parameters = new ArrayList<ParameterNameValue>();

	/** Creates a new instance of TilingSequenceData */
	public TilingSequenceData() {
		name = null;
		version = null;
		groupName = null;
	}

	/** The name of the sequence. */
	public String getName() {
		return name;
	}

	/** The name of the sequence. */
	public void setName(String n) {
		name = n;
	}

	/** The version associated to the sequence. */
	public String getVersion() {
		return version;
	}

	/** The version associated to the sequence. */
	public void setVersion(String v) {
		version = v;
	}

	/** The group name for the sequence. */
	public String getGroupName() {
		return groupName;
	}

	/** The group name for the sequence. */
	public void setGroupName(String g) {
		groupName = g;
	}

	/** The parameter name/value array. */
	public List<ParameterNameValue> getParameters() {
		return parameters;
	}

	public void clearParameters() {
		parameters.clear();
	}

	public void addParameter(ParameterNameValue v) {
		parameters.add(v);
	}

	/** The parameter name/value array. */
	public void addParameters(List<ParameterNameValue> v) {
		parameters.addAll(v);
	}
}
