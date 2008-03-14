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

package affymetrix.fusion;

import affymetrix.calvin.parameter.ParameterNameValue;

/** Defines name/value parameter. */
public class FusionTagValuePair {

	/** The name of the parameter. */
	private String tag;

	/**
	 * Gets the parameter name (tag).
	 * 
	 * @return The name of the parameter.
	 */
	public String getTag() {
		return tag;
	}

	/**
	 * Sets the parameter name (tag).
	 * 
	 * @param t
	 *          The name of the parameter.
	 */
	public void setTag(String t) {
		tag = t;
	}

	/** The value of the parameter. */
	private String value;

	/**
	 * Gets the parameter value.
	 * 
	 * @return The value.
	 */
	public String getValue() {
		return value;
	}

	/**
	 * Sets the parameter value.
	 * 
	 * @param v
	 *          The value.
	 */
	public void setValue(String v) {
		value = v;
	}

	/** Embbedded Calvin parameter type object */
	private ParameterNameValue nvt;

	/** Embbedded Calvin parameter type object */
	public ParameterNameValue getDetailed() {
		return nvt;
	}

	/** Embbedded Calvin parameter type object */
	public void setDetailed(ParameterNameValue detailed) {
		nvt = detailed;
	}

	/** Creates a new instance of FusionTagValuePair */
	public FusionTagValuePair() {
		tag = "";
		value = "";
		nvt = null;
	}

	/**
	 * Creates a new instance of FusionTagValuePair
	 * 
	 * @param p
	 *          The parameter to copy.
	 */
	public FusionTagValuePair(FusionTagValuePair p) {
		tag = p.getTag();
		value = p.getValue();
		nvt = p.getDetailed();
	}

	/**
	 * Makes a copy of the object.
	 * 
	 * @return A copy of the object.
	 */
	public FusionTagValuePair copy() {
		FusionTagValuePair p = new FusionTagValuePair();
		p.setTag(tag);
		p.setValue(value);
		if (nvt != null) {
			p.setDetailed(new ParameterNameValue(nvt));
		}
		return p;
	}
}
