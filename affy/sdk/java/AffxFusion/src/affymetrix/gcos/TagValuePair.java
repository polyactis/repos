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

package affymetrix.gcos;

/** Defines name/value parameter. */
public class TagValuePair {

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

	/** Creates a new instance of TagValuePair */
	public TagValuePair() {
		tag = "";
		value = "";
	}

	/**
	 * Creates a new instance of TagValuePair
	 * 
	 * @param p
	 *          The parameter to copy.
	 */
	public TagValuePair(TagValuePair p) {
		tag = p.getTag();
		value = p.getValue();
	}

	/**
	 * Makes a copy of the object.
	 * 
	 * @return A copy of the object.
	 */
	public TagValuePair copy() {
		TagValuePair p = new TagValuePair();
		p.setTag(tag);
		p.setValue(value);
		return p;
	}
}
