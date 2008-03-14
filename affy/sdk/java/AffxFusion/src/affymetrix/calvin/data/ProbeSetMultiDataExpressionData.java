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

/**
 * Stores data for a expression result of a probe set.
 */
public class ProbeSetMultiDataExpressionData extends ProbeSetMultiDataBase {

	/** The name of the probe set. */
	private String name;

	/** The quantification of the call. */
	private float quantification;

	/**
	 * Creates a new instance of ProbeSetMultiDataExpressionData.
	 */
	public ProbeSetMultiDataExpressionData() {
		name = null;
		quantification = 0.0f;
	}

	/**
	 * The name of the probe set.
	 * 
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * Set the name.
	 * 
	 * @param n
	 *          the name
	 */
	public void setName(String n) {
		name = n;
	}

	/**
	 * Get the quantification.
	 * 
	 * @return the quantification
	 */
	public float getQuantification() {
		return quantification;
	}

	/**
	 * Set the quantification.
	 * 
	 * @param c
	 *          the quantification
	 */
	public void setQuantification(float c) {
		quantification = c;
	}
}
