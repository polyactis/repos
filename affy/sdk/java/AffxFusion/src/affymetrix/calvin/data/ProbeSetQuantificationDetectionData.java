/////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

/** Stores data for a quantification/detection of a probe set. */
public class ProbeSetQuantificationDetectionData {

	/** The name of the probe set. */
	private String name;

	/** The name of the probe set. */
	public String getName() {
		return name;
	}

	/** The name of the probe set. */
	public void setName(String n) {
		name = n;
	}

	/** The probe set id. */
	private int id;

	/** The probe set id. */
	public int getId() {
		return id;
	}

	/** The probe set id. */
	public void setId(int n) {
		id = n;
	}

	/** The quantification associated to the name. */
	private float quantification;

	/** The quantification associated to the name. */
	public float getQuantification() {
		return quantification;
	}

	/** The quantification associated to the name. */
	public void setQuantification(float q) {
		quantification = q;
	}

	/** The detection p-value associated to the name. */
	private float pvalue;

	/** The detection p-value associated to the name. */
	public float getPValue() {
		return pvalue;
	}

	/** The detection p-value associated to the name. */
	public void setPValue(float p) {
		pvalue = p;
	}

	/** Creates a new instance of ProbeSetQuantificationDetectionData */
	public ProbeSetQuantificationDetectionData() {
		name = "";
		id = -1;
		quantification = 0.0f;
		pvalue = 0.0f;
	}
}
