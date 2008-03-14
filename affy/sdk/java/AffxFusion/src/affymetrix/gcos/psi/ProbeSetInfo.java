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

package affymetrix.gcos.psi;

/**
 * This class stores information from a PSI file.
 */
public class ProbeSetInfo {

	/** The name of the probe set. */
	private String probeSetName;

	/**
	 * Sets the name of the probe set.
	 * 
	 * @param name
	 *          The probe set name.
	 */
	public void setProbeSetName(String name) {
		probeSetName = name;
	}

	/**
	 * Gets the name of the probe set.
	 * 
	 * @return The probe set name.
	 */
	public String getProbeSetName() {
		return probeSetName;
	}

	/** The number of probe pairs in the set. */
	private int numberPairs;

	/**
	 * Sets the number of probe pairs in the set.
	 * 
	 * @param n
	 *          The number of pairs.
	 */
	public void setNumberPairs(int n) {
		numberPairs = n;
	}

	/**
	 * Gets the number of probe pairs in the set.
	 * 
	 * @return The number of pairs.
	 */
	public int getNumberPairs() {
		return numberPairs;
	}

	/** Creates a new instance of ProbeSetInfo */
	public ProbeSetInfo() {
		probeSetName = "";
		numberPairs = 0;
	}
}
