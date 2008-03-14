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

/** Stores results for a universal probe set. */
public class UniversalProbeSetResults {

	/** The background value. */
	private float background;

	/**
	 * Gets the background value.
	 * 
	 * @return The background value.
	 */
	public float getBackground() {
		return background;
	}

	/**
	 * Sets the background value.
	 * 
	 * @param b
	 *          The background value.
	 */
	public void setBackground(float b) {
		background = b;
	}

	/** Creates a new instance of UniversalProbeSetResults */
	public UniversalProbeSetResults() {
		background = 0.0f;
	}

}
