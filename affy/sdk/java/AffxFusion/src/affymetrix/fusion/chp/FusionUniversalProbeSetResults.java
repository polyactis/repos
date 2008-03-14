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

import affymetrix.calvin.data.CHPUniversalEntry;
import affymetrix.gcos.chp.UniversalProbeSetResults;

/** Stores results for a universal probe set. */
public class FusionUniversalProbeSetResults {

	/** The GCOS probe set results object. */
	private UniversalProbeSetResults gcosResult;

	/**
	 * Sets the GCOS probe set results object.
	 * 
	 * @param r
	 *          The GCOS probe set results object.
	 */
	public void setGCOSObject(UniversalProbeSetResults r) {
		gcosResult = r;
	}

	/** The calvin probe set results object. */
	private CHPUniversalEntry calvinResult;

	/**
	 * Sets the calvin object.
	 * 
	 * @param r
	 *          The calvin probe set object.
	 */
	public void setCalvinObject(CHPUniversalEntry r) {
		calvinResult = r;
	}

	/** Clears the members. */
	public void clear() {
		gcosResult = null;
		calvinResult = null;
	}

	/**
	 * Gets the background value.
	 * 
	 * @return The background value.
	 */
	public float getBackground() {
		if (gcosResult != null) {
			return gcosResult.getBackground();
		}
		else if (calvinResult != null) {
			return calvinResult.getBackground();
		}
		return 0.0f;
	}

	/** Creates a new instance of UniversalProbeSetResults */
	public FusionUniversalProbeSetResults() {
		clear();
	}

}
