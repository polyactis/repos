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

package affymetrix.fusion.cdf;

import affymetrix.gcos.cdf.CDFProbeSetNames;

/** This class provides storage for the list of probe set names. */
public class FusionCDFProbeSetNames {

	/** The GCOS probe set names object. */
	private CDFProbeSetNames gcosNames;

	/**
	 * The GCOS probe set names object.
	 * 
	 * @param n
	 *          The probe set names object.
	 */
	public void setGCOSObject(final CDFProbeSetNames n) {
		gcosNames = n;
	}

	/**
	 * Gets the probe set name.
	 * 
	 * @param index
	 *          The index to the probe set name of interest.
	 * @return The name of the probe set.
	 */
	public String getName(final int index) {
		if (gcosNames != null) {
			return gcosNames.getName(index);
		}
		return null;
	}

	/** Creates a new instance of FusionCDFProbeSetNames */
	public FusionCDFProbeSetNames() {
		clear();
	}

	/** Clears the members. */
	public void clear() {
		gcosNames = null;
	}
}
