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

import affymetrix.gcos.cdf.CDFQCProbeInformation;

/** This class provides storage for QC probes from a CDF file. */
public class FusionCDFQCProbeInformation {

	/** The gcos object. */
	private CDFQCProbeInformation gcosProbe;

	/**
	 * Sets the GCOS QC probe object.
	 * 
	 * @param p
	 *          The probe object.
	 */
	public void setGCOSObject(final CDFQCProbeInformation p) {
		gcosProbe = p;
	}

	/** Gets the X coordinate */
	public int getX() {
		if (gcosProbe != null) {
			return gcosProbe.getX();
		}
		return 0;
	}

	/** Gets the Y coordinate */
	public int getY() {
		if (gcosProbe != null) {
			return gcosProbe.getY();
		}
		return 0;
	}

	/** Gets the probe length. This value may be 1 for non-synthesized features. */
	public byte getProbeLength() {
		if (gcosProbe != null) {
			return gcosProbe.getProbeLength();
		}
		return 0;
	}

	/** Gets a flag indicating if the probe is a perfect match probe. */
	public boolean isPMProbe() {
		if (gcosProbe != null) {
			return gcosProbe.isPMProbe();
		}
		return false;
	}

	/** Gets a flag indicating if the probe is used for background calculations (blank feature). */
	public boolean isBackground() {
		if (gcosProbe != null) {
			return gcosProbe.isBackground();
		}
		return false;
	}

	/** Creates a new instance of FusionCDFQCProbeInformation */
	public FusionCDFQCProbeInformation() {
		clear();
	}

	/** Clears the members. */
	public void clear() {
		gcosProbe = null;
	}
}
