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

import affymetrix.gcos.cdf.CDFProbeInformation;

/** This class provides storage for an individual probe in a CDF file. */
public class FusionCDFProbeInformation {

	/** The GCOS probe object. */
	private CDFProbeInformation gcosProbe;

	/**
	 * Sets the GCOS probe object.
	 * 
	 * @param p
	 *          The probe object.
	 */
	public void setGCOSObject(final CDFProbeInformation p) {
		gcosProbe = p;
	}

	/** Gets the index of the probes probe pair or quartet (probe list). Also known as the atom position. */
	public int getListIndex() {
		if (gcosProbe != null) {
			return gcosProbe.getListIndex();
		}
		return 0;
	}

	/**
	 * Gets the expos value in the CDF file, this can either be a zero based index equal to the list index or the exon
	 * position.
	 */
	public int getExpos() {
		if (gcosProbe != null) {
			return gcosProbe.getExpos();
		}
		return 0;
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

	/** Gets the probes base at the interrogation position */
	public char getPBase() {
		if (gcosProbe != null) {
			return gcosProbe.getPBase();
		}
		return ' ';
	}

	/** Gets the targets base at the interrogation position */
	public char getTBase() {
		if (gcosProbe != null) {
			return gcosProbe.getTBase();
		}
		return ' ';
	}

	/** Creates a new instance of FusionCDFProbeInformation */
	public FusionCDFProbeInformation() {
		clear();
	}

	/** Clears the members. */
	public void clear() {
		gcosProbe = null;
	}

}
