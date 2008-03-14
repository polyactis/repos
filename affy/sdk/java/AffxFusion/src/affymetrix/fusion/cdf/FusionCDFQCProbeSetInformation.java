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

import affymetrix.gcos.cdf.CDFQCProbeSetInformation;

/** This class provides storage for the probes in a QC probe set. */
public class FusionCDFQCProbeSetInformation {

	/** The gcos object. */
	private CDFQCProbeSetInformation gcosProbeSet;

	/**
	 * Sets the GCOS QC probe set object.
	 * 
	 * @param set
	 *          The probe set object.
	 */
	public void setGCOSObject(final CDFQCProbeSetInformation set) {
		gcosProbeSet = set;
	}

	/** Gets the number of cells in the set. */
	public int getNumCells() {
		if (gcosProbeSet != null) {
			return gcosProbeSet.getNumCells();
		}
		return 0;
	}

	/** The type of QC probes. */
	public short getQCProbeSetType() {
		if (gcosProbeSet != null) {
			return gcosProbeSet.getQCProbeSetType();
		}
		return FusionGeneChipQCProbeSetType.UnknownQCProbeSetType;
	}

	/**
	 * Gets the probe for the given index.
	 * 
	 * @param index
	 *          The zero based index to the probe (cell) vector.
	 * @param info
	 *          The probe information.
	 */
	public void getCell(final int index, final FusionCDFQCProbeInformation info) {
		info.clear();
		if (gcosProbeSet != null) {
			info.setGCOSObject(gcosProbeSet.getCell(index));
		}
	}

	/** Creates a new instance of FusionCDFQCProbeSetInformation */
	public FusionCDFQCProbeSetInformation() {
		clear();
	}

	/** Clears the members. */
	public void clear() {
		gcosProbeSet = null;
	}
}
