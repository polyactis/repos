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

import affymetrix.gcos.cdf.CDFProbeSetInformation;

/** This class provides storage for probe set information from a CDF file. */
public class FusionCDFProbeSetInformation {

	/** The GCOS object. */
	private CDFProbeSetInformation probeSet;

	/**
	 * Sets the gcos object.
	 * 
	 * @param set
	 *          The gcos probe set object.
	 */
	public void setGCOSObject(final CDFProbeSetInformation set) {
		probeSet = set;
	}

	/** The number of lists (atoms) in the probe set. */
	public int getNumLists() {
		if (probeSet != null) {
			return probeSet.getNumLists();
		}
		return 0;
	}

	/** The number of groups (blocks) in the probe set. */
	public int getNumGroups() {
		if (probeSet != null) {
			return probeSet.getNumGroups();
		}
		return 0;
	}

	/** The number of probes in the set. */
	public int getNumCells() {
		if (probeSet != null) {
			return probeSet.getNumCells();
		}
		return 0;
	}

	/** An index for the probe set. */
	public int getIndex() {
		if (probeSet != null) {
			return probeSet.getIndex();
		}
		return 0;
	}

	/** An arbitrary number assigned to the probe set. */
	public int getProbeSetNumber() {
		if (probeSet != null) {
			return probeSet.getProbeSetNumber();
		}
		return 0;
	}

	/** The type of probe set. */
	public int getProbeSetType() {
		if (probeSet != null) {
			return probeSet.getProbeSetType();
		}
		return FusionGeneChipProbeSetType.UnknownProbeSetType;
	}

	/** The direction of the target that the probes are interrogating. */
	public byte getDirection() {
		if (probeSet != null) {
			return probeSet.getDirection();
		}
		return FusionDirectionType.NoDirection;
	}

	/** The number of probes per list. */
	public byte getNumCellsPerList() {
		if (probeSet != null) {
			return probeSet.getNumCellsPerList();
		}
		return 0;
	}

	/**
	 * Gets the group at the given index.
	 * 
	 * @param index
	 *          The zero based index to the group array.
	 * @param info
	 *          The group information.
	 */
	public void getGroup(final int index, final FusionCDFProbeGroupInformation info) {
		info.clear();
		if (probeSet != null) {
			info.setGCOSObject(probeSet.getGroup(index));
		}
	}

	/** Creates a new instance of FusionCDFProbeSetInformation */
	public FusionCDFProbeSetInformation() {
		clear();
	}

	/** Clears the members. */
	public void clear() {
		probeSet = null;
	}
}
