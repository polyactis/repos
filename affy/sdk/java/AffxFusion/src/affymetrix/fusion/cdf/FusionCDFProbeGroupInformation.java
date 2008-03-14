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

import affymetrix.gcos.cdf.CDFProbeGroupInformation;

/** This class provides storage for a group of probes, also known as a block. */
public class FusionCDFProbeGroupInformation {

	/** The GCOS probe group object. */
	private CDFProbeGroupInformation gcosGroup;

	/**
	 * Sets the GCOS probe group object.
	 * 
	 * @param g
	 *          The probe group object.
	 */
	public void setGCOSObject(final CDFProbeGroupInformation g) {
		gcosGroup = g;
	}

	/** The number of probe pairs or probe quartets in the group. */
	public int getNumLists() {
		if (gcosGroup != null) {
			return gcosGroup.getNumLists();
		}
		return 0;
	}

	/** The number of probes in the group. */
	public int getNumCells() {
		if (gcosGroup != null) {
			return gcosGroup.getNumCells();
		}
		return 0;
	}

	/** The first probes list index value. */
	public int getStart() {
		if (gcosGroup != null) {
			return gcosGroup.getStart();
		}
		return 0;
	}

	/** The last probes list index value. */
	public int getStop() {
		if (gcosGroup != null) {
			return gcosGroup.getStop();
		}
		return 0;
	}

	/** The probe set index associated with the group. */
	public int getProbeSetIndex() {
		if (gcosGroup != null) {
			return gcosGroup.getProbeSetIndex();
		}
		return 0;
	}

	/** The group index. */
	public int getGroupIndex() {
		if (gcosGroup != null) {
			return gcosGroup.getGroupIndex();
		}
		return 0;
	}

	/** The name of the group. */
	public String getName() {
		if (gcosGroup != null) {
			return gcosGroup.getName();
		}
		return null;
	}

	/** The number of cells per list (2 for expression and genotyping, 4 for resquencing). */
	public byte getNumCellsPerList() {
		if (gcosGroup != null) {
			return gcosGroup.getNumCellsPerList();
		}
		return 0;
	}

	/** The direction of the target that the probes are interrogating. */
	public byte getDirection() {
		if (gcosGroup != null) {
			return gcosGroup.getDirection();
		}
		return FusionDirectionType.NoDirection;
	}

	/**
	 * Gets the cell information for the given cell index.
	 * 
	 * @param index
	 *          The zero based index to the cell array.
	 * @param cel
	 *          The cell or probe information.
	 */
	public void getCell(final int index, final FusionCDFProbeInformation cel) {
		cel.clear();
		if (gcosGroup != null) {
			cel.setGCOSObject(gcosGroup.getCell(index));
		}
	}

	/** Creates a new instance of FusionCDFProbeGroupInformation */
	public FusionCDFProbeGroupInformation() {
		clear();
	}

	/** Clears the members. */
	public void clear() {
		gcosGroup = null;
	}
}
