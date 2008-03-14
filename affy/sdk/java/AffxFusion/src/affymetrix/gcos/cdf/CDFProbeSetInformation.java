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

package affymetrix.gcos.cdf;

import java.nio.MappedByteBuffer;
import java.util.List;

import affymetrix.gcos.FileIO;
import affymetrix.portability.DataSizes;

/** This class provides storage for probe set information from a CDF file. */
public class CDFProbeSetInformation {

	/** The size of the probe set object in a XDA file. */
	public static final int PROBE_SET_SIZE = (4 + 4 + 4 + 4 + 2 + 1 + 1);

	/** The number of lists (atoms) in the probe set. */
	private int numLists;

	public int getNumLists() {
		return numLists;
	}

	public void setNumLists(int value) {
		numLists = value;
	}

	/** The number of groups (blocks) in the probe set. */
	private int numGroups;

	public int getNumGroups() {
		return numGroups;
	}

	public void setNumGroups(int value) {
		numGroups = value;
	}

	/** The number of probes in the set. */
	private int numCells;

	public int getNumCells() {
		return numCells;
	}

	public void setNumCells(int value) {
		numCells = value;
	}

	/** An index for the probe set. */
	private int setIndex;

	public int getIndex() {
		return setIndex;
	}

	public void setIndex(int value) {
		setIndex = value;
	}

	/** An arbitrary number assigned to the probe set. */
	private int probeSetNumber;

	public int getProbeSetNumber() {
		return probeSetNumber;
	}

	public void setProbeSetNumber(int value) {
		probeSetNumber = value;
	}

	/** The type of probe set. */
	private short probeSetType;

	public int getProbeSetType() {
		return probeSetType;
	}

	public void setProbeSetType(int value) {
		probeSetType = (short)value;
	}

	/** The direction of the target that the probes are interrogating. */
	private byte direction;

	public byte getDirection() {
		return direction;
	}

	public void setDirection(byte value) {
		direction = value;
	}

	/** The number of probes per list. */
	private byte numCellsPerList;

	public byte getNumCellsPerList() {
		return numCellsPerList;
	}

	public void setNumCellsPerList(byte value) {
		numCellsPerList = value;
	}

	/** The groups in the set. */
	private List<CDFProbeGroupInformation> groups;

	/** A mapped byte buffer for XDA files. */
	private MappedByteBuffer xdaBuffer;

	/** The offset from the map buffer to the probe set information. */
	private int offset;

	/**
	 * Gets the group at the given index.
	 * 
	 * @param index
	 *          The zero based index to the group array.
	 * @return The group information.
	 */
	public CDFProbeGroupInformation getGroup(int index) {
		if (groups != null) {
			return groups.get(index);
		}
		else if (xdaBuffer != null) {
			CDFProbeGroupInformation group = new CDFProbeGroupInformation();
			int groupOffset = offset + CDFProbeSetInformation.PROBE_SET_SIZE;
			for (int i = 0; i < index; i++) {
				int cells = FileIO.getInt32(xdaBuffer, groupOffset + DataSizes.INT_SIZE);
				groupOffset += CDFProbeGroupInformation.PROBE_GROUP_SIZE;
				groupOffset += (cells * CDFProbeInformation.PROBE_SIZE);
			}
			group.setMap(xdaBuffer, groupOffset, setIndex, index);
			return group;
		}
		return null;
	}

	public void setGroups(List<CDFProbeGroupInformation> value) {
		groups = value;
	}

	/**
	 * Sets the map and offset for memory mapping.
	 * 
	 * @param buf
	 *          The buffer.
	 * @param o
	 *          The offset to the probe set names.
	 */
	public void setMap(MappedByteBuffer buf, long o, int index) {
		xdaBuffer = buf;
		offset = (int)o;

		int setOffset = offset;

		setIndex = index;

		probeSetType = FileIO.getUInt16(xdaBuffer, setOffset);
		setOffset += DataSizes.SHORT_SIZE;

		direction = FileIO.getUInt8(xdaBuffer, setOffset);
		setOffset += DataSizes.CHAR_SIZE;

		numLists = FileIO.getInt32(xdaBuffer, setOffset);
		setOffset += DataSizes.INT_SIZE;

		numGroups = FileIO.getInt32(xdaBuffer, setOffset);
		setOffset += DataSizes.INT_SIZE;

		numCells = FileIO.getInt32(xdaBuffer, setOffset);
		setOffset += DataSizes.INT_SIZE;

		probeSetNumber = FileIO.getInt32(xdaBuffer, setOffset);
		setOffset += DataSizes.INT_SIZE;

		numCellsPerList = FileIO.getUInt8(xdaBuffer, setOffset);

	}

	/** Creates a new instance of CDFProbeSetInformation */
	public CDFProbeSetInformation() {
		groups = null;
		numLists = 0;
		numGroups = 0;
		numCells = 0;
		setIndex = 0;
		probeSetNumber = 0;
		probeSetType = GeneChipProbeSetType.UnknownProbeSetType;
		direction = DirectionType.NoDirection;
		xdaBuffer = null;
		offset = 0;
	}

}
