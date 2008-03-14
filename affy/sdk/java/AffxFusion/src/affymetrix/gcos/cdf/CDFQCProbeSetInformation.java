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

/** This class provides storage for the probes in a QC probe set. */
public class CDFQCProbeSetInformation {

	/** This is the size of the QC probe set object in a binary CDF file. */
	public static final int QC_PROBE_SET_SIZE = (4 + 2);

	/** The number of probes in the set. */
	private int numCells;

	public int getNumCells() {
		return numCells;
	}

	public void setNumCells(int value) {
		numCells = value;
	}

	/** The type of QC probes. */
	private short qcProbeSetType;

	public short getQCProbeSetType() {
		return qcProbeSetType;
	}

	public void setQCProbeSetType(short value) {
		qcProbeSetType = value;
	}

	/** The array of probes. */
	private List<CDFQCProbeInformation> cells;

	public CDFQCProbeInformation getCell(int index) {
		if (cells != null) {
			return cells.get(index);
		}
		else if (xdaBuffer != null) {
			CDFQCProbeInformation info = new CDFQCProbeInformation();
			int cellOffset = offset + QC_PROBE_SET_SIZE + (index * CDFQCProbeInformation.QC_PROBE_SIZE);

			info.setX(FileIO.getUInt16(xdaBuffer, cellOffset));
			cellOffset += DataSizes.SHORT_SIZE;

			info.setY(FileIO.getUInt16(xdaBuffer, cellOffset));
			cellOffset += DataSizes.SHORT_SIZE;

			info.setProbeLength(FileIO.getUInt8(xdaBuffer, cellOffset));
			cellOffset += DataSizes.CHAR_SIZE;

			info.setPMProbe(FileIO.getUInt8(xdaBuffer, cellOffset) == 1 ? true : false);
			cellOffset += DataSizes.CHAR_SIZE;

			info.setBackground(FileIO.getUInt8(xdaBuffer, cellOffset) == 1 ? true : false);
			return info;
		}
		return null;
	}

	public void setCells(List<CDFQCProbeInformation> value) {
		cells = value;
	}

	/** A mapped byte buffer for XDA files. */
	private MappedByteBuffer xdaBuffer;

	/** The offset from the map buffer to the probe set information. */
	private int offset;

	/**
	 * Sets the map and offset for memory mapping.
	 * 
	 * @param buf
	 *          The buffer.
	 * @param o
	 *          The offset to the probe set names.
	 */
	public void setMap(MappedByteBuffer buf, long o) {
		xdaBuffer = buf;
		offset = (int)o;
		qcProbeSetType = FileIO.getUInt16(xdaBuffer, offset);
		numCells = FileIO.getInt32(xdaBuffer, offset + DataSizes.SHORT_SIZE);
	}

	/** Creates a new instance of CDFQCProbeSetInformation */
	public CDFQCProbeSetInformation() {
		cells = null;
		xdaBuffer = null;
		offset = 0;
		numCells = 0;
		qcProbeSetType = GeneChipQCProbeSetType.UnknownQCProbeSetType;
	}

}
