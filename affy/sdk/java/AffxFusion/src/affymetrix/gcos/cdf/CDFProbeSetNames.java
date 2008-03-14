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
import java.util.Vector;

import affymetrix.gcos.FileIO;

/** This class provides storage for the list of probe set names. */
public class CDFProbeSetNames {

	/** The maximum storage length for a probe set name in the CDF file. */
	public static final int MAX_PROBE_SET_NAME_LENGTH = 64;

	/** Array of probe set names, used if not memory mapping. */
	private Vector<String> probeSetNames = null;

	/** A mapped byte buffer for XDA files. */
	private MappedByteBuffer xdaBuffer;

	/** The offset from the map buffer to the probe set names. */
	private long offset;

	/**
	 * Gets the probe set name.
	 * 
	 * @param index
	 *          The index to the probe set name of interest.
	 * @return The name of the probe set.
	 */
	public String getName(int index) {
		if (probeSetNames.size() > 0) {
			return probeSetNames.get(index);
		}
		else if (xdaBuffer != null) {
			// String name = "";
			int nameOffset = index * MAX_PROBE_SET_NAME_LENGTH;
			return FileIO.getFixedString(xdaBuffer, (int)(offset + nameOffset), MAX_PROBE_SET_NAME_LENGTH);
		}
		return null;
	}

	/**
	 * Sets the size of the vector.
	 * 
	 * @param size
	 *          The size.
	 */
	public void clear() {
		probeSetNames.clear();
	}

	public void addName(String name) {
		probeSetNames.add(name);
	}

	/**
	 * Sets the name of the probe set.
	 * 
	 * @param index
	 *          The index to the probe set array/vector.
	 * @param name
	 *          The probe set name.
	 */
	public void setName(int index, String name) {
		probeSetNames.set(index, name);
	}

	public void setSize(int size) {
		probeSetNames.setSize(size);
	}

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
		offset = o;
	}

	/** Creates a new instance of CDFProbeSetNames */
	public CDFProbeSetNames() {
		xdaBuffer = null;
		offset = 0;
		probeSetNames = new Vector<String>();
	}

	public CDFProbeSetNames(int size) {
		xdaBuffer = null;
		offset = 0;
		probeSetNames = new Vector<String>(size);
	}

}
