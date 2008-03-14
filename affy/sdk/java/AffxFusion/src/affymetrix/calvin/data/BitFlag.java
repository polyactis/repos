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

package affymetrix.calvin.data;

import java.util.BitSet;

/**  */
public class BitFlag {

	public static final byte DEFAULT_DATA_SET_HDR_FLAG = 0;

	private static final int FLAG_SIZE = 16;

	private BitSet flags = new BitSet(FLAG_SIZE);

	public BitFlag() {
		setFlags((short)0);
	}

	public BitFlag(short p) {
		setFlags(p);
	}

	/** Clear all flags. */
	public void clear() {
		flags.clear();
	}

	/**
	 * Get all flags.
	 * 
	 * @return number with flags set.
	 */
	public short getFlags() {
		short val = 0;
		for (int i = 0; i < FLAG_SIZE; i++) {
			if (flags.get(i) == true) {
				val += (1 << i);
			}
		}
		return val;
	}

	/**
	 * Set all flags.
	 * 
	 * @param p
	 *          number with flags set.
	 */
	public void setFlags(short p) {
		flags.clear();
		for (short i = 0; i < FLAG_SIZE / 2; i++) {
			if ((p & 1 << i) >= 1) {
				flags.set(i, true);
			}
			else {
				flags.set(i, false);
			}
		}
	}

	/**
	 * True if the default data set header flag has been set.
	 * 
	 * @return true if the default data set header flag has been set.
	 */
	public boolean hasDefaultDataSetHdr() {
		return flags.get(DEFAULT_DATA_SET_HDR_FLAG);
	}

	/**
	 * Set the default data set header flag.
	 * 
	 * @param p
	 *          true or false.
	 */
	public void setDefaultDataSetHdr(boolean p) {
		flags.set(DEFAULT_DATA_SET_HDR_FLAG, p);
	}

}
