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

/** This class provides storage for QC probes from a CDF file. */
public class CDFQCProbeInformation {

	/** This is the size of the object in a binary CDF file. */
	public static final int QC_PROBE_SIZE = (2 + 2 + 1 + 1 + 1);

	/** The X coordinate of the probe */
	private int x;

	/** Gets the X coordinate */
	public int getX() {
		return x;
	}

	/** Sets the X coordinate */
	public void setX(int value) {
		x = value;
	}

	/** The Y coordinate */
	private int y;

	/** Gets the Y coordinate */
	public int getY() {
		return y;
	}

	/** Sets the Y coordinate */
	public void setY(int value) {
		y = value;
	}

	/** The probe length. This value may be 1 for non-synthesized features */
	private byte plen;

	/** Gets the probe length. This value may be 1 for non-synthesized features. */
	public byte getProbeLength() {
		return plen;
	}

	/** Sets the probe length. */
	public void setProbeLength(byte value) {
		plen = value;
	}

	/** Flag indicating if the probe is a perfect match probe. */
	private byte pmProbe;

	/** Gets a flag indicating if the probe is a perfect match probe. */
	public boolean isPMProbe() {
		return (pmProbe != 0);
	}

	/** Sets a flag indicating if the probe is a perfect match probe. */
	public void setPMProbe(boolean value) {
		pmProbe = (value == true ? (byte)1 : 0);
	}

	/** Flag indicating if the probe is used for background calculations (blank feature). */
	private byte background;

	/** Gets a flag indicating if the probe is used for background calculations (blank feature). */
	public boolean isBackground() {
		return (background != 0);
	}

	/** Sets a flag indicating if the probe is used for background calculations (blank feature). */
	public void setBackground(boolean value) {
		background = (value == true ? (byte)1 : 0);
	}

	/** Creates a new instance of CDFQCProbeInformation */
	public CDFQCProbeInformation() {
		x = 0;
		y = 0;
		plen = 0;
		pmProbe = 0;
		background = 0;
	}

}
