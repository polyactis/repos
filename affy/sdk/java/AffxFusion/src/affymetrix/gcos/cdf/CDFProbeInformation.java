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

/** This class provides storage for an individual probe in a CDF file. */
public class CDFProbeInformation {

	/** The size of the probe object as stored in the XDA CDF file. */
	public static final int PROBE_SIZE = (4 + 4 + 2 + 2 + 1 + 1);

	/** The index of the probes probe pair or quartet (probe list). Also known as the atom position. */
	private int listIndex;

	/** Gets the index of the probes probe pair or quartet (probe list). Also known as the atom position. */
	public int getListIndex() {
		return listIndex;
	}

	/** Sets the index of the probes probe pair or quartet (probe list). Also known as the atom position. */
	public void setListIndex(int value) {
		listIndex = value;
	}

	/**
	 * The expos value in the CDF file, this can either be a zero based index equal to the list index or the exon
	 * position.
	 */
	private int expos;

	/**
	 * Gets the expos value in the CDF file, this can either be a zero based index equal to the list index or the exon
	 * position.
	 */
	public int getExpos() {
		return expos;
	}

	/**
	 * Sets the expos value in the CDF file, this can either be a zero based index equal to the list index or the exon
	 * position.
	 */
	public void setExpos(int value) {
		expos = value;
	}

	/** The X coordinate */
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

	/** The probes base at the interrogation position */
	private char pbase;

	/** Gets the probes base at the interrogation position */
	public char getPBase() {
		return pbase;
	}

	/** Sets the probes base at the interrogation position */
	public void setPBase(char value) {
		pbase = value;
	}

	/** The targets base at the interrogation position */
	private char tbase;

	/** Gets the targets base at the interrogation position */
	public char getTBase() {
		return tbase;
	}

	/** Sets the targets base at the interrogation position */
	public void setTBase(char value) {
		tbase = value;
	}

	/** Creates a new instance of CDFProbeInformation */
	public CDFProbeInformation() {
		listIndex = 0;
		expos = 0;
		x = 0;
		y = 0;
		pbase = ' ';
		tbase = ' ';
	}

}
