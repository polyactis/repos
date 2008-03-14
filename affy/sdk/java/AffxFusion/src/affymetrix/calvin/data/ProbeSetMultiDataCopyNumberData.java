////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
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
////////////////////////////////////////////////////////////////

package affymetrix.calvin.data;

/** Stores data for a copy number result. */
public class ProbeSetMultiDataCopyNumberData extends ProbeSetMultiDataBase {

	/* A value to represent the X chromosome. */
	public static final byte X_CHR = 24;

	/* A value to represent the Y chromosome. */
	public static final byte Y_CHR = 25;

	/* A value to represent the MT chromosome. */
	public static final byte MT_CHR = 26;

	/* A value to represent the no chromosome. */
	public static final byte NO_CHR = -128;

	/* The name of the probe set. */
	private String name;

	/* The chromosome. */
	private byte chr;

	/* The genomic position. */
	private int position;

	/* Creates a new instance of ProbeSetMultiDataCopyNumberData */
	public ProbeSetMultiDataCopyNumberData() {
		name = null;
		chr = 0;
		position = 0;
	}

	/* The name of the probe set. */
	public String getName() {
		return name;
	}

	/* Set the name */
	public void setName(String n) {
		name = n;
	}

	/* Get the chromosome. */
	public byte getChr() {
		return chr;
	}

	/* Set the chromosome. */
	public void setChr(byte c) {
		chr = c;
	}

	/* Get the position. */
	public int getPosition() {
		return position;
	}

	/* Set the position */
	public void setPosition(int p) {
		position = p;
	}

	/*
	 * Convert a string representation of a chromosome to a numeric representation.
	 * 
	 * @param chr The chromosome value. @return A numeric representation of the chromosome value.
	 */
	public static byte chromosomeFromString(String chr) {
		int chrValue = 0;
		try {
			chrValue = Integer.parseInt(chr);
		}
		catch (Exception e) {
			chrValue = 0;
		}
		if (chrValue == 0) {
			if ((chr == "X") || (chr == "x")) {
				chrValue = X_CHR;
			}
			else if ((chr == "Y") || (chr == "y")) {
				chrValue = Y_CHR;
			}
			else if ((chr == "MT") || (chr == "mt")) {
				chrValue = MT_CHR;
			}
			else {
				chrValue = NO_CHR;
			}
		}
		return (byte)chrValue;
	}

	/*
	 * Convert a numeric representation of a chromosome to a string representation.
	 * 
	 * @param chr The chromosome value. @return A string representation of the chromosome value.
	 */
	public static String chromosomeToString(byte chr) {
		if (chr == X_CHR) {
			return "X";
		}
		else if (chr == Y_CHR) {
			return "Y";
		}
		else if (chr == MT_CHR) {
			return "MT";
		}
		else if (chr == NO_CHR) {
			return "-";
		}
		else {
			return Integer.toString(chr);
		}
	}

}
