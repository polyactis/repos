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

/** Stores information from the header section of the CDF file. */
public class CDFFileHeader {
	/** The magic number in an XDA file */
	private int magic;

	/**
	 * Gets the magic number in an XDA file.
	 * 
	 * @return The magic number.
	 */
	public int getMagic() {
		return magic;
	}

	/**
	 * Sets the magic number in an XDA file.
	 * 
	 * @param value
	 *          The magic number.
	 */
	public void setMagic(int value) {
		magic = value;
	}

	/** The version number */
	private int version;

	/**
	 * Gets the version number in an XDA file.
	 * 
	 * @return The version number.
	 */
	public int getVersion() {
		return version;
	}

	/**
	 * Sets the version number in an XDA file.
	 * 
	 * @param value
	 *          The version number.
	 */
	public void setVersion(int value) {
		version = value;
	}

	/** The number of probe sets defined in the file */
	private int numProbeSets;

	/**
	 * Gets the number of probe sets in an XDA file.
	 * 
	 * @return The number of probe sets.
	 */
	public int getNumProbeSets() {
		return numProbeSets;
	}

	/**
	 * Sets the number of probe sets in an XDA file.
	 * 
	 * @param value
	 *          The number of probe sets.
	 */
	public void setNumProbeSets(int value) {
		numProbeSets = value;
	}

	/** The number of QC probe sets in the file */
	private int numQCProbeSets;

	/**
	 * Gets the number of QC probe sets in the file.
	 * 
	 * @return The number of QC probe sets in the file.
	 */
	public int getNumQCProbeSets() {
		return numQCProbeSets;
	}

	/**
	 * Sets the number of QC probe sets in the file.
	 * 
	 * @param value
	 *          The number of QC probe sets in the file.
	 */
	public void setNumQCProbeSets(int value) {
		numQCProbeSets = value;
	}

	/** The reference sequence (used for resequencing arrays only) */
	private String reference;

	/**
	 * Gets the reference sequence.
	 * 
	 * @return The reference sequence.
	 */
	public String getReference() {
		return reference;
	}

	/**
	 * Sets the reference sequence.
	 * 
	 * @param value
	 *          The reference sequence.
	 */
	public void setReference(String value) {
		reference = value;
	}

	/** The number of feature columns in the array */
	private int cols;

	/**
	 * Gets the number of columns on the array.
	 * 
	 * @return The number of columns.
	 */
	public int getCols() {
		return cols;
	}

	/**
	 * Sets the number of columns on the array.
	 * 
	 * @param value
	 *          The cols number.
	 */
	public void setCols(int value) {
		cols = value;
	}

	/** The number of feature rows in the array */
	private int rows;

	/**
	 * Gets the number of rows on the array.
	 * 
	 * @return The number of rows.
	 */
	public int getRows() {
		return rows;
	}

	/**
	 * Sets the number of rows on the array.
	 * 
	 * @param value
	 *          The rows number.
	 */
	public void setRows(int value) {
		rows = value;
	}

	/** Creates a new instance of CDFFileHeader */
	public CDFFileHeader() {
		cols = 0;
		rows = 0;
		magic = 0;
		version = 0;
		numProbeSets = 0;
		numQCProbeSets = 0;
		reference = null;
	}

}
