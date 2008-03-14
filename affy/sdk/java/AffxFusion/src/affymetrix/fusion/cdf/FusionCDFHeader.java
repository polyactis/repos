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

import affymetrix.gcos.cdf.CDFFileHeader;

/** Stores information from the header section of the CDF file. */
public class FusionCDFHeader {

	/** The gcos cdf header object. */
	private CDFFileHeader gcosHeader;

	/**
	 * Sets the GCOS CDF file header object.
	 * 
	 * @param h
	 *          The header object.
	 */
	public void setGCOSObject(final CDFFileHeader h) {
		gcosHeader = h;
	}

	/**
	 * Gets the number of probe sets in an XDA file.
	 * 
	 * @return The number of probe sets.
	 */
	public int getNumProbeSets() {
		if (gcosHeader != null) {
			return gcosHeader.getNumProbeSets();
		}
		return 0;
	}

	/**
	 * Gets the number of QC probe sets in the file.
	 * 
	 * @return The number of QC probe sets in the file.
	 */
	public int getNumQCProbeSets() {
		if (gcosHeader != null) {
			return gcosHeader.getNumQCProbeSets();
		}
		return 0;
	}

	/**
	 * Gets the reference sequence.
	 * 
	 * @return The reference sequence.
	 */
	public String getReference() {
		if (gcosHeader != null) {
			return gcosHeader.getReference();
		}
		return null;
	}

	/**
	 * Gets the number of columns on the array.
	 * 
	 * @return The number of columns.
	 */
	public int getCols() {
		if (gcosHeader != null) {
			return gcosHeader.getCols();
		}
		return 0;
	}

	/**
	 * Gets the number of rows on the array.
	 * 
	 * @return The number of rows.
	 */
	public int getRows() {
		if (gcosHeader != null) {
			return gcosHeader.getRows();
		}
		return 0;
	}

	/** Creates a new instance of FusionCDFHeader */
	public FusionCDFHeader() {
		clear();
	}

	/** Clears the members. */
	public void clear() {
		gcosHeader = null;
	}
}
