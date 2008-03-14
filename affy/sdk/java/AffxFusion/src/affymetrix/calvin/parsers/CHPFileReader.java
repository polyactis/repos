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

package affymetrix.calvin.parsers;

import java.io.FileNotFoundException;
import java.io.IOException;

import affymetrix.calvin.data.CHPData;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;

/** This class reads a CHP data file. It is an interpreter class. */
public class CHPFileReader {

	/** Creates a new instance of CHPFileReader */
	public CHPFileReader() {
		fileName = null;
	}

	/**
	 * Gets the name of the input file.
	 * 
	 * @return The name of the file to read.
	 */
	public String getFilename() {
		return fileName;
	}

	/**
	 * Sets the name of the input file.
	 * 
	 * @param name
	 *          The name of the file to read.
	 */
	public void setFilename(String name) {
		fileName = name;
	}

	/**
	 * Reads the file header of the generic file and reads all the DataPlaneHeader information.
	 * 
	 * @param data
	 *          A reference to a GenericData object that will receive information from the file.
	 * @exception affymetrix_calvin_exceptions::FileNotFoundException
	 *              The file does not exist.
	 * @exception affymetrix_calvin_exceptions::InvalidVersionException
	 *              The file version does not match.
	 * @exception affymetrix_calvin_exceptions::InvalidFileTypeException
	 *              The file is not of the right type.
	 */
	public void read(CHPData data) throws FileNotFoundException, InvalidVersionException, InvalidFileTypeException,
			IOException, UnsignedOutOfLimitsException {
		data.clear();
		GenericFileReader reader = new GenericFileReader();
		if (fileName == null) {
			fileName = data.getFilename();
		}
		reader.setFilename(fileName);
		reader.readHeader(data.getGenericData(), GenericFileReader.ReadHeaderOption.ReadAllHeaders);
	}

	/** Name of the file to read */
	private String fileName;
}
