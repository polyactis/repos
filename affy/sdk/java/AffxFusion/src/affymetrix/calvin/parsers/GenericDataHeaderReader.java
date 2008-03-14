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

import java.io.FileInputStream;
import java.io.IOException;

import affymetrix.calvin.data.GenericDataHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.utils.AffymetrixGuidType;
import affymetrix.portability.UInt;

/** This class reads the GenericDataHeader from a file. */
public class GenericDataHeaderReader {

	/**
	 * Constructor
	 * 
	 * @param fis
	 *          Open stream positioned at the start of the first GenericDataHeader.
	 */
	public GenericDataHeaderReader(FileInputStream fis) {
		fileStream = fis;
	}

	/**
	 * Read the GenericDataHeader and all parent GenericDataHeaders from the input stream.
	 * 
	 * @param gdh
	 *          Reference to the GenericDataHeader object to which to add the GenericDataHeader information.
	 */
	public void read(GenericDataHeader gdh) throws IOException, UnsignedOutOfLimitsException {
		// Read data type identifier
		gdh.setFileTypeId(FileInput.readString8(fileStream));

		// Read file identifier
		AffymetrixGuidType guid = new AffymetrixGuidType();
		guid.setGuid(FileInput.readBlob(fileStream));
		gdh.setFileId(guid);

		// Read file creation time
		gdh.setFileCreationTime(FileInput.readString16(fileStream));

		// Read locale
		gdh.setLocale(FileInput.readString16(fileStream));

		// Read name value pairs
		UInt paramCount = FileInput.readUInt32(fileStream);
		for (int iparam = 0; iparam < paramCount.toLong(); ++iparam) {
			String name = FileInput.readString16(fileStream);
			byte[] mimeValue = FileInput.readBlob(fileStream);
			String type = FileInput.readString16(fileStream);
			ParameterNameValue nvt = new ParameterNameValue(name, mimeValue, type);
			gdh.addNameValParam(nvt);
		}

		// Read number of generic data parent header
		UInt numParents = FileInput.readUInt32(fileStream);

		// Read each parent header in turn - this needs to be recursive
		for (int iparent = 0; iparent < numParents.toLong(); ++iparent) {
			GenericDataHeader parentGDH = new GenericDataHeader();
			read(parentGDH);
			gdh.addParent(parentGDH);
		}
	}

	/** A reference to the file stream. */
	private FileInputStream fileStream;
}
