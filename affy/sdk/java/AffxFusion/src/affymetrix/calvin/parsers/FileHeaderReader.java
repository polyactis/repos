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

import affymetrix.calvin.data.FileHeader;
import affymetrix.calvin.data.GenericDataHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.portability.UInt;

/** This class reads the FileHeader of the generic data file. */
public class FileHeaderReader {

	/**
	 * Constructor
	 * 
	 * @param fs
	 *          Open fstream positioned at the start of the file.
	 * @param fh
	 *          FileHeader object to fill.
	 */
	public FileHeaderReader(FileInputStream fs, FileHeader fh) {
		fileStream = fs;
		header = fh;
	}

	/**
	 * Reads the FileHeader.
	 * 
	 * @exception affymetrix_calvin_exceptions::InvalidVersionException
	 *              The file version does not match.
	 * @exception affymetrix_calvin_exceptions::InvalidFileTypeException
	 *              The file is not of the right type.
	 */
	public void read() throws InvalidVersionException, InvalidFileTypeException, IOException,
			UnsignedOutOfLimitsException {
		readMagicNumber();
		readVersion();
		readDataGroupCnt();
		readFirstDataGroupFilePos();
		readGenericDataHdr();
	}

	/**
	 * Gets the number of DataGroups in the file.
	 * 
	 * @return Number of DataGroups.
	 */
	public UInt getDataGroupCnt() {
		return dataGroupCnt;
	}

	/** Gets the file position of the first DataGroup header. */
	public UInt getFirstDataGroupFilePos() {
		return firstDataGroupFilePos;
	}

	/**
	 * Reads the magic number from the file.
	 * 
	 * @exception affymetrix_calvin_exceptions::InvalidFileTypeException
	 *              The file is not of the right type.
	 */
	protected void readMagicNumber() throws InvalidFileTypeException, IOException {
		// Read magic number
		byte fileMagicNumber = FileInput.readInt8(fileStream);
		if (fileMagicNumber != FileHeader.MAGIC_NUM) {
			throw new InvalidFileTypeException();
		}
	}

	/**
	 * Reads the generic file version number from the file.
	 * 
	 * @exception affymetrix_calvin_exceptions::InvalidVersionException
	 *              The file version does not match.
	 */
	protected void readVersion() throws InvalidVersionException, IOException {
		// Read generic data file version
		byte fileVersionNumber = FileInput.readInt8(fileStream);
		if (fileVersionNumber != FileHeader.VERSION) {
			throw new InvalidVersionException();
		}
	}

	/** Reads the DataGroup count from the file. */
	protected void readDataGroupCnt() throws IOException, UnsignedOutOfLimitsException {
		dataGroupCnt = FileInput.readUInt32(fileStream);
		header.setNumDataGroups(dataGroupCnt);
	}

	/** Reads the file position of the first DataGroup. */
	protected void readFirstDataGroupFilePos() throws IOException, UnsignedOutOfLimitsException {
		firstDataGroupFilePos = FileInput.readUInt32(fileStream);
		header.setFirstDataGroupFilePos(firstDataGroupFilePos);
	}

	/** Reads the GenericDataHeader from the file. */
	protected void readGenericDataHdr() throws IOException, UnsignedOutOfLimitsException {
		// Read all the GenericDataHeader
		GenericDataHeader gdh = new GenericDataHeader();
		GenericDataHeaderReader gdhReader = new GenericDataHeaderReader(fileStream);
		gdhReader.read(gdh);

		// Set the GenericDataHeader in the FileHeader
		header.setGenericDataHdr(gdh);
	}

	/** A reference to the file stream. */
	protected FileInputStream fileStream;

	/** FileHeader object to fill. */
	protected FileHeader header;

	/** Number of DataGroups. */
	UInt dataGroupCnt = UInt.ZERO;

	/** Position of the first DataGroup. */
	UInt firstDataGroupFilePos = UInt.ZERO;
}
