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
import java.io.FileNotFoundException;
import java.io.IOException;

import affymetrix.calvin.data.GenericData;
import affymetrix.calvin.exception.NotImplementedException;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;

/** This class reads a generic data file. It is the top-level generic data file reader. */
public class GenericFileReader {

	/* ! Hint used when opening a file */
	public enum OpenHint {
		All, Sequential, None
	};

	/* ! Indicates how much header information to read */
	public enum ReadHeaderOption {
		ReadAllHeaders, ReadMinDataGroupHeader, ReadNoDataGroupHeader
	};

	/** The name of the input file. */
	private String fileName;

	/** The file stream. */
	private FileInputStream fileStream;

	/** A pointer to the GenericData object */
	private GenericData gendata;

	/** Creates a new instance of GenericFileReader */
	public GenericFileReader() {
		fileName = null;
		fileStream = null;
		gendata = null;
	}

	/**
	 * Gets the name of the input file.
	 * 
	 * @return The name of the input file.
	 */
	public String getFilename() {
		return fileName;
	}

	/**
	 * Sets the name of the input file.
	 * 
	 * @param name
	 *          The name of the input file.
	 */
	public void setFilename(String name) {
		fileName = name;
	}

	/**
	 * Read the file header of the generic file.
	 * 
	 * @param data
	 *          A reference to a GenericData object that will receive header information from the file.
	 * @param option
	 *          Indicates how much DataGroupHeader and DataSetHeader information to read.
	 * @exception affymetrix_calvin_exceptions::FileNotFoundException
	 *              The file does not exist.
	 * @exception affymetrix_calvin_exceptions::InvalidVersionException
	 *              The file version does not match.
	 * @exception affymetrix_calvin_exceptions::InvalidFileTypeException
	 *              The file is not of the right type.
	 */
	public void readHeader(GenericData data, ReadHeaderOption option) throws FileNotFoundException,
			InvalidVersionException, InvalidFileTypeException, IOException, UnsignedOutOfLimitsException {
		openFile();
		switch (option) {
		case ReadNoDataGroupHeader:
			readFileHeaderNoDataGroupHeader(data);
			break;
		case ReadMinDataGroupHeader:
			readFileHeaderMinDP(data);
			break;
		case ReadAllHeaders: // fall through
		default:
			readFileHeader(data);
			break;
		}
		close();
	}

	public void open(GenericData data) throws FileNotFoundException, InvalidVersionException, InvalidFileTypeException,
			NotImplementedException, IOException, UnsignedOutOfLimitsException {
		open(data, OpenHint.All);
	}

	/**
	 * Open the file for reading
	 * 
	 * @param data
	 *          A reference to a GenericData object that will receive header information from the file. Amount of info
	 *          depends on the hint.
	 * @param hint
	 *          A hint on how to open the file.
	 * @exception affymetrix_calvin_exceptions::FileNotFoundException
	 *              The file does not exist.
	 * @exception affymetrix_calvin_exceptions::InvalidVersionException
	 *              The file version does not match.
	 * @exception affymetrix_calvin_exceptions::InvalidFileTypeException
	 *              The file is not of the right type.
	 */
	public void open(GenericData data, OpenHint hint) throws FileNotFoundException, InvalidVersionException,
			InvalidFileTypeException, NotImplementedException, IOException, UnsignedOutOfLimitsException {
		if (hint == OpenHint.All) {
			openFile();
			readFileHeader(data);
			gendata = data;
		}
		else {
			throw new NotImplementedException();
		}
	}

	/**
	 * Gets the number of DataGroups in the file.
	 * 
	 * @return The number of DataGroups in the file.
	 */
	public int getDataGroupCnt() {
		if (gendata != null) {
			return gendata.getDataGroupCnt();
		}
		return 0;
	}

	/**
	 * Closes the file.
	 */
	public void close() {
		if (fileStream != null) {
			try {
				fileStream.close();
			}
			catch (Throwable t) {
			}
		}
		fileStream = null;
	}

	/** Opens the file for reading */
	private void openFile() throws FileNotFoundException {
		try {
			fileStream = new FileInputStream(fileName);
		}
		catch (Throwable t) {
			throw new FileNotFoundException();
		}
	}

	/**
	 * Read the file header and minimize amount of information read from the DataSetHeaders. It does not attempt to read
	 * the complete DataSetHeader. That is deferred until accessed by the DataSet object.
	 * 
	 * @param data
	 *          Reference to the GenericData object to fill.
	 */
	private void readFileHeaderMinDP(GenericData data) throws InvalidVersionException, InvalidFileTypeException,
			IOException, UnsignedOutOfLimitsException {
		data.getHeader().setFilename(fileName);

		FileHeaderReader fhReader = new FileHeaderReader(fileStream, data.getHeader());
		fhReader.read();

		DataGroupHeaderReader dchReader = new DataGroupHeaderReader();
		dchReader.readAllMinimumInfo(fileStream, data.getHeader(), fhReader.getDataGroupCnt());

	}

	/**
	 * Reads the file header of the generic file and reads all the DataSetHeader information.
	 * 
	 * @param data
	 *          Reference to the GenericData object to fill.
	 */
	private void readFileHeader(GenericData data) throws InvalidVersionException, InvalidFileTypeException, IOException,
			UnsignedOutOfLimitsException {
		// Set the file name
		data.getHeader().setFilename(fileName);
		FileHeaderReader hdrReader = new FileHeaderReader(fileStream, data.getHeader());
		hdrReader.read();
		DataGroupHeaderReader grpReader = new DataGroupHeaderReader();
		grpReader.readAll(fileStream, data.getHeader(), hdrReader.getDataGroupCnt());
	}

	/**
	 * Reads the file header of the generic file but does not read any DataGroupHeaders or DataSetHeaders.
	 * 
	 * @param data
	 *          Reference to the GenericData object to fill.
	 */
	private void readFileHeaderNoDataGroupHeader(GenericData data) throws InvalidVersionException,
			InvalidFileTypeException, IOException, UnsignedOutOfLimitsException {
		// Set the file name
		data.getHeader().setFilename(fileName);
		FileHeaderReader hdrReader = new FileHeaderReader(fileStream, data.getHeader());
		hdrReader.read();
	}
}
