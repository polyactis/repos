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

import affymetrix.calvin.data.DataGroupHeader;
import affymetrix.calvin.data.FileHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.portability.UInt;

/** This class reads the all the DataGroupHeader information from a file into a FileHeader object. */
public class DataGroupHeaderReader {

	/** Creates a new instance of DataGroupHeaderReader */
	public DataGroupHeaderReader() {
	}

	/**
	 * Reads all the DataGroupHeaders in a file and the minimum information for each DataSetHeader in every DataGroup.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the first DataGroupHeader in the file.
	 * @param fh
	 *          FileHeader object to fill.
	 * @param dataGroupCnt
	 *          Number of DataGroup in the file.
	 */
	public void readAllMinimumInfo(FileInputStream fileStream, FileHeader fh, UInt dataGroupCnt) throws IOException,
			UnsignedOutOfLimitsException {
		// Get the first data group offset
		UInt nextDataGroupFilePos = fh.getFirstDataGroupFilePos();
		long cnt = dataGroupCnt.toLong();
		for (int i = 0; i < cnt; ++i) {
			// Read the DataGroupHeader
			DataGroupHeader dch = new DataGroupHeader();

			// Move to the indicated position in the file

			fileStream.getChannel().position(nextDataGroupFilePos.toLong());

			nextDataGroupFilePos = readMinimumInfo(fileStream, dch);
			fh.addDataGroupHdr(dch);
		}
	}

	/**
	 * Reads all the DataGroupHeaders in a file and all information for each DataSetHeader in every DataGroup.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the first DataGroupHeader in the file.
	 * @param fh
	 *          FileHeader object to fill.
	 * @param dataGroupCnt
	 *          Number of DataGroup in the file.
	 */
	public void readAll(FileInputStream fileStream, FileHeader fh, UInt dataGroupCnt) throws IOException,
			UnsignedOutOfLimitsException {
		// Get the first data group offset
		UInt nextDataGroupFilePos = fh.getFirstDataGroupFilePos();
		long cnt = dataGroupCnt.toLong();
		for (int i = 0; i < cnt; ++i) {
			// Read the DataGroupHeader
			DataGroupHeader dch = new DataGroupHeader();
			// Move to the indicated position in the file
			fileStream.getChannel().position(nextDataGroupFilePos.toLong());
			nextDataGroupFilePos = read(fileStream, dch);
			fh.addDataGroupHdr(dch);
		}
	}

	/**
	 * Reads the DataGroupHeader and the minimum information for all DataSetHeaders associated with this DataGroupHeader
	 * from the file.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the first DataGroupHeader in the file.
	 * @param dch
	 *          DataGroupHeader object to fill.
	 * @return The file position of the next data group
	 */
	public UInt readMinimumInfo(FileInputStream fileStream, DataGroupHeader dch) throws IOException,
			UnsignedOutOfLimitsException {
		UInt dataSetCnt = readHeader(fileStream, dch);
		// Read the DataSets
		DataSetHeaderReader dphReader = new DataSetHeaderReader();
		dphReader.readAllMinimumInfo(fileStream, dch, dataSetCnt);
		return dch.getNextGroupPos();
	}

	/**
	 * Read the DataGroupHeader and all DataSetHeaders associated with this DataGroupHeader from the file.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the first DataGroupHeader in the file.
	 * @param dch
	 *          DataGroupHeader object to fill.
	 * @return The file position of the next data group
	 */
	public UInt read(FileInputStream fileStream, DataGroupHeader dch) throws IOException, UnsignedOutOfLimitsException {
		UInt dataSetCnt = readHeader(fileStream, dch);

		// Read the DataSets
		DataSetHeaderReader dphReader = new DataSetHeaderReader();
		dphReader.readAll(fileStream, dch, dataSetCnt);
		return dch.getNextGroupPos();
	}

	/**
	 * Reads the DataGroupHeader from the file. Doesn't read all DataSetHeader information.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of a DataGroupHeader in the file.
	 * @param dch
	 *          DataGroupHeader object to fill with the header information.
	 * @return The number of DataSetHeaders associated with the current DataGroupHeader.
	 */
	public UInt readHeader(FileInputStream fileStream, DataGroupHeader dch) throws IOException,
			UnsignedOutOfLimitsException {
		readNextDataGroupFilePos(fileStream, dch);
		readFirstDataSetFilePos(fileStream, dch);
		UInt dataSetCnt = readDataSetCnt(fileStream, dch);
		readDataGroupName(fileStream, dch);
		return dataSetCnt;
	}

	/**
	 * reads the file position of the next DataGroup.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the file position of a DataGroupHeader.
	 * @param dch
	 *          DataGroupHeader object in which to write the file position.
	 */
	protected void readNextDataGroupFilePos(FileInputStream fileStream, DataGroupHeader dch) throws IOException,
			UnsignedOutOfLimitsException {
		dch.setNextGroupPos(FileInput.readUInt32(fileStream));
	}

	/**
	 * Reads the file position of the first DataSet associated with the current DataGroup.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the file position of the first DataSetHeader in the DataGroup.
	 * @param dch
	 *          DataGroupHeader object to which to add the DataSetHeader information.
	 */
	protected void readFirstDataSetFilePos(FileInputStream fileStream, DataGroupHeader dch) throws IOException,
			UnsignedOutOfLimitsException {
		dch.setDataSetPos(FileInput.readUInt32(fileStream));
	}

	/**
	 * Reads the number of DataSets associated with the current DataGroup.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the DataGroupHeader count.
	 * @param dch
	 *          DataGroupHeader object in which to add the DataSet count.
	 */
	protected UInt readDataSetCnt(FileInputStream fileStream, DataGroupHeader dch) throws IOException,
			UnsignedOutOfLimitsException {
		return FileInput.readUInt32(fileStream);

	}

	/**
	 * Reads the number of DataGroup name.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the DataGroupHeader name.
	 * @param dch
	 *          DataGroupHeader object to which to add the DataGroup name.
	 */
	protected void readDataGroupName(FileInputStream fileStream, DataGroupHeader dch) throws IOException {
		dch.setName(FileInput.readString16(fileStream));

	}

}
