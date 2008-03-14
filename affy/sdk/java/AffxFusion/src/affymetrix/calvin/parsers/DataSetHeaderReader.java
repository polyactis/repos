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

import affymetrix.calvin.data.ColumnInfo;
import affymetrix.calvin.data.DataGroupHeader;
import affymetrix.calvin.data.DataSetHeader;
import affymetrix.calvin.data.ColumnInfo.DataSetColumnTypes;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.portability.UInt;

/** This class reads the all DataSetHeader information associated with a DataGroup from a file. */
public class DataSetHeaderReader {

	/** Creates a new instance of DataSetHeaderReader */
	public DataSetHeaderReader() {
	}

	/**
	 * Reads the minimum DataSetHeader information for all DataSets associated with a DataGroup.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the first DataSetHeader in a DataGroup.
	 * @param dch
	 *          DataGroupHeader object to which to add the DataSetHeader information.
	 * @param dataSetCnt
	 *          Number of DataSets in the DataGroup.
	 */
	public void readAllMinimumInfo(FileInputStream fileStream, DataGroupHeader dch, UInt dataSetCnt) throws IOException,
			UnsignedOutOfLimitsException {
		// Get the first dataSet offset
		UInt nextDataSetFilePos = dch.getDataSetPos();
		long cnt = dataSetCnt.toLong();
		for (int i = 0; i < cnt; ++i) {
			DataSetHeader dph = new DataSetHeader();

			// Move to the indicated position in the file
			fileStream.getChannel().position(nextDataSetFilePos.toLong());
			nextDataSetFilePos = readMinimumInfo(fileStream, dph);

			// Add the DataSetHeader to the file header
			dch.addDataSetHdr(dph);
		}
	}

	/**
	 * Reads the complete DataSetHeader information for all DataSets associated with a DataGroup.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the first DataSetHeader in a DataGroup.
	 * @param dch
	 *          DataGroupHeader object to which to add the DataSetHeader information.
	 * @param dataSetCnt
	 *          Number of DataSets in the DataGroup.
	 */
	public void readAll(FileInputStream fileStream, DataGroupHeader dch, UInt dataSetCnt) throws IOException,
			UnsignedOutOfLimitsException {

		// Get the first dataSet offset
		UInt nextDataSetFilePos = dch.getDataSetPos();

		long cnt = dataSetCnt.toLong();
		for (int i = 0; i < cnt; ++i) {
			DataSetHeader dph = new DataSetHeader();

			// Move to the indicated position in the file
			fileStream.getChannel().position(nextDataSetFilePos.toLong());

			nextDataSetFilePos = read(fileStream, dph);

			// Add the DataSetHeader to the file header
			dch.addDataSetHdr(dph);
		}
	}

	/**
	 * Reads the minimum DataSetHeader information.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the DataSetHeader.
	 * @param dsh
	 *          Reference to the DataSetHeader object to fill.
	 * @return The file position of the next DataSet.
	 */
	public UInt readMinimumInfo(FileInputStream fileStream, DataSetHeader dsh) throws IOException,
			UnsignedOutOfLimitsException {
		readDataSetStartFilePos(fileStream, dsh);
		readDataFilePos(fileStream, dsh);
		UInt nextDataSetFilePos = readNextDataSetFilePos(fileStream, dsh);
		readName(fileStream, dsh);
		return nextDataSetFilePos;
	}

	/**
	 * Reads the complete DataSetHeader information.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the DataSetHeader.
	 * @param dsh
	 *          Reference to the DataSetHeader object to fill.
	 * @return The file position of the next DataSet.
	 */
	public UInt read(FileInputStream fileStream, DataSetHeader dsh) throws IOException, UnsignedOutOfLimitsException {
		readDataSetStartFilePos(fileStream, dsh);
		readDataFilePos(fileStream, dsh);
		UInt nextDataSetFilePos = readNextDataSetFilePos(fileStream, dsh);
		readName(fileStream, dsh);
		readParameters(fileStream, dsh);
		readColumns(fileStream, dsh);
		readRowCount(fileStream, dsh);
		return nextDataSetFilePos;

	}

	/**
	 * Read the file position of the start of the DataSet.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the DataSetHeader.
	 * @param dsh
	 *          Reference to the DataSetHeader object to fill.
	 */
	protected void readDataSetStartFilePos(FileInputStream fileStream, DataSetHeader dsh) {
		try {
			long pos = fileStream.getChannel().position();
			dsh.setHeaderStartFilePos(new UInt(pos));
		}
		catch (Throwable t) {
		}
	}

	/**
	 * Read the file position to the start of the data.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the data file position.
	 * @param dsh
	 *          Reference to the DataSetHeader object to fill.
	 */
	protected void readDataFilePos(FileInputStream fileStream, DataSetHeader dsh) throws IOException,
			UnsignedOutOfLimitsException {
		dsh.setDataStartFilePos(FileInput.readUInt32(fileStream));
	}

	/**
	 * Read the file position to the next DataSet.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the next DataSet file position.
	 * @param dsh
	 *          Reference to the DataSetHeader object to fill.
	 * @return The file position of the next data set.
	 */
	protected UInt readNextDataSetFilePos(FileInputStream fileStream, DataSetHeader dsh) throws IOException,
			UnsignedOutOfLimitsException {
		UInt nextDataSetFilePos = FileInput.readUInt32(fileStream);
		dsh.setNextSetFilePos(nextDataSetFilePos);
		return nextDataSetFilePos;
	}

	/**
	 * Read the DataSetHeader name.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the DataSetHeader name.
	 * @param dsh
	 *          Reference to the DataSetHeader object to fill.
	 */
	protected void readName(FileInputStream fileStream, DataSetHeader dsh) throws IOException,
			UnsignedOutOfLimitsException {
		String name = FileInput.readString16(fileStream);
		dsh.setName(name);
	}

	/**
	 * Read the parameter list (name-value-type).
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the DataSetHeader parameter list count.
	 * @param dsh
	 *          Reference to the DataSetHeader object to fill.
	 */
	protected void readParameters(FileInputStream fileStream, DataSetHeader dsh) throws IOException,
			UnsignedOutOfLimitsException {
		UInt params = FileInput.readUInt32(fileStream);
		long paramsCnt = params.toLong();
		for (int i = 0; i < paramsCnt; ++i) {
			String paramName = FileInput.readString16(fileStream);
			byte[] mimeValue = FileInput.readBlob(fileStream);
			String paramType = FileInput.readString16(fileStream);
			ParameterNameValue nvt = new ParameterNameValue(paramName, mimeValue, paramType);
			dsh.addNameValParam(nvt);
		}
	}

	/**
	 * Read column information.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the DataSetHeader column count.
	 * @param dsh
	 *          Reference to the DataSetHeader object to fill.
	 */
	protected void readColumns(FileInputStream fileStream, DataSetHeader dsh) throws IOException,
			UnsignedOutOfLimitsException {
		// Read the number of columns
		UInt columns = FileInput.readUInt32(fileStream);
		long columnCnt = columns.toLong();
		for (int icol = 0; icol < columnCnt; ++icol) {
			// Read the name
			String name = FileInput.readString16(fileStream);
			// Read the type
			int type = FileInput.readInt8(fileStream);

			// Read the size
			int size = FileInput.readInt32(fileStream);

			dsh.addColumn(new ColumnInfo(name, DataSetColumnTypes.values()[type], size));
		}
	}

	/**
	 * Read the number of rows.
	 * 
	 * @param fileStream
	 *          Open fstream positioned at the start of the DataSetHeader row count.
	 * @param dsh
	 *          Reference to the DataSetHeader object to fill.
	 */
	protected void readRowCount(FileInputStream fileStream, DataSetHeader dsh) throws IOException,
			UnsignedOutOfLimitsException {
		int numRows = FileInput.readInt32(fileStream);
		dsh.setRowCnt(numRows);
	}

}
