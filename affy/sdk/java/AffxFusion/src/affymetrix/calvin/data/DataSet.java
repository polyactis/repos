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

package affymetrix.calvin.data;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel.MapMode;
import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parsers.FileInput;
import affymetrix.portability.UByte;
import affymetrix.portability.UInt;
import affymetrix.portability.UShort;

/** This class provides methods to access the data of a DataSet. */
public class DataSet {

	/** name of the file containing the data dataGroup. */
	// private String fileName;
	/** copy of the DataSetHeader */
	private DataSetHeader header;

	/** pointer to the mapped data, doesn't account for allocation granularity. */
	private MappedByteBuffer mappedData;

	/** The file input stream. */
	private FileInputStream fileInputStream;

	/**
	 * Array of column byte offsets. Updated when the file is opened. There are columns + 1 elements
	 */
	private List<Integer> columnByteOffsets = new ArrayList<Integer>();

	/** Byte offset to the start of the view */
	private UInt mapStart = UInt.ZERO;

	/** Number of bytes mapped to the view */
	private UInt mapLen = new UInt();

	/**
	 * Constructor. Use this ctor to access the data using memory-mapping. On Windows, memory-mapping will be restricted
	 * to 200MB view of the DataSet data.
	 * 
	 * @param fileName_
	 *          The name of the generic file to access.
	 * @param header_
	 *          The DataSetHeader of the DataSet to access.
	 * @param handle_
	 *          A handle to the file mapping object
	 */
	public DataSet(String fileName_, DataSetHeader header_, FileInputStream handle_) {
		// fileName = fileName_;
		header = header_;
		mappedData = null;
		fileInputStream = handle_;
	}

	/**
	 * Method to release memory held by this object. Closes object before deleting.
	 */
	public void delete() {
		close();
	}

	/**
	 * Method to open the DataSet to access the data.
	 * 
	 * @return true if successful
	 */
	public boolean open() throws IOException, UnsignedOutOfLimitsException {
		updateColumnByteOffsets();
		return openMM();
	}

	/**
	 * Open the DataSet using memory-mapping
	 * 
	 * @return True if the DataSet was successully mapped.
	 */
	private boolean openMM() throws IOException, UnsignedOutOfLimitsException {
		mapLen.set(header.getDataSize());
		mapStart = header.getDataStartFilePos();
		mappedData = fileInputStream.getChannel().map(MapMode.READ_ONLY, mapStart.toLong(), mapLen.toLong());
		mappedData.order(ByteOrder.BIG_ENDIAN);
		return true;
	}

	/** Updates the columnByteOffsets member. */
	private void updateColumnByteOffsets() {
		columnByteOffsets.clear();
		int accum = 0;
		int cols = header.getColumnCnt();
		for (int col = 0; col < cols; ++col) {
			columnByteOffsets.add(accum);
			accum += header.getColumnInfo(col).getSize();
		}
		columnByteOffsets.add(accum);
	}

	/** Method to close the DataSet. */
	public void close() {
		mappedData = null;
		fileInputStream = null;
	}

	/**
	 * Method to get a reference to the DataSetHeader
	 * 
	 * @return A reference to the DataSetHeader.
	 */
	public DataSetHeader getHeader() {
		return header;
	}

	/** Return the number of rows in the DataSet. */
	public int getRows() {
		return header.getRowCnt();
	}

	/** Return the number of columns in the DataSet. */
	public int getCols() {
		return header.getColumnCnt();
	}

	/**
	 * Determines if the DataSet is open
	 * 
	 * @return true if the DataSet is open
	 */
	boolean isOpen() {
		return (mappedData != null);
	}

	/**
	 * Provides access to single data elements
	 * 
	 * @param row
	 *          Row index.
	 * @param col
	 *          Column index.
	 * @return Reference to the data type to fill with the data.
	 */
	public byte getDataByte(int row, int col) throws UnsignedOutOfLimitsException {
		if (header.getRowCnt() == 0) {
			return 0;
		}
		mappedData.position((int)filePosition(row, col).toLong());
		return FileInput.readInt8(mappedData);
	}

	public UByte getDataUByte(int row, int col) throws UnsignedOutOfLimitsException {
		if (header.getRowCnt() == 0) {
			return UByte.ZERO;
		}
		mappedData.position((int)filePosition(row, col).toLong());
		return FileInput.readUInt8(mappedData);
	}

	public short getDataShort(int row, int col) throws UnsignedOutOfLimitsException {
		if (header.getRowCnt() == 0) {
			return 0;
		}
		mappedData.position((int)filePosition(row, col).toLong());
		return FileInput.readInt16(mappedData);
	}

	public UShort getDataUShort(int row, int col) throws UnsignedOutOfLimitsException {
		if (header.getRowCnt() == 0) {
			return UShort.ZERO;
		}
		mappedData.position((int)filePosition(row, col).toLong());
		return FileInput.readUInt16(mappedData);
	}

	public int getDataInt(int row, int col) throws UnsignedOutOfLimitsException {
		if (header.getRowCnt() == 0) {
			return 0;
		}
		mappedData.position((int)filePosition(row, col).toLong());
		return FileInput.readInt32(mappedData);
	}

	public UInt getDataUInt(int row, int col) throws UnsignedOutOfLimitsException {
		if (header.getRowCnt() == 0) {
			return UInt.ZERO;
		}
		mappedData.position((int)filePosition(row, col).toLong());
		return FileInput.readUInt32(mappedData);
	}

	public float getDataFloat(int row, int col) throws UnsignedOutOfLimitsException {
		if (header.getRowCnt() == 0) {
			return 0.0f;
		}
		mappedData.position((int)filePosition(row, col).toLong());
		return FileInput.readFloat(mappedData);
	}

	public String getDataString16(int row, int col) throws IOException, UnsignedOutOfLimitsException {
		if (header.getRowCnt() == 0) {
			return null;
		}
		mappedData.position((int)filePosition(row, col).toLong());
		return FileInput.readString16(mappedData);
	}

	public String getDataString8(int row, int col) throws IOException, UnsignedOutOfLimitsException {
		if (header.getRowCnt() == 0) {
			return null;
		}
		mappedData.position((int)filePosition(row, col).toLong());
		return FileInput.readString8(mappedData);
	}

	/**
	 * Return the bytes per row.
	 * 
	 * @return Bytes in a row.
	 */
	private int bytesPerRow() {
		return columnByteOffsets.get(header.getColumnCnt());
	}

	public UInt filePosition(int rowStart, int col) throws UnsignedOutOfLimitsException {
		return filePosition(rowStart, col, 1);
	}

	/**
	 * Returns the address of a data element given a row and column. Ensures that data from rowStart to rowCount+rowStart
	 * are mapped unless that is larger than the mapped window.
	 * 
	 * @param rowStart
	 *          Row index
	 * @param col
	 *          Column index
	 * @param rowCount
	 *          The number of rows to ensure are mapped starting at rowStart
	 * @return Pointer to the data element at rowStart
	 */
	public UInt filePosition(int rowStart, int col, int rowCount) throws UnsignedOutOfLimitsException {
		UInt startByte = new UInt(header.getDataStartFilePos());
		startByte.add((bytesPerRow() * rowStart) + columnByteOffsets.get(col));
		startByte.subtract(mapStart);
		return startByte;
	}
}
