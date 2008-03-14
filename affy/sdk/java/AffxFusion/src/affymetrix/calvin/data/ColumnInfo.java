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

import affymetrix.portability.DataSizes;

/** Base class for the varous columns. */
public class ColumnInfo {

	public enum DataSetColumnTypes {
		ByteColType, UByteColType, ShortColType, UShortColType, IntColType, UIntColType, FloatColType, ASCIICharColType, UnicodeCharColType
	}

	/** Name of the column. */
	private String name;

	/** Type of data in this column */
	private DataSetColumnTypes type;

	/** size of an individual element in bytes */
	private int size;

	/** number of elements in column */
	private int len;

	/** overhead size in bytes */
	private int overhead;

	/**
	 * Constructor - used only by dervied types
	 * 
	 * @param name_
	 *          Name of the column.
	 * @param type_
	 *          Type of data in the column.
	 * @param size_
	 *          Size of each element of in the column.
	 * @param len_
	 *          Number of elements of type in the column.
	 * @param overhead_
	 *          Number of extra bytes in the column
	 */
	protected ColumnInfo(String name_, DataSetColumnTypes type_, int size_, int len_, int overhead_) {
		name = name_;
		type = type_;
		size = size_;
		len = len_;
		overhead = overhead_;
	}

	/**
	 * Constructor - used by the file read operation
	 * 
	 * @param name_
	 *          Name of the column.
	 * @param type_
	 *          Type of data in the column.
	 * @param totalSize
	 *          Total size of the colum in bytes.
	 */
	public ColumnInfo(String name_, DataSetColumnTypes type_, int totalSize) {
		name = name_;
		type = type_;
		size = totalSize;
		len = 1;
		overhead = 0;
		if ((type == DataSetColumnTypes.ASCIICharColType) || (type == DataSetColumnTypes.UnicodeCharColType)) {
			overhead = 4;
			if (type == DataSetColumnTypes.UnicodeCharColType) {
				size = DataSizes.SHORT_SIZE;
				len = (totalSize - overhead) / size;
			}
			else if (type == DataSetColumnTypes.ASCIICharColType) {
				size = DataSizes.CHAR_SIZE;
				len = (totalSize - overhead) / size;
			}
		}
	}

	/**
	 * Equiality operator
	 * 
	 * @param p
	 *          object to compare against
	 */
	public boolean equals(ColumnInfo p) {
		return (name.equals(p.name) && (type == p.type) && (size == p.size));
	}

	/**
	 * Returns the type of the data in the column
	 * 
	 * @return Returns the type of the data in the column
	 */
	public DataSetColumnTypes getColumnType() {
		return type;
	}

	/**
	 * Returns the total size of the column in bytes.
	 * 
	 * @return Size in bytes of the column
	 */
	public int getSize() {
		return size * len + overhead;
	}

	/**
	 * Returns the number of elements of type in the column
	 * 
	 * @return Number of elements of type in the column
	 */
	int getLength() {
		return len;
	}

	/**
	 * Get the column name.
	 * 
	 * @return The column name.
	 */
	public String getName() {
		return name;
	}
}
