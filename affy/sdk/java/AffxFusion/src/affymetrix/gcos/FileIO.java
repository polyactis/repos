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

package affymetrix.gcos;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;

import affymetrix.portability.DataSizes;

/**
 * 
 * @author ljevon
 */
public class FileIO {

	/** Creates a new instance of FileIO */
	public FileIO() {
	}

	/**
	 * Reads an integer from a little endian file.
	 * 
	 * @param fis
	 *          The file stream.
	 * @return The integer.
	 */
	public static int readInt32(FileInputStream fis) {
		byte[] b = new byte[DataSizes.INT_SIZE];
		try {
			fis.read(b);
			return ByteBuffer.wrap(b).order(ByteOrder.LITTLE_ENDIAN).getInt();
		}
		catch (Throwable t) {
			return 0;
		}
	}

	/**
	 * Reads a float from a little endian file.
	 * 
	 * @param fis
	 *          The file stream.
	 * @return The float.
	 */
	public static float readFloat(FileInputStream fis) {
		byte[] b = new byte[DataSizes.FLOAT_SIZE];
		try {
			fis.read(b);
			return ByteBuffer.wrap(b).order(ByteOrder.LITTLE_ENDIAN).getFloat();
		}
		catch (Throwable t) {
			return 0.0f;
		}
	}

	/**
	 * Reads a short from a little endian file.
	 * 
	 * @param fis
	 *          The file stream.
	 * @return The short.
	 */
	public static short readInt16(FileInputStream fis) {
		byte[] b = new byte[DataSizes.SHORT_SIZE];
		try {
			fis.read(b);
			return ByteBuffer.wrap(b).order(ByteOrder.LITTLE_ENDIAN).getShort();
		}
		catch (Throwable t) {
			return 0;
		}
	}

	/**
	 * Reads an unsigned short from a little endian file.
	 * 
	 * @param fis
	 *          The file stream.
	 * @return The value as an integer.
	 */
	public static int readUInt16(FileInputStream fis) {
		byte[] b = new byte[DataSizes.SHORT_SIZE];
		try {
			fis.read(b);
			return ByteBuffer.wrap(b).order(ByteOrder.LITTLE_ENDIAN).getShort();
		}
		catch (Throwable t) {
			return 0;
		}
	}

	/**
	 * Reads a byte from a little endian file.
	 * 
	 * @param fis
	 *          The file stream.
	 * @return The byte.
	 */
	public static byte readInt8(FileInputStream fis) {
		try {
			return (byte)fis.read();
		}
		catch (Throwable t) {
			return 0;
		}
	}

	/**
	 * Reads a string from a little endian file.
	 * 
	 * @param fis
	 *          The file stream.
	 * @return The string.
	 */
	public static String readString(FileInputStream fis) {
		try {
			int n = readInt32(fis);
			if (n > 0) {
				return readFixedString(fis, n);
			}
		}
		catch (Throwable t) {
		}
		return null;
	}

	/**
	 * Read a string of fixed length.
	 * 
	 * @param fis
	 *          The file stream.
	 * @param len
	 *          The length of the string.
	 * @return The string.
	 */
	public static String readFixedString(FileInputStream fis, int len) {
		try {
			byte[] b = new byte[len];
			fis.read(b);
			String str = new String();
			for (int i = 0; i < len; i++) {
				if (b[i] == 0) {
					break;
				}
				str += (char)b[i];
			}
			return str;
		}
		catch (Throwable t) {
		}
		return null;
	}

	/**
	 * Reads the next line of text.
	 * 
	 * @param b
	 *          The buffered reader.
	 * @return The next line from the file or null if not found.
	 */
	public static String readNextLine(BufferedReader b) {
		String str;
		try {
			while ((str = b.readLine()) != null) {
				if (str.length() > 0) {
					return str;
				}
			}
		}
		catch (Throwable t) {
		}
		return null;
	}

	/**
	 * Reads the next floating point value from the buffer.
	 * 
	 * @param buffer
	 *          The buffered reader.
	 * @param offset
	 *          The offset to the file.
	 * @return The next floating point value.
	 */
	public static float getFloat(MappedByteBuffer buffer, int offset) {
		ByteBuffer b = (ByteBuffer)buffer.position(offset);
		return b.getFloat();
	}

	/**
	 * Reads the next integer from the buffer.
	 * 
	 * @param buffer
	 *          The buffered reader.
	 * @param offset
	 *          The offset to the file.
	 * @return The next integer value.
	 */
	public static int getInt32(MappedByteBuffer buffer, int offset) {
		ByteBuffer b = (ByteBuffer)buffer.position(offset);
		return b.getInt();
	}

	/**
	 * Reads the next short from the buffer.
	 * 
	 * @param buffer
	 *          The buffered reader.
	 * @param offset
	 *          The offset to the file.
	 * @return The next short value.
	 */
	public static short getUInt16(MappedByteBuffer buffer, int offset) {
		ByteBuffer b = (ByteBuffer)buffer.position(offset);
		return b.getShort();
	}

	/**
	 * Reads the next short from the buffer.
	 * 
	 * @param buffer
	 *          The buffered reader.
	 * @param offset
	 *          The offset to the file.
	 * @return The next short value.
	 */
	public static short getInt16(MappedByteBuffer buffer, int offset) {
		ByteBuffer b = (ByteBuffer)buffer.position(offset);
		return b.getShort();
	}

	/**
	 * Reads the next byte from the buffer.
	 * 
	 * @param buffer
	 *          The buffered reader.
	 * @param offset
	 *          The offset to the file.
	 * @return The next byte value.
	 */
	public static byte getUInt8(MappedByteBuffer buffer, int offset) {
		ByteBuffer b = (ByteBuffer)buffer.position(offset);
		return b.get();
	}

	/**
	 * Reads the next byte from the buffer.
	 * 
	 * @param buffer
	 *          The buffered reader.
	 * @param offset
	 *          The offset to the file.
	 * @return The next byte value.
	 */
	public static byte getInt8(MappedByteBuffer buffer, int offset) {
		ByteBuffer b = (ByteBuffer)buffer.position(offset);
		return b.get();
	}

	/**
	 * Read a string of fixed length.
	 * 
	 * @param buffer
	 *          The buffered reader.
	 * @param offset
	 *          The offset to the file.
	 * @param len
	 *          The string length.
	 * @return The string.
	 */
	public static String getFixedString(MappedByteBuffer buffer, int offset, int len) {
		try {
			String str = "";
			byte bvalue;
			ByteBuffer b = (ByteBuffer)buffer.position(offset);
			for (int i = 0; i < len; i++) {
				bvalue = b.get();
				if (bvalue == 0) {
					break;
				}
				str += (char)bvalue;
			}
			return str;
		}
		catch (Throwable t) {
			return "";
		}
	}

}
