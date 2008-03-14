/////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
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
package affymetrix.calvin.writers;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import affymetrix.calvin.utils.IOUtils;
import affymetrix.portability.DataSizes;
import affymetrix.portability.UByte;
import affymetrix.portability.UInt;
import affymetrix.portability.UShort;

public class FileOutput {

	// / <summary>
	// / Write an 8 bit string of fixed size to a big endian file.
	// / </summary>
	// / <param name="input">The out file stream.</param>
	// / <param name="length">The length of the string.</param>
	public static void writeString8(BufferedFileOutputStream output, String value, int maxLnChar) throws IOException {
		byte[] buf = value.getBytes(IOUtils.ASCII_CHARSET);
		writeString(output, buf, value.length(), maxLnChar);
	}

	/**
	 * Write an 8-bit string.
	 * 
	 * @param output
	 *          output file
	 * @param value
	 *          the value to write
	 */
	public static void writeString8(BufferedFileOutputStream output, String value) throws IOException {
		writeString8(output, value, value.length());
	}

	/**
	 * Write a unicode string.
	 * 
	 * @param output
	 *          the output file
	 * @param value
	 *          the value to be written
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public static void writeString16(BufferedFileOutputStream output, String value) throws IOException {
		writeString16(output, value, value.length());
	}

	/**
	 * Write unicode string.
	 * 
	 * @param output
	 *          output file
	 * @param value
	 *          string value
	 * @param maxLn
	 *          length in bytes to be appended to the value
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public static void writeString16(BufferedFileOutputStream output, String value, int maxLnChar) throws IOException {
		byte[] buf = value.getBytes(IOUtils.UNICODE_CHARSET);
		writeString(output, buf, value.length(), maxLnChar * 2);
	}

	private static void writeString(BufferedFileOutputStream output, byte[] value, int strLn, int maxLnByte)
			throws IOException {
		// how much padding will we need?
		int padLn;
		writeInt32(output, strLn);
		if ((value == null) || (value.length == 0)) {
			padLn = maxLnByte;
		}
		else {
			if (maxLnByte < value.length) {
				output.write(value, 0, maxLnByte);
				return;
			}
			output.write(value);
			padLn = maxLnByte - value.length;
		}
		// padding
		while (padLn > 0) {
			output.write((byte)0);
			padLn--;
		}
	}

	public static void writeBlob(BufferedFileOutputStream output, byte[] value) throws IOException {
		writeBlob(output, value, value.length);
	}

	public static void writeBlob(BufferedFileOutputStream output, byte[] value, int length) throws IOException {
		if (length < value.length) {
			length = value.length;
		}
		writeInt32(output, length);
		output.write(value);
	}

	// / <summary>
	// / Writes a 32 bit unsigned integer to a big endian file.
	// / </summary>
	// / <param name="output">The output file stream.</param>
	// / <param name="value">The value to write to the file.</param>
	public static void writeUInt32(BufferedFileOutputStream output, UInt value) throws IOException {
		output.write(value.getBytes());
	}

	public static void writeInt32(BufferedFileOutputStream output, int value) throws IOException {
		ByteBuffer buf = ByteBuffer.allocate(DataSizes.INT_SIZE);
		buf.putInt(value);
		output.write(buf.order(ByteOrder.BIG_ENDIAN).array());
	}

	public static void writeInt16(BufferedFileOutputStream output, short value) throws IOException {
		ByteBuffer buf = ByteBuffer.allocate(DataSizes.SHORT_SIZE);
		buf.putShort(value);
		output.write(buf.order(ByteOrder.BIG_ENDIAN).array());
	}

	public static void writeUInt16(BufferedFileOutputStream output, UShort value) throws IOException {
		output.write(value.getBytes());
	}

	public static void writeInt8(BufferedFileOutputStream output, byte value) throws IOException {
		ByteBuffer buf = ByteBuffer.allocate(DataSizes.CHAR_SIZE);
		buf.put(value);
		output.write(buf.order(ByteOrder.BIG_ENDIAN).array());
	}

	public static void writeUInt8(BufferedFileOutputStream output, UByte value) throws IOException {
		output.write(value.getBytes());
	}

	public static void writeFloat(BufferedFileOutputStream output, float value) throws IOException {
		ByteBuffer buf = ByteBuffer.allocate(DataSizes.FLOAT_SIZE);
		buf.putFloat(value);
		output.write(buf.order(ByteOrder.BIG_ENDIAN).array());
	}
}
