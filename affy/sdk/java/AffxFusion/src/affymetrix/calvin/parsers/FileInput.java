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
import java.io.UnsupportedEncodingException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.portability.DataSizes;
import affymetrix.portability.UByte;
import affymetrix.portability.UInt;
import affymetrix.portability.UShort;

/** Provides functions for reading atom data items from a file or stream. */
public class FileInput {

	/**
	 * Reads an 8 bit integer from a big endian file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return The integer read from the file.
	 */
	public static byte readInt8(FileInputStream fis) throws IOException {
		return (byte)fis.read();
	}

	/**
	 * Reads a 16 bit integer from a big endian file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return The integer read from the file.
	 */
	public static short readInt16(FileInputStream fis) throws IOException {
		byte[] b = new byte[DataSizes.SHORT_SIZE];
		fis.read(b);
		return ByteBuffer.wrap(b).order(ByteOrder.BIG_ENDIAN).getShort();

	}

	/**
	 * Reads a 32 bit integer from a big endian file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return The integer read from the file.
	 */
	public static int readInt32(FileInputStream fis) throws IOException {
		byte[] b = new byte[DataSizes.INT_SIZE];
		fis.read(b);
		return ByteBuffer.wrap(b).order(ByteOrder.BIG_ENDIAN).getInt();
	}

	/**
	 * Reads an 8 bit unsigned integer from a big endian file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return The integer read from the file.
	 */
	public static UByte readUInt8(FileInputStream fis) throws IOException, UnsignedOutOfLimitsException {
		byte[] b = new byte[DataSizes.SHORT_SIZE];
		fis.read(b, 1, DataSizes.CHAR_SIZE);
		return new UByte(ByteBuffer.wrap(b).order(ByteOrder.BIG_ENDIAN).getShort());
	}

	/**
	 * Reads a 16 bit unsigned integer from a big endian file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return The integer read from the file.
	 */
	public static UShort readUInt16(FileInputStream fis) throws IOException, UnsignedOutOfLimitsException {
		byte[] b = new byte[DataSizes.INT_SIZE];
		fis.read(b, 2, DataSizes.SHORT_SIZE);
		return new UShort(ByteBuffer.wrap(b).order(ByteOrder.BIG_ENDIAN).getInt());

	}

	/**
	 * Reads a 32 bit unsigned integer from a big endian file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return The integer read from the file.
	 */
	public static UInt readUInt32(FileInputStream fis) throws IOException, UnsignedOutOfLimitsException {
		byte[] b = new byte[8];
		fis.read(b, 4, DataSizes.INT_SIZE);
		return new UInt(ByteBuffer.wrap(b).order(ByteOrder.BIG_ENDIAN).getLong());

	}

	/**
	 * Reads a 32 bit floating point number from a big endian file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return The floating point number read from the file.
	 */
	public static float readFloat(FileInputStream fis) throws IOException {
		byte[] b = new byte[DataSizes.FLOAT_SIZE];
		fis.read(b);
		return ByteBuffer.wrap(b).order(ByteOrder.BIG_ENDIAN).getFloat();
	}

	/**
	 * Reads an 8 bit integer from a big endian file stream (memory map pointer).
	 * 
	 * @param b
	 *          The input file stream.
	 * @return The integer read from the file stream.
	 */
	public static byte readInt8(ByteBuffer b) {
		return b.get();
	}

	/**
	 * Reads a 16 bit integer from a big endian file stream (memory map pointer).
	 * 
	 * @param b
	 *          The input file stream.
	 * @return The integer read from the file stream.
	 */
	public static short readInt16(ByteBuffer b) {
		return b.getShort();
	}

	/**
	 * Reads a 32 bit integer from a big endian file stream (memory map pointer).
	 * 
	 * @param b
	 *          The input file stream.
	 * @return The integer read from the file stream.
	 */
	public static int readInt32(ByteBuffer b) {
		return b.getInt();
	}

	/**
	 * Reads an 8 bit unsigned integer from a big endian file stream (memory map pointer).
	 * 
	 * @param b
	 *          The input file stream.
	 * @return The integer read from the file stream.
	 */
	public static UByte readUInt8(ByteBuffer b) throws UnsignedOutOfLimitsException {
		return new UByte((short)(0xFF & (int)b.get()));
	}

	/**
	 * Reads a 16 bit unsigned integer from a big endian file stream (memory map pointer).
	 * 
	 * @param b
	 *          The input file stream.
	 * @return The integer read from the file stream.
	 */
	public static UShort readUInt16(ByteBuffer b) throws UnsignedOutOfLimitsException {
		return new UShort(0xFFFF & (int)b.getShort());
	}

	/**
	 * Reads a 32 bit unsigned integer from a big endian file stream (memory map pointer).
	 * 
	 * @param b
	 *          The input file stream.
	 * @return The integer read from the file stream.
	 */
	public static UInt readUInt32(ByteBuffer b) throws UnsignedOutOfLimitsException {
		return new UInt(0xFFFFFFFFL & (long)b.getInt());
	}

	/**
	 * Reads a 32 bit floating point number from a big endian file stream (memory map pointer).
	 * 
	 * @param b
	 *          The input file stream.
	 * @return The floating point number read from the file stream.
	 */
	public static float readFloat(ByteBuffer b) {
		return b.getFloat();
	}

	/**
	 * Reads a 16 bit unicode string of fixed size from a big endian file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @param len
	 *          The length of the string.
	 * @return The string read from the file stream.
	 */
	public static String readString16(FileInputStream fis, int len) throws IOException {
		byte[] b = new byte[DataSizes.SHORT_SIZE * len];
		fis.read(b);
		return new String(b, "UTF-16BE").trim();
	}

	/**
	 * Reads a 16 bit unicode string from a big endian file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return The string read from the file stream.
	 */
	public static String readString16(FileInputStream fis) throws IOException {
		int len = readInt32(fis);
		return readString16(fis, len);
	}

	/**
	 * Reads an 8 bit string of fixed size from a big endian file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @param len
	 *          The length of the string.
	 * @return The string read from the file stream.
	 */
	public static String readString8(FileInputStream fis, int len) throws IOException {
		byte[] b = new byte[len];
		fis.read(b);
		return new String(b, "US-ASCII").trim();
	}

	/**
	 * Reads an 8 bit string from a big endian file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return The string read from the file stream.
	 */
	public static String readString8(FileInputStream fis) throws IOException {
		int len = readInt32(fis);
		return readString8(fis, len);
	}

	/**
	 * Reads a 16 bit unicode string of fixed size from a big endian file stream (memory map pointer).
	 * 
	 * @param b
	 *          The input file stream.
	 * @param len
	 *          The length of the string.
	 * @return The string read from the file stream.
	 */
	public static String readString16(ByteBuffer b, int len) throws UnsupportedEncodingException {
		byte[] buf = new byte[len];
		b.get(buf);
		return new String(buf, "UTF-16BE").trim();
	}

	/**
	 * Reads a 16 bit unicode string from a big endian file stream (memory map pointer).
	 * 
	 * @param b
	 *          The input file stream.
	 * @return The string read from the file stream.
	 */
	public static String readString16(ByteBuffer b) throws IOException {
		int len = readInt32(b);
		return readString16(b, len);
	}

	/**
	 * Reads an 8 bit string of fixed size from a big endian file stream (memory map pointer).
	 * 
	 * @param b
	 *          The input file stream.
	 * @param len
	 *          The length of the string.
	 * @return The string read from the file stream.
	 */
	public static String readString8(ByteBuffer b, int len) throws UnsupportedEncodingException {
		byte[] buf = new byte[len];
		b.get(buf);
		return new String(buf, "US-ASCII").trim();
	}

	/**
	 * Reads an 8 bit string from a big endian file stream (memory map pointer).
	 * 
	 * @param b
	 *          The input file stream.
	 * @return The string read from the file stream.
	 */
	public static String readString8(ByteBuffer b) throws IOException {
		int len = readInt32(b);
		return readString8(b, len);
	}

	/**
	 * Reads a blob from a file.
	 */
	public static byte[] readBlob(FileInputStream fis) throws IOException {
		int size = FileInput.readInt32(fis);
		byte[] value = new byte[size];
		fis.read(value);
		return value;
	}
}
