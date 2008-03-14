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

package affymetrix.calvin.parameter;

import java.io.UnsupportedEncodingException;

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.utils.IOUtils;
import affymetrix.portability.DataSizes;
import affymetrix.portability.UByte;
import affymetrix.portability.UInt;
import affymetrix.portability.UShort;

/**
 * A class to hold a parameter name/value/type. This class will convert several built-in types between the MIME string
 * and their native types.
 */
public class ParameterNameValue {

	// /** an 8 bit integer. */
	// public static final int Int8Type = 0;
	//
	// /** an 8 bit unsigned integer. */
	// public static final int UInt8Type = 1;
	//
	// /** a 16 bit integer. */
	// public static final int Int16Type = 2;
	//
	// /** a 16 bit unsigned integer. */
	// public static final int UInt16Type = 3;
	//
	// /** a 32 bit integer. */
	// public static final int Int32Type = 4;
	//
	// /** a 32 bit unsigned integer. */
	// public static final int UInt32Type = 5;
	//
	// /** a 32 bit floating point. */
	// public static final int FloatType = 6;
	//
	// /** a 16 bit character. */
	// public static final int TextType = 7;
	//
	// /** an 8 bit character. */
	// public static final int AsciiType = 8;
	//
	// /** an 8 bit integer. */
	// public static final int UnknownType = 9;

	/* built-in parameter types. */
	public enum ParameterType {
		Int8Type, /* an 8 bit integer. */
		UInt8Type, /* an 8 bit unsigned integer. */
		Int16Type, /* a 16 bit integer. */
		UInt16Type, /* a 16 bit unsigned integer. */
		Int32Type, /* a 32 bit integer. */
		UInt32Type, /* a 32 bit unsigned integer. */
		FloatType, /* a 32 bit floating point. */
		TextType, /* a 16 bit character. */
		AsciiType, /* an 8 bit character. */
		UnknownType
		/* an 8 bit integer. */
	};

	/** The name of the parameter */
	private String name;

	/** The MIME type of the parameter */
	private String type;

	/** The MIME value of the parameter */
	private MIMEValue value = new MIMEValue();

	/** Creates a new instance of ParameterNameValue */
	public ParameterNameValue() {
		name = IOUtils.EMPTY;
		type = IOUtils.EMPTY;
	}

	/**
	 * Constructor. Useful when reading in from a file
	 * 
	 * @param n
	 *          Parameter name
	 * @param mimeValue
	 *          MIME encoded value in a buffer
	 * @param mimeValueSize
	 *          The size in bytes of the MIME encoded value.
	 * @param mimeType
	 *          The value type.
	 */
	public ParameterNameValue(String n, byte[] mimeValue, String mimeType) {
		name = n;
		value.setValue(mimeValue);
		type = mimeType;
	}

	/**
	 * Constructor. Useful when reading in from a file
	 * 
	 * @param n
	 *          Parameter name
	 * @param mimeValue
	 *          MIME encoded value
	 * @param mimeType
	 */
	public ParameterNameValue(String n, MIMEValue mimeValue, String mimeType) {
		this(n, mimeValue.getBytes(), mimeType);
	}

	/**
	 * A copy constructor.
	 * 
	 * @param param
	 *          The parameter to copy.
	 */
	public ParameterNameValue(ParameterNameValue param) {
		this(param.getName(), param.getMIMEValue(), param.getMIMEType());
	}

	/**
	 * An equality operator. Compares the parameter name.
	 * 
	 * @param param
	 *          The parameter to compare.
	 * @return True if the parameter names are the same.
	 */
	boolean equals(ParameterNameValue param) {
		return equals(param.getName());
	}

	/**
	 * An equality operator. Compares the parameter name.
	 * 
	 * @param n
	 *          A parameter name to compare.
	 * @return True if the parameter names are the same.
	 */
	public boolean equals(String n) {
		return name.equals(n);
	}

	/**
	 * get the name of the parameter
	 * 
	 * @return The name of the parameter
	 */
	public String getName() {
		return name;
	}

	/**
	 * set the name of the parameter.
	 * 
	 * @param v
	 *          The name of the parameter.
	 */
	public void setName(String v) {
		name = v;
	}

	/**
	 * get the parameter type
	 * 
	 * @return The parameter type of the object.
	 */
	public ParameterType getParameterType() {
		int ln = ParameterType.UnknownType.ordinal();
		for (int i = 0; i < ln; ++i) {
			if (type.equals(MIMEValue.TYPE_TABLE[i])) {
				return ParameterType.values()[i];
			}
		}
		return ParameterType.UnknownType;
	}

	/**
	 * gets value as a int8_t
	 * 
	 * @return Value of the parameter as an int8_t
	 */
	public byte getValueInt8() {
		return (byte)valueToInt(value);
	}

	/**
	 * sets the value as an int8_t
	 * 
	 * @param v
	 */
	public void setValueInt8(byte v) {
		type = MIMEValue.Int8MIMEType;
		intToValue(v, value);
	}

	/**
	 * gets value as a byte
	 * 
	 * @return Value of the parameter as an byte
	 */
	public UByte getValueUInt8() throws UnsignedOutOfLimitsException {
		return new UByte((short)(0xFF & (int)value.getBytes()[3]));
	}

	/**
	 * sets the value as a byte
	 * 
	 * @param v
	 */
	public void setValueUInt8(UByte v) {
		type = MIMEValue.UInt8MIMEType;
		intToValue(v.toShort(), value);
	}

	/**
	 * gets value as a int16_t
	 * 
	 * @return Value of the parameter as an int16_t
	 */
	public short getValueInt16() {
		return (short)valueToInt(value);
	}

	/**
	 * sets the value as an int16_t
	 * 
	 * @param v
	 */
	public void setValueInt16(short v) {
		type = MIMEValue.Int16MIMEType;
		intToValue(v, value);
	}

	/**
	 * gets value as a u_int16_t
	 * 
	 * @return Value of the parameter as an u_int16_t
	 */
	public UShort getValueUInt16() throws UnsignedOutOfLimitsException {
		byte[] buf = new byte[2];
		System.arraycopy(value.getBytes(), 2, buf, 0, 2);
		return new UShort(buf);
	}

	/**
	 * sets the value as a u_int16_t
	 * 
	 * @param v
	 */
	public void setValueUInt16(UShort v) {
		type = MIMEValue.UInt16MIMEType;
		intToValue(v.toInt(), value);
	}

	/**
	 * gets value as a int
	 * 
	 * @return Value of the parameter as an int
	 */
	public int getValueInt32() {
		return valueToInt(value);
	}

	/**
	 * sets the value as an int
	 * 
	 * @param v
	 */
	public void setValueInt32(int v) {
		type = MIMEValue.Int32MIMEType;
		intToValue(v, value);
	}

	/**
	 * gets value as a int
	 * 
	 * @return Value of the parameter as an int
	 */
	public UInt getValueUInt32() throws UnsignedOutOfLimitsException {
		return new UInt(value.getBytes());
	}

	/**
	 * sets the value as a int
	 * 
	 * @param v
	 */
	public void setValueUInt32(UInt v) {
		type = MIMEValue.UInt32MIMEType;
		value.setValue(v.getBytes());
	}

	/**
	 * gets value as a float
	 * 
	 * @return Value of the parameter as an float
	 */
	public float getValueFloat() {
		int ival = valueToInt(value);
		return Float.intBitsToFloat(ival);
	}

	/**
	 * sets the value as a float
	 * 
	 * @param v
	 */
	public void setValueFloat(float v) {
		type = MIMEValue.FloatMIMEType;
		intToValue(Float.floatToIntBits(v), value);
	}

	/**
	 * Gets the value as a wstring
	 * 
	 * @param val
	 *          The MIME value object.
	 * @return Value of the parameter as an wstring
	 */
	protected static String getTextFromMIMEValue(MIMEValue val) {
		return getStringFromMIMEValue(val, "UTF-16BE");
	}

	/**
	 * gets value as a wstring
	 * 
	 * @return Value of the parameter as an wstring
	 */
	public String getValueText() {
		return getTextFromMIMEValue(value);
	}

	/**
	 * sets the value as a text type.
	 * 
	 * @param v
	 *          String representation of the value.
	 */
	public void setValueText(String v) {
		setValueText(v, 0);
	}

	public void setValueText(String v, int reserve) {
		type = MIMEValue.TextMIMEType;
		value.setValue(getBytes(v, "UTF-16BE", reserve));
	}

	protected static byte[] getBytes(String value, String encoding, int reserve) {
		byte[] result = null;
		try {
			result = value.getBytes(encoding);
			if (reserve > result.length) {
				byte[] buf = new byte[reserve];
				System.arraycopy(result, 0, buf, 0, result.length);
				result = buf;
			}
		}
		catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		}
		return result;
	}

	/**
	 * gets value as a string
	 * 
	 * @return Value of the parameter as an string
	 */
	protected static String getAsciiFromMIMEValue(MIMEValue val) {
		return getStringFromMIMEValue(val, "US-ASCII");
	}

	private static String getStringFromMIMEValue(MIMEValue val, String encoding) {
		String result = null;
		try {
			result = new String(val.getBytes(), encoding);
		}
		catch (UnsupportedEncodingException e) {

		}
		return result.trim();
	}

	/**
	 * gets value as a string
	 * 
	 * @return Value of the parameter as an string
	 */
	public String getValueAscii() {
		return getAsciiFromMIMEValue(value);
	}

	/**
	 * sets the value as a text type.
	 * 
	 * @param v
	 *          String representation of the value.
	 */
	public void setValueAscii(String v) {
		setValueAscii(v, 0);
	}

	public void setValueAscii(String v, int reserve) {
		type = MIMEValue.AsciiMIMEType;
		value.setValue(getBytes(v, "US-ASCII", reserve));
	}

	// Raw MIME methods
	/**
	 * Returns the mime type without interpretation.
	 * 
	 * @return String of the MIME type
	 */
	public String getMIMEType() {
		return type;
	}

	/**
	 * sets the mime type without attempting to interpret it.
	 * 
	 * @param v
	 *          The mime type.
	 */
	public void setMIMEType(String v) {
		type = v;
	}

	/**
	 * Returns the mime value without interpretation.
	 * 
	 * @return MIME encoded string.
	 */
	public MIMEValue getMIMEValue() {
		return value;
	}

	/**
	 * sets the mime value without attempting to interpret it.
	 * 
	 * @param v
	 *          MIME encoded string.
	 */
	public void setMIMEValue(MIMEValue v) {
		value.setValue(v.getBytes());
	}

	/**
	 * Converts known types to a string
	 * 
	 * @return A string representation of the value
	 */
	@Override
	public String toString() {
		try {
			switch (getParameterType()) {
			case Int8Type:
				return String.valueOf(getValueInt8());

			case Int16Type:
				return String.valueOf(getValueInt16());

			case Int32Type:
				return String.valueOf(getValueInt32());

			case UInt8Type:
				return String.valueOf(getValueUInt8());

			case UInt16Type:
				return String.valueOf(getValueUInt16());

			case UInt32Type:
				return String.valueOf(getValueUInt32().toLong());

			case FloatType:
				return String.valueOf(getValueFloat());

			case TextType:
				return getValueText();

			case AsciiType:
				return getValueAscii();
			}
		}
		catch (UnsignedOutOfLimitsException e) {
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * Converts a value to an integer.
	 * 
	 * @param val
	 *          The mime representation.
	 * @return The integer representation.
	 */
	protected static int valueToInt(MIMEValue val) {
		int ival = 0;
		byte[] buf = val.getBytes();
		if (buf[0] < 0) {
			ival = ((buf[0] + 256) << 24);
		}
		else {
			ival = (buf[0] << 24);
		}
		if (buf[1] < 0) {
			ival += ((buf[1] + 256) << 16);
		}
		else {
			ival += (buf[1] << 16);
		}
		if (buf[2] < 0) {
			ival += ((buf[2] + 256) << 8);
		}
		else {
			ival += (buf[2] << 8);
		}
		if (buf[3] < 0) {
			ival += (buf[3] + 256);
		}
		else {
			ival += buf[3];
		}
		return ival;
	}

	/**
	 * Converts an integer to a value.
	 * 
	 * @param v
	 *          The integer representation.
	 * @param val
	 *          The value.
	 */
	protected static void intToValue(int v, MIMEValue val) {
		byte[] buf = new byte[DataSizes.INT_SIZE];
		buf[0] = (byte)(v >> 24);
		buf[1] = (byte)(v >> 16);
		buf[2] = (byte)(v >> 8);
		buf[3] = (byte)v;
		val.setValue(buf);
	}

}
