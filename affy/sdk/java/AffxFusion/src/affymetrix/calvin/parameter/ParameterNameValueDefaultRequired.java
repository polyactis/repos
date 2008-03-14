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

import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.utils.IOUtils;

/** A class to hold a name/value/type/default/required attributes. */
public class ParameterNameValueDefaultRequired extends ParameterNameValue {

	public static enum ParameterValueType {
		NoParameterType, IntegerParameterType, FloatParameterType, TextParameterType, DateParameterType, TimeParameterType, DateTimeParameterType, ControlSingleParameterType, ControlMultiParameterType
	};

	public static final String[] TYPE_TABLE = { IOUtils.EMPTY, "Int", "Float", "String", "Date", "Time", "DateTime",
			"SingleControl", "MultiControl" };

	/**
	 * Converts the type to a string.
	 * 
	 * @param t
	 *          The type.
	 * @return The string representation.
	 */
	public static String ParameterValueTypeToString(ParameterValueType t) {
		return TYPE_TABLE[t.ordinal()];
	}

	/**
	 * Converts the string to a type.
	 * 
	 * @param str
	 *          The string representation.
	 * @return The type.
	 */
	public static int ParameterValueTypeFromString(String str) {
		for (int i = 0; i < ParameterValueType.ControlMultiParameterType.ordinal(); i++) {
			if (str.equals(TYPE_TABLE[i])) {
				return i;
			}
		}
		return ParameterValueType.NoParameterType.ordinal();
	}

	/** The MIME value of the default parameter */
	private MIMEValue defaultValue = new MIMEValue();

	/** A flag to indicate if a default exist. */
	private boolean hasDefault = false;

	/** A flag to indicate if the parameter is required. */
	private boolean required = false;

	/** A list of parameter values for controlled vocabulary. */
	private List<String> controlled = new ArrayList<String>();

	/** A list of multi-selected controlled values. */
	private List<String> controlMultiValues = new ArrayList<String>();

	/** The type of value stored. */
	private ParameterValueType valueType = ParameterValueType.NoParameterType;

	/** Sets the parameter type. */
	private void setParameterType() {
		ParameterType paramType = getParameterType();
		switch (paramType) {
		case Int8Type:
		case UInt8Type:
		case Int16Type:
		case UInt16Type:
		case Int32Type:
		case UInt32Type:
			valueType = ParameterValueType.IntegerParameterType;
			break;

		case FloatType:
			valueType = ParameterValueType.FloatParameterType;
			break;

		case TextType:
		case AsciiType:
			valueType = ParameterValueType.TextParameterType;
			break;

		default:
			break;
		}
	}

	/**
	 * Gets a vector of multi-selected controlled values.
	 * 
	 * @return The multi-selected controlled values.
	 */
	public List<String> getControlMultiValues() {
		return controlMultiValues;
	}

	/** The type of value stored. */
	public ParameterValueType getValueType() {
		return valueType;
	}

	/** The type of value stored. */
	public void setValueType(ParameterValueType v) {
		valueType = v;
	}

	/**
	 * Gets a vector of parameter values for controlled vocabulary.
	 * 
	 * @return The vector of parameter values for controlled vocabulary.
	 */
	public List<String> getControlledVocabulary() {
		return controlled;
	}

	/**
	 * Returns the mime value without interpretation.
	 * 
	 * @return MIME encoded string.
	 */
	public MIMEValue getDefaultMIMEValue() {
		return defaultValue;
	}

	/**
	 * Gets the required flag.
	 * 
	 * @return The required flag.
	 */
	public boolean getRequiredFlag() {
		return required;
	}

	/**
	 * Sets the required flag.
	 * 
	 * @param r
	 *          The required flag.
	 */
	public void setRequiredFlag(boolean r) {
		required = r;
	}

	/**
	 * Gets a flag to indicate if a default exist.
	 * 
	 * @return The flag indicating if a default value exists.
	 */
	public boolean getHasDefault() {
		return hasDefault;
	}

	/**
	 * Sets the flag to indicate if a default exist.
	 * 
	 * @param d
	 *          The flag indicating if a default value exists.
	 */
	public void setHasDefault(boolean d) {
		hasDefault = d;
	}

	public ParameterNameValueDefaultRequired() {
		super();
	}

	/**
	 * MIMEructor. Useful when reading in from a file
	 * 
	 * @param n
	 *          Parameter name
	 * @param mimeValue
	 *          MIME encoded value
	 * @param mimeType
	 *          The MIME type.
	 * @param defaultMimeValue
	 *          The default MIME encoded value in a buffer.
	 * @param req
	 *          Flag to indicate if the parameter is required.
	 */
	public ParameterNameValueDefaultRequired(String n, MIMEValue mimeValue, String mimeType, MIMEValue defaultMimeValue,
			boolean req) {
		this(n, mimeValue.getBytes(), mimeType, defaultMimeValue.getBytes(), req);
	}

	/**
	 * MIMEructor. Useful when reading in from a file
	 * 
	 * @param n
	 *          Parameter name
	 * @param mimeValue
	 *          MIME encoded value in a buffer
	 * @param mimeValueSize
	 *          The size in bytes of the MIME encoded value.
	 * @param mimeType
	 *          The MIME type.
	 * @param defaultMimeValue
	 *          The default MIME encoded value in a buffer.
	 * @param defaultMimeValueSize
	 *          The size in bytes of the default MIME encoded value.
	 * @param req
	 *          Flag to indicate if the parameter is required.
	 */
	public ParameterNameValueDefaultRequired(String n, byte[] mimeValue, String mimeType, byte[] defaultMimeValue,
			boolean req) {
		super(n, mimeValue, mimeType);
		hasDefault = true;
		defaultValue.setValue(defaultMimeValue);
		required = req;
		setParameterType();
	}

	/**
	 * Converts default value to a string
	 * 
	 * @return A string representation of the default value
	 */
	public String defaultToString() {
		switch (getParameterType()) {
		case Int8Type:
			return String.valueOf(getDefaultValueInt8());

		case Int16Type:
			return String.valueOf(getDefaultValueInt16());

		case Int32Type:
			return String.valueOf(getDefaultValueInt32());

		case UInt8Type:
			return String.valueOf(getDefaultValueUInt8());

		case UInt16Type:
			return String.valueOf(getDefaultValueUInt16());

		case UInt32Type:
			return String.valueOf(getDefaultValueUInt32());

		case FloatType:
			return String.valueOf(getDefaultValueFloat());

		case TextType:
			return getDefaultValueText();

		case AsciiType:
			return getDefaultValueAscii();
		}
		return null;
	}

	/**
	 * Gets default value as a byte
	 * 
	 * @return Value of the parameter as an byte
	 * @exception affymetrix_calvin_exceptions::ParameterMismatchException
	 *              Parameter is not an byte
	 */
	public byte getDefaultValueInt8() {
		return (byte)valueToInt(defaultValue);
	}

	/**
	 * Sets the default value as a byte
	 * 
	 * @param v
	 *          The value
	 */
	public void setDefaultValueInt8(byte v) {
		hasDefault = true;
		super.setMIMEType(MIMEValue.Int8MIMEType);
		intToValue(v, defaultValue);
	}

	/**
	 * Gets the default value as a u_byte
	 * 
	 * @return Value of the parameter as an u_byte
	 * @exception affymetrix_calvin_exceptions::ParameterMismatchException
	 *              Parameter is not a u_byte
	 */
	public byte getDefaultValueUInt8() {
		return (byte)valueToInt(defaultValue);
	}

	/**
	 * Sets the default value as a u_byte
	 * 
	 * @param v
	 *          The value
	 */
	public void setDefaultValueUInt8(byte v) {
		hasDefault = true;
		super.setMIMEType(MIMEValue.UInt8MIMEType);
		intToValue(v, defaultValue);
	}

	/**
	 * Gets the default value as a short
	 * 
	 * @return Value of the parameter as an short
	 * @exception affymetrix_calvin_exceptions::ParameterMismatchException
	 *              Parameter is not an short
	 */
	public short getDefaultValueInt16() {
		return (short)valueToInt(defaultValue);
	}

	/**
	 * Sets the default value as an short
	 * 
	 * @param v
	 *          The value
	 */
	public void setDefaultValueInt16(short v) {
		hasDefault = true;
		super.setMIMEType(MIMEValue.Int16MIMEType);
		intToValue(v, defaultValue);
	}

	/**
	 * Gets the default value as a u_short
	 * 
	 * @return Value of the parameter as an u_short
	 * @exception affymetrix_calvin_exceptions::ParameterMismatchException
	 *              Parameter is not a u_short
	 */
	public short getDefaultValueUInt16() {
		return (short)valueToInt(defaultValue);
	}

	/**
	 * Sets the default value as a u_short
	 * 
	 * @param v
	 *          The value
	 */
	public void setDefaultValueUInt16(short v) {
		hasDefault = true;
		super.setMIMEType(MIMEValue.UInt16MIMEType);
		intToValue(v, defaultValue);
	}

	/**
	 * Gets the default value as a int
	 * 
	 * @return Value of the parameter as an int
	 * @exception affymetrix_calvin_exceptions::ParameterMismatchException
	 *              Parameter is not an int
	 */
	public int getDefaultValueInt32() {
		return valueToInt(defaultValue);
	}

	/**
	 * Sets the default value as an int
	 * 
	 * @param v
	 *          The value
	 */
	public void setDefaultValueInt32(int v) {
		hasDefault = true;
		super.setMIMEType(MIMEValue.Int32MIMEType);
		intToValue(v, defaultValue);
	}

	/**
	 * Gets the default value as a int
	 * 
	 * @return Value of the parameter as an int
	 * @exception affymetrix_calvin_exceptions::ParameterMismatchException
	 *              Parameter is not a int
	 */
	public int getDefaultValueUInt32() {
		return valueToInt(defaultValue);
	}

	/**
	 * Sets the default value as a int
	 * 
	 * @param v
	 *          The value
	 */
	public void setDefaultValueUInt32(int v) {
		hasDefault = true;
		super.setMIMEType(MIMEValue.UInt32MIMEType);
		intToValue(v, defaultValue);
	}

	/**
	 * Gets the default value as a float
	 * 
	 * @return Value of the parameter as an float
	 * @exception affymetrix_calvin_exceptions::ParameterMismatchException
	 *              Parameter is not a float
	 */
	public float getDefaultValueFloat() {
		int ival = valueToInt(defaultValue);
		return Float.intBitsToFloat(ival);
	}

	/**
	 * Sets the default value as a float
	 * 
	 * @param v
	 *          The value
	 */
	public void setDefaultValueFloat(float v) {

		hasDefault = true;
		super.setMIMEType(MIMEValue.FloatMIMEType);
		intToValue(Float.floatToRawIntBits(v), defaultValue);
	}

	/**
	 * Gets the default value as a wstring
	 * 
	 * @return Value of the parameter as an wstring
	 * @exception affymetrix_calvin_exceptions::ParameterMismatchException
	 *              Parameter is not a text type
	 */
	public String getDefaultValueText() {
		return getTextFromMIMEValue(defaultValue);
	}

	/**
	 * Sets the default value as a text type.
	 * 
	 * @param v
	 *          String representation of the default value.
	 */
	public void setDefaultValueText(String v) {
		super.setMIMEType(MIMEValue.TextMIMEType);
		defaultValue.setValue(getBytes(v, "UTF-16BE", 0));
	}

	/**
	 * Gets the default value as a string
	 * 
	 * @return Value of the parameter as an string
	 * @exception affymetrix_calvin_exceptions::ParameterMismatchException
	 *              Parameter is not a text type
	 */
	public String getDefaultValueAscii() {
		return getAsciiFromMIMEValue(defaultValue);
	}

	/**
	 * Sets the default value as a text type.
	 * 
	 * @param v
	 *          String representation of the default value.
	 */
	public void setDefaultValueAscii(String v) {
		super.setMIMEType(MIMEValue.AsciiMIMEType);
		defaultValue.setValue(getBytes(v, "US-ASCII", 0));
	}
}
