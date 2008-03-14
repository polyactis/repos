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

package affymetrix.calvin.exception;

/** Base exception class for the Calvin File SDK. */
public class CalvinException extends Exception {

	public static final long serialVersionUID = 0;

	/**
	 * Source in the WIN32 use case is used for registry lookup to determine the resource file. The message code is used
	 * incojuction with source to lookup the specific message string based on code.
	 */
	protected String sourceName;

	/** Message for the exception */
	protected String errorDescription;

	/** Date/Time stamp of exception */
	protected String timeStamp;

	/** File the exception occured in */
	protected String fileName;

	/** Line number that the exception occured */
	protected short lineNumber;

	/** Message code to be used by client for logic or to use for lookup in localized resource file */
	protected long errorCode;

	/**
	 * Constructor Default constructor. Initializes all variables to 0 if numeric, "" if string and the date to the
	 * current date and time.
	 */
	public CalvinException() {
		errorCode = 0;
		lineNumber = 0;
		fileName = "";
		timeStamp = java.util.Calendar.getInstance().toString();
		errorDescription = "";
		sourceName = "";
	}

	/**
	 * Constructor Constructs a CalvinException object.
	 * 
	 * @param _Source
	 *          Typically will be the application name. Can also be the file sdk object type E.g. "Array File Reader".
	 * @param _Description
	 *          A brief message which describes what happened.
	 * @param _TimeStamp
	 *          Time and date the exception occured.
	 * @param _FileName
	 *          File name of the source file the exception occured. The __FILE__ can be used to determine this
	 *          information.
	 * @param _LineNumber
	 *          Line number in the source file that generated the exception. The __LINE__ can be used to determine this
	 *          information.
	 * @param _ErrorCode
	 *          An numeric value the is unique to this error/exception type
	 */
	public CalvinException(String _Source, String _Description, String _TimeStamp, String _FileName, short _LineNumber,
			long _ErrorCode) {
		sourceName = _Source;
		errorDescription = _Description;
		lineNumber = _LineNumber;
		fileName = _FileName;
		timeStamp = _TimeStamp;
		errorCode = _ErrorCode;
	}

	/**
	 * Constructor Constructs a CalvinException object Note: The time date stamp will be automatically set to the current
	 * time.
	 * 
	 * @param _Source
	 *          Typically will be the application name. Can also be the file sdk object type E.g. "Array File Reader".
	 * @param _Description
	 *          A brief message which describes what happened.
	 * @param _ErrorCode
	 *          An numeric value the is unique to this error/exception type
	 */
	public CalvinException(String _Source, String _Description, long _ErrorCode) {
		sourceName = _Source;
		errorDescription = _Description;
		lineNumber = 0;
		fileName = "";
		timeStamp = java.util.Calendar.getInstance().toString();
		errorCode = _ErrorCode;
	}

	/**
	 * Constructor Constructs a CalvinException object Note: The time date stamp will be automatically set to the current
	 * time.
	 * 
	 * @param _Description
	 *          A brief message which describes what happened.
	 * @param _ErrorCode
	 *          An numeric value the is unique to this error/exception type
	 */
	public CalvinException(String _Description, long _ErrorCode) {
		sourceName = "";
		errorDescription = _Description;
		lineNumber = 0;
		fileName = "";
		timeStamp = java.util.Calendar.getInstance().toString();
		errorCode = _ErrorCode;
	}

	/**
	 * Constructor Constructs a CalvinException object Note: The time date stamp will be automatically set to the current
	 * time.
	 * 
	 * @param _ErrorCode
	 *          An numeric value the is unique to this error/exception type
	 */
	public CalvinException(long _ErrorCode) {
		sourceName = "";
		errorDescription = "";
		lineNumber = 0;
		fileName = "";
		timeStamp = java.util.Calendar.getInstance().toString();
		errorCode = _ErrorCode;
	}

	/**
	 * The source name associated with the exception.
	 * 
	 * @return The source name associated with the exception.
	 */
	public String getSource() {
		return sourceName;
	}

	/**
	 * The source name associated with the exception.
	 * 
	 * @param value
	 *          Source name associated with the exception.
	 */
	public void setSource(String value) {
		sourceName = value;
	}

	/**
	 * The error description associated with the exception.
	 * 
	 * @return The error description associated with the exception.
	 */
	public String getDescription() {
		return errorDescription;
	}

	/**
	 * The error description associated with the exception.
	 * 
	 * @param value
	 *          Error description associated with the exception.
	 */
	public void setDescription(String value) {
		errorDescription = value;
	}

	/**
	 * The error time stamp associated with the exception.
	 * 
	 * @return The error time stamp associated with the exception.
	 */
	public String getTimeStamp() {
		return timeStamp;
	}

	/**
	 * The error time stamp associated with the exception.
	 * 
	 * @param value
	 *          Error time stamp associated with the exception.
	 */
	public void setTimeStamp(String value) {
		timeStamp = value;
	}

	/**
	 * The error source file name associated with the exception.
	 * 
	 * @return The source file name associated with the exception.
	 */
	public String getSourceFile() {
		return fileName;
	}

	/**
	 * The error source file name associated with the exception.
	 * 
	 * @param value
	 *          Source file name associated with the exception.
	 */
	public void setSourceFile(String value) {
		fileName = value;
	}

	/**
	 * The error source line number associated with the exception.
	 * 
	 * @return The source line number associated with the exception.
	 */
	public short getLineNumber() {
		return lineNumber;
	}

	/**
	 * The error source line number associated with the exception.
	 * 
	 * @param value
	 *          Source line number associated with the exception.
	 */
	public void setLineNumber(short value) {
		lineNumber = value;
	}

	/**
	 * The error code associated with the exception.
	 * 
	 * @return The error code associated with the exception.
	 */
	public long getErrorCode() {
		return errorCode;
	}

	/**
	 * The error code associated with the exception.
	 * 
	 * @param value
	 *          Error code associated with the exception.
	 */
	public void setErrorCode(long value) {
		errorCode = value;
	}

	/**
	 * Returns a string describing the exception
	 * 
	 * @return Returns a string describing the exception.
	 */
	@Override
	public String toString() {
		return "";
	}

}
