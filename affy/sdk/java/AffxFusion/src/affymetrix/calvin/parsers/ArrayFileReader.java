package affymetrix.calvin.parsers;

import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.XMLReaderFactory;

import affymetrix.calvin.array.ArrayData;
import affymetrix.calvin.exception.SAXArrayStopParsingException;
import affymetrix.calvin.utils.AffymetrixGuidType;

public class ArrayFileReader {
	public ArrayFileReader() {
	}

	public boolean read(String fileName, ArrayData arrayData) {
		return read(fileName, arrayData, false);
	}

	/**
	 * Read the entire file using the XML SAX parser.
	 */
	boolean read(String fileName, ArrayData arrayData, boolean headerOnly) {
		arrayData.clear();

		// // Initialize the XML4C2 system
		// try
		// {
		// XMLPlatformUtils::Initialize();
		// }
		// catch ( XMLException)
		// {
		// return false;
		// }

		boolean status = false;
		// SAXParser parser = new SAXParser();
		// parser->setValidationScheme(SAXParser::Val_Never);
		// parser->setLoadExternalDTD(false);
		// parser->setDoNamespaces(false);
		// parser->setDoSchema(false);
		// parser->setValidationSchemaFullChecking(false);
		XMLReader reader = null;
		try {
			reader = XMLReaderFactory.createXMLReader();
		}
		catch (SAXException se) {
			se.printStackTrace();
		}

		SAXArrayHandler handler = new SAXArrayHandler(arrayData, headerOnly);
		reader.setContentHandler(handler);

		// SAXArrayHandlers handler(arrayData, headerOnly);
		// parser->setDocumentHandler(handler);
		// parser->setErrorHandler(handler);

		try {
			reader.parse(fileName);
			// int errorCount = parser->getErrorCount();
			// if (errorCount == 0)
			// {
			// status = true;
			// fileVersionNumber = handler.FileVersionNumber();
			// }
		}
		catch (SAXArrayStopParsingException se) {
			status = true;
			// fileVersionNumber = handler.FileVersionNumber();
		}
		catch (Exception e) {
			e.printStackTrace();
			status = false;
		}
		// delete parser;
		// XMLPlatformUtils::Terminate();

		return status;
	}

	/**
	 * Check if the data type matches what is in the file.
	 */
	public static boolean isFileType(String fileName, AffymetrixGuidType dataTypeId) {
		return (getDataTypeIdentifier(fileName) == dataTypeId);
	}

	/**
	 * Read just the first few entries to determine if this file is of the right type. Check the magic number, version
	 * number and data type identifier. If they all match then this is the right type of file.
	 */
	public static AffymetrixGuidType getDataTypeIdentifier(String fileName) {
		ArrayFileReader reader = new ArrayFileReader();
		ArrayData arrayData = new ArrayData();
		reader.read(fileName, arrayData, true);
		return arrayData.getDataTypeIdentifier();
	}
}
