package affymetrix.calvin.writers;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.Iterator;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import affymetrix.calvin.array.ArrayAttributes;
import affymetrix.calvin.array.ArrayData;
import affymetrix.calvin.array.ArrayId;
import affymetrix.calvin.array.ArrayMedia;
import affymetrix.calvin.array.CreateStep;
import affymetrix.calvin.array.PATAssignmentMethod;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.parameter.ParameterNameValueDefaultRequired;
import affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.ParameterValueType;
import affymetrix.calvin.parsers.SAXArrayHandler;
import affymetrix.calvin.utils.AffymetrixGuidType;
import affymetrix.calvin.utils.IOUtils;

public class ArrayFileWriter {

	/** The expected version number */
	public static final String ARRAY_SET_FILE_VERSION_NUMBER = "1.0";

	/** An identifier to the type of data stored in the file */
	protected AffymetrixGuidType dataTypeIdentifier = null;

	/**
	 * Initialize the class.
	 */
	public ArrayFileWriter() {
		dataTypeIdentifier = new AffymetrixGuidType(ArrayId.ARRAY_FILE_TYPE_IDENTIFIER);
	}

	/**
	 * Write the entire file, the header and body.
	 */
	public boolean write(String filename, ArrayData arrayData) {
		// // Initialize the XML4C2 system.
		// try {
		// XMLPlatformUtils.Initialize();
		// }
		// catch (XMLException e) {
		// return false;
		// }

		boolean status = false;
		BufferedOutputStream bufStream = null;
		try {
			// // Create a DOM implementation object and create the document type for it.
			// DOMImplementation impl = DOMImplementationRegistry.newInstance().getDOMImplementation("LS");
			// // DOMDocumentType* dt = impl.createDocumentType(ToXMLCh(ARRAY_FILE_ELEMENT), 0, ToXMLCh(ARRAY_FILE_DTD));
			// Document doc = impl.createDocument(null, null, null);
			// doc.setXmlStandalone(true);
			// // doc.appendChild(dt);
			//
			// // Create the serializer.
			// // DOMWriter theSerializer = impl.createDOMWriter();
			// // theSerializer.setEncoding(ToXMLCh(ARRAY_FILE_ENCODING));

			// Create an empty XML Document
			DocumentBuilder docBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
			Document doc = docBuilder.newDocument();

			// ArrayFile element
			Element arrayElement = createArrayElement(arrayData, doc /* , dataTypeIdentifier */);

			// PhysicalArrays element
			addPhysicalArrays(arrayData, doc, arrayElement);

			// UserAttributes element
			addUserAttributes(arrayData, doc, arrayElement);

			// Add the array element to the document.
			doc.appendChild(arrayElement);

			// Write the file
			TransformerFactory tFactory = TransformerFactory.newInstance();
			Transformer transformer = tFactory.newTransformer();
			transformer.setOutputProperty(OutputKeys.INDENT, "yes");
			transformer.setOutputProperty(OutputKeys.STANDALONE, "no");
			transformer.setOutputProperty(OutputKeys.ENCODING, "UTF-16");
			DOMSource source = new DOMSource(doc);
			bufStream = new BufferedOutputStream(new FileOutputStream(filename));
			transformer.transform(source, new StreamResult(bufStream));
			bufStream.flush();
			status = true;
		}
		catch (Exception e) {
			e.printStackTrace();
			status = false;
		}
		finally {
			if (bufStream != null) {
				try {
					bufStream.close();
				}
				catch (IOException e) {
				}
			}
		}
		return status;
	}

	/**
	 * Create the array element. This is the top level parent element.
	 */
	private Element createArrayElement(ArrayData arrayData, Document doc) {
		// ArrayFile element
		Element arrayElement = doc.createElement(SAXArrayHandler.ARRAY_FILE_ELEMENT);
		arrayElement
				.setAttribute(SAXArrayHandler.ARRAY_FILE_ELEMENT_TYPE_ATTRIBUTE, ArrayId.ARRAY_SET_FILE_TYPE_IDENTIFIER);
		arrayElement.setAttribute(SAXArrayHandler.ARRAY_FILE_ELEMENT_VERSION_ATTRIBUTE, ARRAY_SET_FILE_VERSION_NUMBER);
		arrayElement.setAttribute(SAXArrayHandler.ARRAY_FILE_ELEMENT_ID_ATTRIBUTE, guidToXmlString(arrayData
				.getArraySetFileIdentifier()));
		arrayElement.setAttribute(SAXArrayHandler.ARRAY_FILE_ELEMENT_PROJECT_ATTRIBUTE, arrayData.getInitialProject());
		arrayElement.setAttribute(SAXArrayHandler.ARRAY_FILE_ELEMENT_CREATE_DATE_TIME_ATTRIBUTE, arrayData
				.getCreationDateTime());
		arrayElement.setAttribute(SAXArrayHandler.ARRAY_FILE_ELEMENT_CREATED_BY_ATTRIBUTE, arrayData.getCreatedBy());

		String step = CreateStep.toString(arrayData.getCreatedStep());
		arrayElement.setAttribute(SAXArrayHandler.ARRAY_FILE_ELEMENT_CREATED_STEP_ATTRIBUTE, step);
		return arrayElement;
	}

	/**
	 * Add the physical arrays to the document if they exist.
	 */
	private void addPhysicalArrays(ArrayData arrayData, Document doc, Element arrayElement) {
		if (arrayData.getPhysicalArraysAttributes().size() > 0) {
			Element physicalArraysElement = doc.createElement(SAXArrayHandler.PHYSICAL_ARRAYS_ELEMENT);
			Iterator<ArrayAttributes> it = arrayData.getPhysicalArraysAttributes().iterator();
			while (it.hasNext()) {
				ArrayAttributes att = it.next();
				Element physicalArrayElement = doc.createElement(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT);
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_ARRAY_TYPE_ATTRIBUTE,
						ArrayId.ARRAY_TYPE_IDENTIFIER);
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_ID_ATTRIBUTE, guidToXmlString(att
						.getIdentifier()));
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_NAME_ATTRIBUTE, att.getArrayName());
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_BARCODE_ATTRIBUTE, att
						.getArrayBarcode());
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_TYPE_ATTRIBUTE, ArrayMedia
						.toString(att.getMedia()));
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_ROW_ATTRIBUTE, String.valueOf(att
						.getMediaRow()));
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_COL_ATTRIBUTE, String.valueOf(att
						.getMediaCol()));
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_MEDIA_FILE_NAME_ATTRIBUTE, att
						.getMediaFileName());
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_MEDIA_FILE_GUID_ATTRIBUTE,
						guidToXmlString(att.getMediaFileGUID()));
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_LIB_PACKAGE_NAME_ATTRIBUTE, att
						.getLibraryPackageName());
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_MASTERFILE_GUID_ATTRIBUTE,
						guidToXmlString(att.getMasterFileId()));
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_CREATED_BY_ATTRIBUTE, att
						.getCreatedBy());
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_CREATION_DATE_ATTRIBUTE, att
						.getCreationDateTime());
				physicalArrayElement.setAttribute(SAXArrayHandler.ARRAY_FILE_ELEMENT_CREATED_STEP_ATTRIBUTE, CreateStep
						.toString(att.getCreatedStep()));
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_MASTERFILE_ATTRIBUTE, att
						.getMasterFile());
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_PAT_ASSIGNMENT_ATTRIBUTE,
						PATAssignmentMethod.toString(att.getPATAssignment()));
				physicalArrayElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ELEMENT_COMMENT_ATTRIBUTE, att.getComment());

				Iterator<ParameterNameValue> attrIt = att.getAttributes().iterator();
				while (attrIt.hasNext()) {
					ParameterNameValue param = attrIt.next();
					Element paramElement = doc.createElement(SAXArrayHandler.PHYSICAL_ARRAY_ATTRIBUTE_ELEMENT);
					paramElement.setAttribute(SAXArrayHandler.PHYSICAL_ARRAY_ATTRIBUTE_ELEMENT_NAME_ATTRIBUTE, param.getName());
					paramElement.setTextContent(param.getValueText());
					physicalArrayElement.appendChild(paramElement);
				}
				physicalArraysElement.appendChild(physicalArrayElement);
			}
			arrayElement.appendChild(physicalArraysElement);
		}
	}

	/**
	 * Add the user attributes to the document if they exist.
	 */
	private void addUserAttributes(ArrayData arrayData, Document doc, Element arrayElement) {
		if (arrayData.getUserAttributes().size() > 0) {
			Element userAttributesElement = doc.createElement(SAXArrayHandler.USER_ATTRIBUTES_ELEMENT);

			Iterator<ParameterNameValueDefaultRequired> paramIt = arrayData.getUserAttributes().iterator();
			while (paramIt.hasNext()) {
				ParameterNameValueDefaultRequired param = paramIt.next();
				Element paramElement = doc.createElement(SAXArrayHandler.USER_ATTRIBUTES_ATTRIBUTE_ELEMENT);
				paramElement.setAttribute(SAXArrayHandler.USER_ATTRIBUTES_ATTRIBUTE_ELEMENT_NAME_ATTRIBUTE, param.getName());
				paramElement.setAttribute(SAXArrayHandler.USER_ATTRIBUTES_ATTRIBUTE_ELEMENT_TYPE_ATTRIBUTE,
						ParameterNameValueDefaultRequired.ParameterValueTypeToString(param.getValueType()));

				//required attribute
				paramElement.setAttribute(SAXArrayHandler.USER_ATTRIBUTES_ATTRIBUTE_ELEMENT_REQUIRED_ATTRIBUTE, Boolean
						.toString(param.getRequiredFlag()));

				if (param.getHasDefault()) {
					paramElement.setAttribute(SAXArrayHandler.USER_ATTRIBUTES_ATTRIBUTE_ELEMENT_DEFAULT_ATTRIBUTE, String
							.valueOf(param.getDefaultValueFloat()));
				}
				if (param.getValueType() == ParameterValueType.ControlMultiParameterType) {
					Iterator<String> valIt = param.getControlMultiValues().iterator();
					while (valIt.hasNext()) {
						Element valueElement = doc.createElement(SAXArrayHandler.USER_ATTRIBUTES_ATTRIBUTE_VALUE_ELEMENT);
						valueElement.setTextContent(valIt.next());
						paramElement.appendChild(valueElement);
					}
				}
				else {
					Element valueElement = doc.createElement(SAXArrayHandler.USER_ATTRIBUTES_ATTRIBUTE_VALUE_ELEMENT);
					valueElement.setTextContent(param.toString());
					paramElement.appendChild(valueElement);
				}
				if (param.getControlledVocabulary().size() > 0) {
					Iterator<String> controlIt = param.getControlledVocabulary().iterator();
					while (controlIt.hasNext()) {
						Element controlElement = doc.createElement(SAXArrayHandler.USER_ATTRIBUTES_ATTRIBUTE_CONTROL_ELEMENT);
						controlElement.setAttribute(SAXArrayHandler.USER_ATTRIBUTES_ATTRIBUTE_CONTROL_ELEMENT_VALUE_ATTRIBUTE,
								controlIt.next());
						paramElement.appendChild(controlElement);
					}
				}
				userAttributesElement.appendChild(paramElement);
			}
			arrayElement.appendChild(userAttributesElement);
		}
	}

	private String guidToXmlString(AffymetrixGuidType guid) {
		byte[] g = guid.getGuid();
		for (int i = 0; i < g.length; i++) {
			if (g[i] == 0) {
				byte[] buf = new byte[i];
				System.arraycopy(g, 0, buf, 0, buf.length);
				try {
					return new String(buf, IOUtils.ASCII_CHARSET);
				}
				catch (UnsupportedEncodingException e) {
					e.printStackTrace();
				}
			}
		}
		return null;
	}
}
