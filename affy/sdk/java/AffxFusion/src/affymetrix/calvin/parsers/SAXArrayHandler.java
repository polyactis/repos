package affymetrix.calvin.parsers;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import affymetrix.calvin.array.ArrayAttributes;
import affymetrix.calvin.array.ArrayData;
import affymetrix.calvin.array.ArrayMedia;
import affymetrix.calvin.array.CreateStep;
import affymetrix.calvin.array.PATAssignmentMethod;
import affymetrix.calvin.exception.SAXArrayStopParsingException;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.parameter.ParameterNameValueDefaultRequired;
import affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.ParameterValueType;
import affymetrix.calvin.utils.AffymetrixGuidType;
import affymetrix.calvin.utils.IOUtils;

public class SAXArrayHandler extends DefaultHandler {
	/* ! The name of the DTD for the array file. */
	public static final String ARRAY_FILE_DTD = "ArraySetAndTemplateFile.dtd";

	/* ! The encoding to use for array files. */
	public static final String ARRAY_FILE_ENCODING = "UTF-16";

	/* ! The name of the element that contains the array file id, type and version data. */
	public static final String ARRAY_FILE_ELEMENT = "ArraySetFile";

	/* ! The attribute name of the ID field for the array file. */
	public static final String ARRAY_FILE_ELEMENT_ID_ATTRIBUTE = "GUID";

	/* ! The attribute name of the type field for the array file. */
	public static final String ARRAY_FILE_ELEMENT_TYPE_ATTRIBUTE = "Type";

	/* ! The attribute name of the version field for the array file. */
	public static final String ARRAY_FILE_ELEMENT_VERSION_ATTRIBUTE = "Version";

	/* ! The attribute name of the original project field for the array file. */
	public static final String ARRAY_FILE_ELEMENT_PROJECT_ATTRIBUTE = "OriginalProjectName";

	/* ! The attribute name of the creation date field for the array file. */
	public static final String ARRAY_FILE_ELEMENT_CREATE_DATE_TIME_ATTRIBUTE = "CreatedDateTime";

	/* ! The attribute name of the create by field for the array file. */
	public static final String ARRAY_FILE_ELEMENT_CREATED_BY_ATTRIBUTE = "CreatedBy";

	/* ! The attribute name of the create step field for the array file. */
	public static final String ARRAY_FILE_ELEMENT_CREATED_STEP_ATTRIBUTE = "CreatedStep";

	/* ! The name of the element that contains the list of physical arrays. */
	public static final String PHYSICAL_ARRAYS_ELEMENT = "PhysicalArrays";

	/* ! The name of the element that contains the attributes of a single physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT = "PhysicalArray";

	/* ! The attribute name of the type field for a physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT_ARRAY_TYPE_ATTRIBUTE = "Type";

	/* ! The attribute name of the ID field for a physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT_ID_ATTRIBUTE = "GUID";

	/* ! The attribute name of the array name field for a physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT_NAME_ATTRIBUTE = "ArrayName";

	/* ! The attribute name of the barcode field for a physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT_BARCODE_ATTRIBUTE = "AffyBarcode";

	/* ! The attribute name of the type field for a physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT_TYPE_ATTRIBUTE = "MediaType";

	/* ! The attribute name of the row field for a physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT_ROW_ATTRIBUTE = "MediaRow";

	/* ! The attribute name of the col field for a physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT_COL_ATTRIBUTE = "MediaCol";

	/* ! The attribute name of the media file name. */
	public static final String PHYSICAL_ARRAY_ELEMENT_MEDIA_FILE_NAME_ATTRIBUTE = "MediaFileName";

	/* ! The attribute name of the media file guid. */
	public static final String PHYSICAL_ARRAY_ELEMENT_MEDIA_FILE_GUID_ATTRIBUTE = "MediaFileGUID";

	/* ! The attribute name of the library file package name. */
	public static final String PHYSICAL_ARRAY_ELEMENT_LIB_PACKAGE_NAME_ATTRIBUTE = "LibraryPackageName";

	/* ! The attribute name of the master file field for a physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT_MASTERFILE_ATTRIBUTE = "MasterFileName";

	/* ! The attribute name of the master file guid field for a physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT_MASTERFILE_GUID_ATTRIBUTE = "MasterFileGUID";

	/* ! The attribute name of the pat assignment field for a physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT_PAT_ASSIGNMENT_ATTRIBUTE = "PATAssignmentMethod";

	/* ! The attribute name of the creation date field for a physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT_CREATION_DATE_ATTRIBUTE = "CreatedDateTime";

	/* ! The attribute name of the creation user field for a physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT_CREATED_BY_ATTRIBUTE = "CreatedBy";

	/* ! The attribute name of the comment field for a physical array. */
	public static final String PHYSICAL_ARRAY_ELEMENT_COMMENT_ATTRIBUTE = "Comment";

	/* ! The name of the element that contains a single attribute of a single physical array. */
	public static final String PHYSICAL_ARRAY_ATTRIBUTE_ELEMENT = "ArrayAttribute";

	/* ! The attribute name of the name field for a physical array attribute. */
	public static final String PHYSICAL_ARRAY_ATTRIBUTE_ELEMENT_NAME_ATTRIBUTE = "Name";

	/* ! The attribute name of the value field for a physical array attribute. */
	public static final String PHYSICAL_ARRAY_ATTRIBUTE_ELEMENT_VALUE_ATTRIBUTE = "Value";

	/* ! The name of the element that contains the user attributes. */
	public static final String USER_ATTRIBUTES_ELEMENT = "UserAttributes";

	/* ! The name of the element that contains a single user attribute. */
	public static final String USER_ATTRIBUTES_ATTRIBUTE_ELEMENT = "UserAttribute";

	/* ! The name of the element that contains a single user attribute. */
	public static final String USER_ATTRIBUTES_ATTRIBUTE_VALUE_ELEMENT = "UserAttributeValue";

	/* ! The attribute name of the name field for a user attribute. */
	public static final String USER_ATTRIBUTES_ATTRIBUTE_ELEMENT_NAME_ATTRIBUTE = "Name";

	/* ! The attribute name of the type field for a user attribute. */
	public static final String USER_ATTRIBUTES_ATTRIBUTE_ELEMENT_TYPE_ATTRIBUTE = "Type";

	/* ! The attribute name of the default value field for a user attribute. */
	public static final String USER_ATTRIBUTES_ATTRIBUTE_ELEMENT_DEFAULT_ATTRIBUTE = "DefaultValue";

	/* ! The attribute name of the required field for a user attribute. */
	public static final String USER_ATTRIBUTES_ATTRIBUTE_ELEMENT_REQUIRED_ATTRIBUTE = "Required";

	/* ! The attribute name of the value field for a user attribute. */
	public static final String USER_ATTRIBUTES_ATTRIBUTE_ELEMENT_VALUE_ATTRIBUTE = "Value";

	/* ! The name of the element that contains a controlled vocabulary. */
	public static final String USER_ATTRIBUTES_ATTRIBUTE_CONTROL_ELEMENT = "Control";

	/* ! The attribute name of the value field for a controlled vocabulary. */
	public static final String USER_ATTRIBUTES_ATTRIBUTE_CONTROL_ELEMENT_VALUE_ATTRIBUTE = "Value";

	/* ! Enumerants to hold the elements in an array file. */
	public enum ArrayFileElements {
		ARRAY_FILE, PHYSICAL_ARRAYS, PHYSICAL_ARRAY, PHYSICAL_ARRAY_ATTRIBUTES, USER_ATTRIBUTES, USER_ATTRIBUTES_ATTRIBUTE, USER_ATTRIBUTES_ATTRIBUTE_VALUE, USER_ATTRIBUTES_ATTRIBUTE_CONTROL
	}

	/* ! A pointer to the array object. */
	private ArrayData arrayData = null;

	/* ! A flag used to indicate that the header line should only be read. */
	private boolean readHeaderOnly;

	/* ! The parent element that is currently being processed. */
	private ArrayFileElements currentElement = null;

	/* ! The files version number. */
	// private String fileVersionNumber = null;
	/*
	 * Store the array data and set the starting element to the head.
	 */
	public SAXArrayHandler(ArrayData data, boolean hdrOnly) {
		arrayData = data;
		currentElement = ArrayFileElements.ARRAY_FILE;
		readHeaderOnly = hdrOnly;
	}

	/*
	 * Back up the current element.
	 */
	public void endElement(String uri, String localName, String qName) throws SAXException {
		moveCurrentElementBack(localName);
	}

	public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException {
		if (!moveCurrentElementForward(localName)) {
			return;
		}
		Map<String, String> atts = new HashMap<String, String>();
		int len = attributes.getLength();
		for (int i = 0; i < len; i++) {
			atts.put(attributes.getLocalName(i), attributes.getValue(i));
		}
		switch (currentElement) {
		case ARRAY_FILE:
			storeArrayFileAttributes(atts);
			break;

		case PHYSICAL_ARRAY:
			storePhysicalArrayElementAttributes(atts);
			break;

		case PHYSICAL_ARRAY_ATTRIBUTES:
			storePhysicalArrayAttribute(atts);
			break;

		case USER_ATTRIBUTES_ATTRIBUTE:
			storeUserAttribute(atts);
			break;

		case USER_ATTRIBUTES_ATTRIBUTE_CONTROL:
			storeUserAttributeControl(atts);
			break;

		default:
			break;
		}
		if (readHeaderOnly) {
			throw new SAXArrayStopParsingException();
		}
	}

	/*
	 * Back up the current element.
	 */
	private void moveCurrentElementBack(String name) {
		if (IOUtils.isNullOrEmpty(name)) {
			return;
		}
		if (name.equals(PHYSICAL_ARRAY_ATTRIBUTE_ELEMENT)
				&& (currentElement == ArrayFileElements.PHYSICAL_ARRAY_ATTRIBUTES)) {
			currentElement = ArrayFileElements.PHYSICAL_ARRAY;
		}
		else if (name.equals(PHYSICAL_ARRAY_ELEMENT)) {
			currentElement = ArrayFileElements.PHYSICAL_ARRAYS;
		}
		else if (name.equals(PHYSICAL_ARRAYS_ELEMENT)) {
			currentElement = ArrayFileElements.ARRAY_FILE;
		}
		else if (name.equals(USER_ATTRIBUTES_ATTRIBUTE_CONTROL_ELEMENT)
				&& (currentElement == ArrayFileElements.USER_ATTRIBUTES_ATTRIBUTE_CONTROL)) {
			currentElement = ArrayFileElements.USER_ATTRIBUTES_ATTRIBUTE;
		}
		else if (name.equals(USER_ATTRIBUTES_ATTRIBUTE_VALUE_ELEMENT)
				&& (currentElement == ArrayFileElements.USER_ATTRIBUTES_ATTRIBUTE_VALUE)) {
			currentElement = ArrayFileElements.USER_ATTRIBUTES_ATTRIBUTE;
		}
		else if (name.equals(USER_ATTRIBUTES_ATTRIBUTE_ELEMENT)
				&& (currentElement == ArrayFileElements.USER_ATTRIBUTES_ATTRIBUTE_VALUE)) {
			currentElement = ArrayFileElements.USER_ATTRIBUTES;
		}
		else if (name.equals(USER_ATTRIBUTES_ELEMENT)) {
			currentElement = ArrayFileElements.ARRAY_FILE;
		}
	}

	/*
	 * Advance the current element.
	 */
	private boolean moveCurrentElementForward(String name) {
		if (IOUtils.isNullOrEmpty(name)) {
			return false;
		}
		if (name.equals(ARRAY_FILE_ELEMENT)) {
			currentElement = ArrayFileElements.ARRAY_FILE;
		}
		else if (name.equals(PHYSICAL_ARRAYS_ELEMENT)) {
			currentElement = ArrayFileElements.PHYSICAL_ARRAYS;
		}
		else if (name.equals(PHYSICAL_ARRAY_ELEMENT)) {
			currentElement = ArrayFileElements.PHYSICAL_ARRAY;
		}
		else if (name.equals(PHYSICAL_ARRAY_ATTRIBUTE_ELEMENT) && (currentElement == ArrayFileElements.PHYSICAL_ARRAY)) {
			currentElement = ArrayFileElements.PHYSICAL_ARRAY_ATTRIBUTES;
		}
		else if (name.equals(USER_ATTRIBUTES_ELEMENT)) {
			currentElement = ArrayFileElements.USER_ATTRIBUTES;
		}
		else if (name.equals(USER_ATTRIBUTES_ATTRIBUTE_ELEMENT) && (currentElement == ArrayFileElements.USER_ATTRIBUTES)) {
			currentElement = ArrayFileElements.USER_ATTRIBUTES_ATTRIBUTE;
		}
		else if (name.equals(USER_ATTRIBUTES_ATTRIBUTE_VALUE_ELEMENT)
				&& (currentElement == ArrayFileElements.USER_ATTRIBUTES_ATTRIBUTE)) {
			currentElement = ArrayFileElements.USER_ATTRIBUTES_ATTRIBUTE_VALUE;
		}
		else if (name.equals(USER_ATTRIBUTES_ATTRIBUTE_CONTROL_ELEMENT)
				&& (currentElement == ArrayFileElements.USER_ATTRIBUTES_ATTRIBUTE)) {
			currentElement = ArrayFileElements.USER_ATTRIBUTES_ATTRIBUTE_CONTROL;
		}
		else if (name.equals(ARRAY_FILE_ELEMENT) || name.equals(PHYSICAL_ARRAYS_ELEMENT)
				|| name.equals(PHYSICAL_ARRAY_ELEMENT) || name.equals(PHYSICAL_ARRAY_ATTRIBUTE_ELEMENT)
				|| name.equals(USER_ATTRIBUTES_ELEMENT) || name.equals(USER_ATTRIBUTES_ATTRIBUTE_ELEMENT)
				|| name.equals(USER_ATTRIBUTES_ATTRIBUTE_VALUE_ELEMENT)
				|| name.equals(USER_ATTRIBUTES_ATTRIBUTE_CONTROL_ELEMENT)) {
		}
		else {
			return false;
		}
		return true;
	}

	@Override
	public void characters(char[] ch, int start, int length) throws SAXException {
		// String str = XMLChToWString(chars);
		String str = String.valueOf(ch, start, length);
		if (currentElement == ArrayFileElements.PHYSICAL_ARRAY_ATTRIBUTES) {
			int last = arrayData.getPhysicalArraysAttributes().size() - 1;
			List<ParameterNameValue> arrayAtts = arrayData.getPhysicalArraysAttributes().get(last).getAttributes();
			last = arrayAtts.size() - 1;
			arrayAtts.get(last).setValueText(str);
		}
		else if (currentElement == ArrayFileElements.USER_ATTRIBUTES_ATTRIBUTE_VALUE) {
			int last = arrayData.getUserAttributes().size() - 1;
			ParameterNameValueDefaultRequired param = arrayData.getUserAttributes().get(last);
			switch (param.getValueType()) {
			case IntegerParameterType:
				param.setValueInt32(Integer.parseInt(str));
				break;

			case FloatParameterType:
				param.setValueFloat(Float.parseFloat(str));
				break;

			case TextParameterType:
			case DateParameterType:
			case TimeParameterType:
			case DateTimeParameterType:
			case ControlSingleParameterType:
				param.setValueText(str);
				break;

			case ControlMultiParameterType:
				param.getControlMultiValues().add(str);
				break;

			default:
				break;
			}
		}
	}

	/*
	 * Store the user attribute.
	 */
	private void storeUserAttribute(Map<String, String> attributes) {
		ParameterNameValueDefaultRequired param = new ParameterNameValueDefaultRequired();
		param.setName(attributes.get(USER_ATTRIBUTES_ATTRIBUTE_ELEMENT_NAME_ATTRIBUTE));
		String value = attributes.get(USER_ATTRIBUTES_ATTRIBUTE_ELEMENT_TYPE_ATTRIBUTE);
		String defValue = attributes.get(USER_ATTRIBUTES_ATTRIBUTE_ELEMENT_DEFAULT_ATTRIBUTE);

		if (value.equals(ParameterNameValueDefaultRequired
				.ParameterValueTypeToString(ParameterValueType.IntegerParameterType))) {
			param.setValueInt32(0);
			param.setValueType(ParameterValueType.IntegerParameterType);
			if (defValue.length() > 0) {
				param.setDefaultValueInt32(Integer.parseInt(defValue));
			}
		}
		else if (value.equals(ParameterNameValueDefaultRequired
				.ParameterValueTypeToString(ParameterValueType.FloatParameterType))) {
			param.setValueFloat(0.0f);
			param.setValueType(ParameterValueType.FloatParameterType);
			if (!IOUtils.isNullOrEmpty(defValue)) {
				param.setDefaultValueFloat(Float.parseFloat(defValue));
			}
		}
		else {
			if (value.equals(ParameterNameValueDefaultRequired
					.ParameterValueTypeToString(ParameterValueType.TextParameterType))) {
				param.setValueType(ParameterValueType.TextParameterType);
			}
			else if (value.equals(ParameterNameValueDefaultRequired
					.ParameterValueTypeToString(ParameterNameValueDefaultRequired.ParameterValueType.DateParameterType))) {
				param.setValueType(ParameterValueType.DateParameterType);
			}
			else if (value.equals(ParameterNameValueDefaultRequired
					.ParameterValueTypeToString(ParameterNameValueDefaultRequired.ParameterValueType.TimeParameterType))) {
				param.setValueType(ParameterValueType.TimeParameterType);
			}
			else if (value.equals(ParameterNameValueDefaultRequired
					.ParameterValueTypeToString(ParameterNameValueDefaultRequired.ParameterValueType.DateTimeParameterType))) {
				param.setValueType(ParameterValueType.DateTimeParameterType);
			}
			else if (value.equals(ParameterNameValueDefaultRequired
					.ParameterValueTypeToString(ParameterNameValueDefaultRequired.ParameterValueType.ControlSingleParameterType))) {
				param.setValueType(ParameterValueType.ControlSingleParameterType);
			}
			else if (value.equals(ParameterNameValueDefaultRequired
					.ParameterValueTypeToString(ParameterNameValueDefaultRequired.ParameterValueType.ControlMultiParameterType))) {
				param.setValueType(ParameterValueType.ControlMultiParameterType);
			}
			param.setValueText(IOUtils.EMPTY);
			if (!IOUtils.isNullOrEmpty(defValue)) {
				param.setDefaultValueText(defValue);
			}
		}
		value = attributes.get(USER_ATTRIBUTES_ATTRIBUTE_ELEMENT_REQUIRED_ATTRIBUTE);
		if (value != null && value.equalsIgnoreCase(String.valueOf(true))) {
			param.setRequiredFlag(true);
		}
		else {
			param.setRequiredFlag(false);
		}
		arrayData.getUserAttributes().add(param);
	}

	/*
	 * Store the control value for a user attribute.
	 */
	private void storeUserAttributeControl(Map<String, String> attributes) {
		List<ParameterNameValueDefaultRequired> ua = arrayData.getUserAttributes();
		int last = ua.size() - 1;
		ParameterNameValueDefaultRequired param = ua.get(last);
		param.getControlledVocabulary().add(attributes.get(USER_ATTRIBUTES_ATTRIBUTE_CONTROL_ELEMENT_VALUE_ATTRIBUTE));
	}

	/*
	 * Store the array attribute.
	 */
	void storePhysicalArrayAttribute(Map<String, String> attributes) {
		int last = arrayData.getPhysicalArraysAttributes().size() - 1;
		List<ParameterNameValue> arrAttrs = arrayData.getPhysicalArraysAttributes().get(last).getAttributes();
		ParameterNameValue p = new ParameterNameValue();
		p.setName(attributes.get(PHYSICAL_ARRAY_ATTRIBUTE_ELEMENT_NAME_ATTRIBUTE));
		// p.setValueText(attributes.get(PHYSICAL_ARRAY_ATTRIBUTE_ELEMENT_VALUE_ATTRIBUTE));
		arrAttrs.add(p);
	}

	/*
	 * Create a new entry in the physical array and store the ID.
	 */
	void storePhysicalArrayElementAttributes(Map<String, String> attributes) {

		ArrayAttributes arrAttrs = new ArrayAttributes();
		arrAttrs.setIdentifier(new AffymetrixGuidType(attributes.get(PHYSICAL_ARRAY_ELEMENT_ID_ATTRIBUTE)));

		String str = attributes.get(PHYSICAL_ARRAY_ELEMENT_ROW_ATTRIBUTE);
		if (!IOUtils.isNullOrEmpty(str)) {
			arrAttrs.setMediaRow(Integer.parseInt(str));
		}

		str = attributes.get(PHYSICAL_ARRAY_ELEMENT_COL_ATTRIBUTE);
		if (!IOUtils.isNullOrEmpty(str)) {
			arrAttrs.setMediaCol(Integer.parseInt(str));
		}

		str = attributes.get(PHYSICAL_ARRAY_ELEMENT_MEDIA_FILE_NAME_ATTRIBUTE);
		if (!IOUtils.isNullOrEmpty(str)) {
			arrAttrs.setMediaFileName(str);
		}

		str = attributes.get(PHYSICAL_ARRAY_ELEMENT_MEDIA_FILE_GUID_ATTRIBUTE);
		if (!IOUtils.isNullOrEmpty(str)) {
			arrAttrs.setMediaFileGUID(new AffymetrixGuidType(str));
		}

		str = attributes.get(PHYSICAL_ARRAY_ELEMENT_LIB_PACKAGE_NAME_ATTRIBUTE);
		if (!IOUtils.isNullOrEmpty(str)) {
			arrAttrs.setLibraryPackageName(str);
		}

		str = attributes.get(PHYSICAL_ARRAY_ELEMENT_NAME_ATTRIBUTE);
		if (!IOUtils.isNullOrEmpty(str)) {
			arrAttrs.setArrayName(str);
		}

		str = attributes.get(PHYSICAL_ARRAY_ELEMENT_BARCODE_ATTRIBUTE);
		if (!IOUtils.isNullOrEmpty(str)) {
			arrAttrs.setArrayBarcode(str);
		}

		str = attributes.get(PHYSICAL_ARRAY_ELEMENT_TYPE_ATTRIBUTE);
		arrAttrs.setMedia(ArrayMedia.toArrayMediaType(str));

		str = attributes.get(PHYSICAL_ARRAY_ELEMENT_MASTERFILE_ATTRIBUTE);
		if (!IOUtils.isNullOrEmpty(str)) {
			arrAttrs.setMasterFile(str);
			;
		}
		str = attributes.get(PHYSICAL_ARRAY_ELEMENT_PAT_ASSIGNMENT_ATTRIBUTE);
		if (!IOUtils.isNullOrEmpty(str)) {
			arrAttrs.setPATAssignment(PATAssignmentMethod.toPATAssignmentMethodType(str));
		}
		str = attributes.get(PHYSICAL_ARRAY_ELEMENT_MASTERFILE_GUID_ATTRIBUTE);
		if (!IOUtils.isNullOrEmpty(str)) {
			arrAttrs.setMasterFileId(new AffymetrixGuidType(str));
		}
		arrAttrs.setCreationDateTime(attributes.get(PHYSICAL_ARRAY_ELEMENT_CREATION_DATE_ATTRIBUTE));
		arrAttrs.setCreatedBy(attributes.get(PHYSICAL_ARRAY_ELEMENT_CREATED_BY_ATTRIBUTE));
		arrAttrs.setComment(attributes.get(PHYSICAL_ARRAY_ELEMENT_COMMENT_ATTRIBUTE));

		str = attributes.get(ARRAY_FILE_ELEMENT_CREATED_STEP_ATTRIBUTE);
		if (!IOUtils.isNullOrEmpty(str)) {
			arrAttrs.setCreatedStep(CreateStep.toCreateStepType(str));
		}

		arrayData.getPhysicalArraysAttributes().add(arrAttrs);
	}

	/*
	 * Store the attributes to the array data and member variables.
	 */
	void storeArrayFileAttributes(Map<String, String> attributes) {
		String str = attributes.get(ARRAY_FILE_ELEMENT_TYPE_ATTRIBUTE);
		if (!IOUtils.isNullOrEmpty(str)) {
			arrayData.setDataTypeIdentifier(new AffymetrixGuidType(str));
		}

		// fileVersionNumber = attributes.get(ARRAY_FILE_ELEMENT_VERSION_ATTRIBUTE);

		arrayData.setArraySetFileIdentifier(new AffymetrixGuidType(attributes.get(ARRAY_FILE_ELEMENT_ID_ATTRIBUTE)));
		arrayData.setInitialProject(attributes.get(ARRAY_FILE_ELEMENT_PROJECT_ATTRIBUTE));
		arrayData.setCreationDateTime(attributes.get(ARRAY_FILE_ELEMENT_CREATE_DATE_TIME_ATTRIBUTE));
		arrayData.setCreatedBy(attributes.get(ARRAY_FILE_ELEMENT_CREATED_BY_ATTRIBUTE));
		str = attributes.get(ARRAY_FILE_ELEMENT_CREATED_STEP_ATTRIBUTE);
		arrayData.setCreatedStep(CreateStep.toCreateStepType(str));
	}
}
