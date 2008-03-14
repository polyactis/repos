package affymetrix.calvin.writers;

import java.util.Iterator;

import junit.framework.TestCase;
import affymetrix.calvin.array.ArrayAttributes;
import affymetrix.calvin.array.ArrayData;
import affymetrix.calvin.array.ArrayId;
import affymetrix.calvin.array.ArrayMedia.ArrayMediaType;
import affymetrix.calvin.array.CreateStep.CreateStepType;
import affymetrix.calvin.array.PATAssignmentMethod.PATAssignmentMethodType;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.parameter.ParameterNameValueDefaultRequired;
import affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.ParameterValueType;
import affymetrix.calvin.parsers.ArrayFileReader;
import affymetrix.calvin.utils.AffymetrixGuidType;

public class ArrayFileWriterTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testWriteFile() {
		ArrayData array1 = new ArrayData();

		array1.setCreatedBy("ljevon");
		array1.setCreatedStep(CreateStepType.ScanningStep);
		array1.setCreationDateTime("datetime");
		array1.setArraySetFileIdentifier(new AffymetrixGuidType("id"));
		array1.setInitialProject("proj1");

		ArrayAttributes attr = new ArrayAttributes();
		attr.setIdentifier(new AffymetrixGuidType("array1-id"));
		attr.setArrayBarcode("123");
		attr.setArrayName("name");
		attr.setComment("comment goes here");
		attr.setCreatedBy("me");
		attr.setCreatedStep(CreateStepType.JobOrderServerStep);
		attr.setCreationDateTime("now");
		attr.setMasterFile("master");
		attr.setMasterFileId(new AffymetrixGuidType("456"));
		attr.setMedia(ArrayMediaType.PlateOrStripMedia);
		attr.setMediaRow(1);
		attr.setMediaCol(2);
		attr.setPATAssignment(PATAssignmentMethodType.AffyBarcodeAssignment);
		attr.setMediaFileGUID(new AffymetrixGuidType("mediaguid"));
		attr.setMediaFileName("mediafile");
		attr.setLibraryPackageName("libpackage");
		array1.getPhysicalArraysAttributes().add(attr);

		ParameterNameValue param1 = new ParameterNameValue();
		param1.setName("array-att-name-1");
		param1.setValueText("array-att-value-1");
		attr.getAttributes().add(param1);

		ParameterNameValue param2 = new ParameterNameValue();
		param2.setName("array-att-name-2");
		param2.setValueText("array-att-value-2");
		attr.getAttributes().add(param2);

		ParameterNameValueDefaultRequired vparam1 = new ParameterNameValueDefaultRequired();
		vparam1.setName("user-att-name-1");
		vparam1.setValueText("user-att-value-1");
		vparam1.setRequiredFlag(false);
		vparam1.setHasDefault(false);
		vparam1.setValueType(ParameterValueType.TextParameterType);
		vparam1.getControlledVocabulary().clear();
		array1.getUserAttributes().add(vparam1);

		ParameterNameValueDefaultRequired vparam2 = new ParameterNameValueDefaultRequired();
		vparam2.setName("user-att-name-2");
		vparam2.setValueFloat(1f);
		vparam2.setDefaultValueFloat(2f);
		vparam2.setHasDefault(true);
		vparam2.setRequiredFlag(true);
		vparam2.setValueType(ParameterValueType.FloatParameterType);
		vparam2.getControlledVocabulary().clear();
		array1.getUserAttributes().add(vparam2);

		ParameterNameValueDefaultRequired vparam3 = new ParameterNameValueDefaultRequired();
		vparam3.setName("user-att-name-3");
		vparam3.getControlMultiValues().add("one");
		vparam3.getControlMultiValues().add("two");
		vparam3.setRequiredFlag(false);
		vparam3.setHasDefault(false);
		vparam3.setValueType(ParameterValueType.ControlMultiParameterType);
		vparam3.getControlledVocabulary().add("one");
		vparam3.getControlledVocabulary().add("two");
		vparam3.getControlledVocabulary().add("three");
		array1.getUserAttributes().add(vparam3);

		// Write the file.
		String testFile = "test_file_array.xml";
		ArrayFileWriter writer = new ArrayFileWriter();
		writer.write(testFile, array1);

		// Read the file.
		ArrayFileReader reader = new ArrayFileReader();
		ArrayData array2 = new ArrayData();
		reader.read(testFile, array2);

		assertEquals(array2.getCreatedBy(), "ljevon");
		assertEquals(array2.getCreatedStep(), CreateStepType.ScanningStep);
		assertEquals(array2.getCreationDateTime(), "datetime");
		assertEquals(array2.getArraySetFileIdentifier().toString(), new AffymetrixGuidType("id").toString());
		assertEquals(array2.getInitialProject(), "proj1");
		assertEquals(array2.getDataTypeIdentifier().toString(), new AffymetrixGuidType(
				ArrayId.ARRAY_SET_FILE_TYPE_IDENTIFIER).toString());

		// Array attributes.
		attr = array2.getPhysicalArraysAttributes().get(0);
		assertEquals(array2.getPhysicalArraysAttributes().size(), 1);
		assertEquals(attr.getIdentifier().toString(), new AffymetrixGuidType("array1-id").toString());
		assertEquals(attr.getArrayBarcode(), "123");
		assertEquals(attr.getArrayName(), "name");
		assertEquals(attr.getComment(), "comment goes here");
		assertEquals(attr.getCreatedBy(), "me");
		assertEquals(attr.getCreatedStep(), CreateStepType.JobOrderServerStep);
		assertEquals(attr.getCreationDateTime(), "now");
		assertEquals(attr.getMasterFile(), "master");
		assertEquals(attr.getMasterFileId().toString(), new AffymetrixGuidType("456").toString());
		assertEquals(attr.getMedia(), ArrayMediaType.PlateOrStripMedia);
		assertEquals(attr.getMediaRow(), 1);
		assertEquals(attr.getMediaCol(), 2);
		assertEquals(attr.getPATAssignment(), PATAssignmentMethodType.AffyBarcodeAssignment);
		assertEquals(attr.getMediaFileGUID().toString(), new AffymetrixGuidType("mediaguid").toString());
		assertEquals(attr.getMediaFileName(), "mediafile");
		assertEquals(attr.getLibraryPackageName(), "libpackage");

		assertEquals(array2.getPhysicalArraysAttributes().get(0).getAttributes().size(), 2);
		ParameterNameValue arrayAttr = array2.getPhysicalArraysAttributes().get(0).getAttributes().get(0);
		assertEquals(arrayAttr.getName(), "array-att-name-1");
		assertEquals(arrayAttr.getValueText(), "array-att-value-1");
		arrayAttr = array2.getPhysicalArraysAttributes().get(0).getAttributes().get(1);
		assertEquals(arrayAttr.getName(), "array-att-name-2");
		assertEquals(arrayAttr.getValueText(), "array-att-value-2");

		// User attributes
		assertEquals(array2.getUserAttributes().size(), 3);
		Iterator<ParameterNameValueDefaultRequired> userAttIt = array2.getUserAttributes().iterator();
		ParameterNameValueDefaultRequired userAttr = userAttIt.next();
		assertEquals(userAttr.getName(), "user-att-name-1");
		assertEquals(userAttr.getValueText(), "user-att-value-1");
		assertFalse(userAttr.getRequiredFlag());
		assertEquals(userAttr.getValueType(), ParameterValueType.TextParameterType);
		assertEquals(userAttr.getControlledVocabulary().size(), 0);

		userAttr = userAttIt.next();
		assertEquals(userAttr.getName(), "user-att-name-2");
		assertEquals(userAttr.getValueFloat(), 1, 0.0001);
		assertEquals(userAttr.getDefaultValueFloat(), 2, 0.0001);
		assertTrue(userAttr.getRequiredFlag());
		assertEquals(userAttr.getValueType(), ParameterValueType.FloatParameterType);
		assertEquals(userAttr.getControlledVocabulary().size(), 0);

		userAttr = userAttIt.next();
		assertEquals(userAttr.getName(), "user-att-name-3");
		assertFalse(userAttr.getRequiredFlag());
		assertEquals(userAttr.getValueType(), ParameterValueType.ControlMultiParameterType);

		assertEquals(userAttr.getControlMultiValues().size(), 2);
		Iterator<String> cIt = userAttr.getControlMultiValues().iterator();
		assertEquals(cIt.next(), "one");
		assertEquals(cIt.next(), "two");

		assertEquals(userAttr.getControlledVocabulary().size(), 3);
		cIt = userAttr.getControlledVocabulary().iterator();
		assertEquals(cIt.next(), "one");
		assertEquals(cIt.next(), "two");
		assertEquals(cIt.next(), "three");
	}
}
