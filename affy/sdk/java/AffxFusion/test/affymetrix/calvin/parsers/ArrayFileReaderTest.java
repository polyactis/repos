package affymetrix.calvin.parsers;

import java.io.File;
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
import affymetrix.calvin.utils.AffymetrixGuidType;

public class ArrayFileReaderTest extends TestCase {

	/** The identifier for an array file. */
	private static final String ARRAY_SET_FILE_TYPE_IDENTIFIER = "affymetrix-calvin-arraysetfile";

	private String testNonExistantFile = null;// = "../data/file_does_not_exist";

	private String testDataFileInvalid = null;// = "../data/test.file.data_header_only";

	private String testDataFile = null;// = "../data/ArraySetFile.xml";

	protected void setUp() throws Exception {
		super.setUp();
		testNonExistantFile = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\file_does_not_exist")
				.getCanonicalPath();
		testDataFileInvalid = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\test.file.data_header_only")
				.getCanonicalPath();
		testDataFile = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\ArraySetFile.xml").getCanonicalPath();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	private void checkArrayData(boolean headerOnly) {
		ArrayFileReader reader = new ArrayFileReader();
		ArrayData array = new ArrayData();
		String name = testDataFile;
		assertTrue(reader.read(name, array, headerOnly));
		assertEquals(array.getDataTypeIdentifier(), ArrayId.ARRAY_SET_FILE_TYPE_IDENTIFIER);
		assertEquals(array.getArraySetFileIdentifier(), "432-432-432-432");
		assertEquals(array.getCreatedStep(), CreateStepType.ArrayRegistrationStep);
		assertEquals(array.getInitialProject(), "my_project");
		assertEquals(array.getCreationDateTime(), "8/12/2005 9:00AM");
		assertEquals(array.getCreatedBy(), "ljevon");
		if (headerOnly) {
			return;
		}
		assertEquals(array.getPhysicalArraysAttributes().size(), 1);
		ArrayAttributes atts = array.getPhysicalArraysAttributes().get(0);
		assertEquals(atts.getIdentifier(), "123-123-123-123");
		assertEquals(atts.getArrayName(), "mychip");
		assertEquals(atts.getArrayBarcode(), "@1234567890");
		assertEquals(atts.getMedia(), ArrayMediaType.PlateOrStripMedia);
		assertEquals(atts.getMediaRow(), 1);
		assertEquals(atts.getMediaCol(), 12);
		assertEquals(atts.getMasterFile(), "Test3.master");
		assertEquals(atts.getMasterFileId().toString(), new AffymetrixGuidType("123-123-123-123").toString());
		assertEquals(atts.getPATAssignment(), PATAssignmentMethodType.AffyBarcodeAssignment);
		assertEquals(atts.getCreationDateTime(), "8/12/2005 10:00AM");
		assertEquals(atts.getCreatedBy(), "ljevon");
		assertEquals(atts.getCreatedStep(), CreateStepType.ScanningStep);
		assertEquals(atts.getComment(), "here is a comment");
		assertEquals(atts.getLibraryPackageName(), "libpackage");
		assertEquals(atts.getMediaFileGUID(), "mediaguid");
		assertEquals(atts.getMediaFileName(), "mediafile");

		assertEquals(atts.getAttributes().size(), 2);
		Iterator<ParameterNameValue> paramIt = atts.getAttributes().iterator();
		ParameterNameValue p = paramIt.next();
		assertEquals(p.getName(), "SampleID");
		assertEquals(p.getValueText(), "433232");
		p = paramIt.next();
		assertEquals(p.getName(), "LIMS System");
		assertEquals(p.getValueText(), "Nautilis");

		assertEquals(array.getUserAttributes().size(), 6);
		Iterator<ParameterNameValueDefaultRequired> userIt = array.getUserAttributes().iterator();
		ParameterNameValueDefaultRequired pr = userIt.next();
		assertEquals(pr.getName(), "Species");
		assertEquals(pr.getValueText(), "Homo Sapien");
		assertTrue(pr.getRequiredFlag());
		assertEquals(pr.getControlledVocabulary().size(), 0);
		assertEquals(pr.getValueType(), ParameterValueType.TextParameterType);
		pr = userIt.next();
		assertEquals(pr.getName(), "Individual");
		assertEquals(pr.getValueText(), "me");
		assertFalse(pr.getRequiredFlag());
		assertEquals(pr.getControlledVocabulary().size(), 0);
		assertEquals(pr.getValueType(), ParameterValueType.TextParameterType);
		pr = userIt.next();
		assertEquals(pr.getName(), "When");
		assertEquals(pr.getValueText(), "now");
		assertFalse(pr.getRequiredFlag());
		assertEquals(pr.getControlledVocabulary().size(), 0);
		assertEquals(pr.getValueType(), ParameterValueType.TimeParameterType);
		pr = userIt.next();
		assertEquals(pr.getName(), "How much");
		double eps = 1e-5;
		assertEquals(pr.getValueFloat(), 1.0, eps);
		assertFalse(pr.getRequiredFlag());
		assertEquals(pr.getControlledVocabulary().size(), 0);
		assertEquals(pr.getValueType(), ParameterValueType.FloatParameterType);
		pr = userIt.next();
		assertEquals(pr.getName(), "Sex");
		assertEquals(pr.getValueText(), "");
		assertFalse(pr.getRequiredFlag());
		assertEquals(pr.getValueType(), ParameterValueType.ControlMultiParameterType);

		Iterator<String> cit = pr.getControlMultiValues().iterator();
		String s = cit.next();
		assertEquals(s, "Male");
		s = cit.next();
		assertEquals(s, "Other");

		assertEquals(pr.getControlledVocabulary().size(), 4);
		cit = pr.getControlledVocabulary().iterator();
		s = cit.next();
		assertEquals(s, "Male");
		s = cit.next();
		assertEquals(s, "Female");
		s = cit.next();
		assertEquals(s, "Unknown");
		s = cit.next();
		assertEquals(s, "Other");
		pr = userIt.next();
		assertEquals(pr.getName(), "Age");
		assertEquals(pr.getValueText(), "30 something");
		assertFalse(pr.getRequiredFlag());
		assertEquals(pr.getValueType(), ParameterValueType.TextParameterType);
		assertEquals(pr.getControlledVocabulary().size(), 0);
	}

	public void readTest() {
		checkArrayData(false);
	}

	public void readHeaderDataOnlyTest() {
		checkArrayData(true);
	}

	public void readWhenFileDoesNotExistTest() {
		ArrayFileReader reader = new ArrayFileReader();
		ArrayData array = new ArrayData();
		String name = testNonExistantFile;
		assertFalse(reader.read(name, array));
	}

	public void readWhenFileIsNotValidTest() {
		ArrayFileReader reader = new ArrayFileReader();
		ArrayData array = new ArrayData();
		String name = testDataFileInvalid;
		assertFalse(reader.read(name, array));
	}

	public void isFileTypeTest() {
		assertFalse(ArrayFileReader.isFileType(testNonExistantFile, new AffymetrixGuidType(ARRAY_SET_FILE_TYPE_IDENTIFIER)));
		assertFalse(ArrayFileReader.isFileType(testDataFileInvalid, new AffymetrixGuidType(ARRAY_SET_FILE_TYPE_IDENTIFIER)));
		assertTrue(ArrayFileReader.isFileType(testDataFile, new AffymetrixGuidType(ARRAY_SET_FILE_TYPE_IDENTIFIER)));
	}

	public void dataTypeIdentifierTest() {
		String name = testDataFile;
		AffymetrixGuidType guid = ArrayFileReader.getDataTypeIdentifier(name);
		assertEquals(guid.toString(), ARRAY_SET_FILE_TYPE_IDENTIFIER);
	}
}
