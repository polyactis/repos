/*
 * FusionCHPQuantificationDetectionDataTest.java
 * JUnit based test
 *
 * Created on December 27, 2005, 8:30 AM
 */

package affymetrix.fusion.chp;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.data.ProbeSetQuantificationDetectionData;
import affymetrix.calvin.parameter.ParameterNameValue;

/**
 * 
 * @author ljevon
 */
public class FusionCHPQuantificationDetectionDataTest extends TestCase {

	public FusionCHPQuantificationDetectionDataTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() throws Exception {
		TestSuite suite = new TestSuite(FusionCHPQuantificationDetectionDataTest.class);

		FusionCHPQuantificationDetectionData.registerReader();
		FusionCHPLegacyData.registerReader();

		return suite;
	}

	public void testFileId() throws Exception {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_quantification_detection_file");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		assertEquals(chp != null, true);
		FusionCHPQuantificationDetectionData sigChp = FusionCHPQuantificationDetectionData.fromBase(chp);
		assertTrue(sigChp != null);
		assertTrue(Arrays.equals(sigChp.getFileId().getGuid(), "0000065535-1152158444-0000019912-0000001869-0000011538"
				.getBytes("US-ASCII")));
		assertTrue(Arrays.equals(sigChp.getGenericData().getFileIdentifier().getGuid(),
				"0000065535-1152158444-0000019912-0000001869-0000011538".getBytes("US-ASCII")));
	}

	public void testReadNonQuantification() throws Exception {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_expression_file");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		assertTrue(chp != null);
		FusionCHPQuantificationDetectionData sigChp = FusionCHPQuantificationDetectionData.fromBase(chp);
		assertTrue(sigChp == null);
	}

	public void testRead() throws Exception {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_quantification_detection_file");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		assertTrue(chp != null);
		FusionCHPQuantificationDetectionData sigChp = FusionCHPQuantificationDetectionData.fromBase(chp);
		assertTrue(sigChp != null);

		assertEquals(sigChp.getAlgName(), "sig");
		assertEquals(sigChp.getAlgVersion(), "1.0");
		assertEquals(sigChp.getArrayType(), "test3");
		assertEquals(sigChp.getEntryCount(), 2);

		List<ParameterNameValue> params = sigChp.getAlgParams();
		assertEquals(params.size(), 1);
		ParameterNameValue param = params.get(0);
		assertEquals(param.getName(), "an1");
		assertEquals(param.getValueText(), "av1");

		params = sigChp.getSummaryParams();
		assertEquals(params.size(), 1);
		param = params.get(0);
		assertEquals(param.getName(), "sn1");
		assertEquals(param.getValueText(), "sv1");

		ProbeSetQuantificationDetectionData e;
		e = sigChp.getQuantificationDetectionEntry(0);
		assertEquals(e.getName(), "abc");
		assertEquals(e.getQuantification(), 10.0f, 0.0001f);
		assertEquals(e.getPValue(), 0.1f, 0.0001f);
		e = sigChp.getQuantificationDetectionEntry(1);
		assertEquals(e.getName(), "xyz");
		assertEquals(e.getQuantification(), 20.0f, 0.0001f);
		assertEquals(e.getPValue(), 0.2f, 0.0001f);

	}

}
