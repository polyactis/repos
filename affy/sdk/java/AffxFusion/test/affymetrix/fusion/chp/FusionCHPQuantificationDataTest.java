/*
 * FusionCHPQuantificationDataTest.java
 * JUnit based test
 *
 * Created on December 27, 2005, 8:30 AM
 */

package affymetrix.fusion.chp;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.data.ProbeSetQuantificationData;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.utils.IOUtils;

/**
 * 
 * @author ljevon
 */
public class FusionCHPQuantificationDataTest extends TestCase {

	public FusionCHPQuantificationDataTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() throws Exception {
		TestSuite suite = new TestSuite(FusionCHPQuantificationDataTest.class);

		FusionCHPQuantificationData.registerReader();
		FusionCHPLegacyData.registerReader();

		return suite;
	}

	public void testFileId() throws Exception {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_quantification_file");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		assertTrue(chp != null);
		FusionCHPQuantificationData sigChp = FusionCHPQuantificationData.fromBase(chp);
		assertTrue(sigChp != null);
		assertTrue(Arrays.equals(sigChp.getFileId().getGuid(), "0000065535-1152133546-0000000153-0000003902-0000014604"
				.getBytes(IOUtils.ASCII_CHARSET)));
		assertTrue(Arrays.equals(sigChp.getGenericData().getFileIdentifier().getGuid(),
				"0000065535-1152133546-0000000153-0000003902-0000014604".getBytes(IOUtils.ASCII_CHARSET)));
	}

	public void testReadNonQuantification() throws Exception {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_expression_file");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		assertTrue(chp != null);
		FusionCHPQuantificationData sigChp = FusionCHPQuantificationData.fromBase(chp);
		assertTrue(sigChp == null);
	}

	public void testRead() throws IOException, UnsignedOutOfLimitsException {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_quantification_file");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		assertTrue(chp != null);
		FusionCHPQuantificationData sigChp = FusionCHPQuantificationData.fromBase(chp);
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

		ProbeSetQuantificationData e;
		e = sigChp.getQuantificationEntry(0);
		assertEquals(e.getName(), "abc");
		assertEquals(e.getQuantification(), 10.0f, 0.0001f);
		e = sigChp.getQuantificationEntry(1);
		assertEquals(e.getName(), "xyz");
		assertEquals(e.getQuantification(), 20.0f, 0.0001f);

	}
}
