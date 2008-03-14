/*
 * FusionCHPTilingDataTest.java
 * JUnit based test
 *
 * Created on November 30, 2005, 9:07 PM
 */

package affymetrix.fusion.chp;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import junit.framework.TestCase;
import affymetrix.calvin.data.CHPTilingEntry;
import affymetrix.calvin.data.TilingSequenceData;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.utils.AffymetrixGuidType;

/**
 * 
 * @author ljevon
 */
public class FusionCHPTilingDataTest extends TestCase {

	public FusionCHPTilingDataTest(String testName) throws Exception {
		super(testName);
		FusionCHPLegacyData.registerReader();
		FusionCHPGenericData.registerReader();
		FusionCHPTilingData.registerReader();
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public void testFileId() throws Exception {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_tiling_file");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		assertTrue(chp != null);
		FusionCHPTilingData tileChp = FusionCHPTilingData.fromBase(chp);
		assertTrue(tileChp != null);
		AffymetrixGuidType guid = new AffymetrixGuidType("0000039321-1131034775-0000032391-0000005436-0000004827");
		assertTrue(Arrays.equals(chp.getFileId().getGuid(), guid.getGuid()));
	}

	public void testReadNonTiling() throws Exception {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_expression_file");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		assertTrue(chp != null);
		FusionCHPTilingData tileChp = FusionCHPTilingData.fromBase(chp);
		assertTrue(tileChp == null);
	}

	public void testRead() throws Exception {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_tiling_file");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		assertTrue(chp != null);
		FusionCHPTilingData tileChp = FusionCHPTilingData.fromBase(chp);
		assertTrue(tileChp != null);

		assertTrue(tileChp.getNumberSequences() == 2);
		assertTrue(tileChp.getAlgName().compareTo("tile") == 0);
		assertTrue(tileChp.getAlgVersion().compareTo("1.0") == 0);

		List<ParameterNameValue> params = tileChp.getAlgParams();
		assertTrue(params.size() == 1);
		ParameterNameValue param = params.get(0);
		assertTrue(param.getName().compareTo("p1") == 0);
		assertTrue(param.getValueText().compareTo("v1") == 0);

		double eps = 0.00001;
		CHPTilingEntry e;
		TilingSequenceData seq;

		tileChp.openTilingSequenceDataSet(0);
		seq = tileChp.getTilingSequenceData();

		assertTrue(seq.getName().compareTo("n1") == 0);
		assertTrue(seq.getGroupName().compareTo("g1") == 0);
		assertTrue(seq.getVersion().compareTo("v1") == 0);

		params = seq.getParameters();
		assertTrue(params.size() == 1);
		param = params.get(0);
		assertTrue(param.getName().compareTo("seq1_p1") == 0);
		assertTrue(param.getValueText().compareTo("seq1_v1") == 0);

		assertTrue(tileChp.getTilingSequenceEntryCount(0) == 2);
		e = tileChp.getTilingSequenceEntry(0);
		assertTrue(e.getPosition() == 10);
		assertEquals(e.getValue(), 10.0f, eps);

		e = tileChp.getTilingSequenceEntry(1);
		assertTrue(e.getPosition() == 20);
		assertEquals(e.getValue(), 20.0f, eps);

		tileChp.openTilingSequenceDataSet(1);
		seq = tileChp.getTilingSequenceData();

		assertTrue(seq.getName().compareTo("n2") == 0);
		assertTrue(seq.getGroupName().compareTo("g2") == 0);
		assertTrue(seq.getVersion().compareTo("v2") == 0);
		params = seq.getParameters();
		assertTrue(params.size() == 1);
		param = params.get(0);
		assertTrue(param.getName().compareTo("seq2_p1") == 0);
		assertTrue(param.getValueText().compareTo("seq2_v1") == 0);

		assertTrue(tileChp.getTilingSequenceEntryCount(1) == 3);
		e = tileChp.getTilingSequenceEntry(0);
		assertTrue(e.getPosition() == 11);
		assertEquals(e.getValue(), 11.0f, eps);

		e = tileChp.getTilingSequenceEntry(1);
		assertTrue(e.getPosition() == 21);
		assertEquals(e.getValue(), 21.0f, eps);

		e = tileChp.getTilingSequenceEntry(2);
		assertTrue(e.getPosition() == 31);
		assertEquals(e.getValue(), 31.0f, eps);
	}
}
