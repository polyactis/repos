/*
 * CHPFileDataTest.java
 * JUnit based test
 *
 * Created on October 13, 2005, 11:00 PM
 */

package affymetrix.gcos.chp;

import java.io.File;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.gcos.TagValuePair;

/**
 * 
 * @author ljevon
 */
public class CHPFileDataTest extends TestCase {

	public CHPFileDataTest(String testName) throws Exception {
		super(testName);
		EXP_COMP_CHP_FILE = new File("..\\..\\..\\CPPTest\\data\\test.exp.comp.chp").getCanonicalPath();
		EXP_ABS_CHP_FILE = new File("..\\..\\..\\CPPTest\\data\\test.exp.abs.chp").getCanonicalPath();
		NO_FILE = new File("..\\..\\..\\CPPTest\\data\\no_file.chp").getCanonicalPath();
		TAG_V11_FILE = new File("..\\..\\..\\CPPTest\\data\\tag_v11.chp").getCanonicalPath();
		TAG_XDA_FILE = new File("..\\..\\..\\CPPTest\\data\\tag_XDA.chp").getCanonicalPath();
		RESEQ_V1_FILE = new File("..\\..\\..\\CPPTest\\data\\reseq_xda_v1.chp").getCanonicalPath();
		RESEQ_V2_FILE = new File("..\\..\\..\\CPPTest\\data\\reseq_xda_v2.chp").getCanonicalPath();
		RESEQ_OLD_FILE = new File("..\\..\\..\\CPPTest\\data\\reseq_v13.chp").getCanonicalPath();
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(CHPFileDataTest.class);

		return suite;
	}

	private String EXP_COMP_CHP_FILE;

	private String EXP_ABS_CHP_FILE;

	private String NO_FILE;

	private String TAG_V11_FILE;

	private String TAG_XDA_FILE;

	private String RESEQ_V1_FILE;

	private String RESEQ_V2_FILE;

	private String RESEQ_OLD_FILE;

	public void testPropertyPaths() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(EXP_COMP_CHP_FILE);
		assertEquals(chp.getFileName(), EXP_COMP_CHP_FILE);
	}

	public void testReadHeaderFail() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(NO_FILE);
		assertFalse(chp.readHeader());
	}

	public void testReadFail() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(NO_FILE);
		assertFalse(chp.read());
	}

	public void testExists() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(NO_FILE);
		assertFalse(chp.exists());
		chp.setFileName(EXP_COMP_CHP_FILE);
		assertTrue(chp.exists());
		chp.setFileName(EXP_ABS_CHP_FILE);
		assertTrue(chp.exists());
	}

	public void testIsXDACompatibleFile() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(EXP_COMP_CHP_FILE);
		assertTrue(chp.isXDACompatibleFile());
		chp.setFileName(EXP_ABS_CHP_FILE);
		assertTrue(chp.isXDACompatibleFile());
	}

	public void testReadExpComp() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(EXP_COMP_CHP_FILE);
		assertTrue(chp.read());
		assertEquals(chp.getFileName(), EXP_COMP_CHP_FILE);

		// Accessors for header information.
		// assertEquals( chp.getHeader().getVersionNumber(), 1 );
		assertEquals(chp.getHeader().getCols(), 120);
		assertEquals(chp.getHeader().getRows(), 120);
		assertEquals(chp.getHeader().getNumProbeSets(), 8);
		assertEquals(chp.getHeader().getChipType(), "TestExon");
		assertEquals(chp.getHeader().getAlgName(), "Plier");
		assertEquals(chp.getHeader().getAlgorithmParameter("p1"), "1.1.1");
		assertEquals(chp.getHeader().getAlgorithmParameter("p2"), "2.1.1");
		assertEquals(chp.getHeader().getSummaryParameter("cp1"), "1.2.1");
		assertEquals(chp.getHeader().getSummaryParameter("cp2"), "2.2.1");
		assertEquals(chp.getHeader().getAlgVersion(), "1.0");
		assertEquals(chp.getHeader().getParentCellFile(), "TestExon.cel");
		assertEquals(chp.getHeader().getProgID(), "test_id");

		assertEquals(chp.getHeader().getBackgroundZoneInfo().getNumberZones(), 5);
		assertEquals(chp.getHeader().getBackgroundZoneInfo().getSmoothFactor(), 3.8f, 0.00001f);
		assertEquals(chp.getHeader().getBackgroundZone(0, 0).getBackground(), 3.2f, 0.00001f);
		assertEquals(chp.getHeader().getBackgroundZone(88, 33).getBackground(), 8.9f, 0.00001f);
		assertEquals(chp.getHeader().getBackgroundZone(0, 60).getBackground(), 13.1f, 0.00001f);
		assertEquals(chp.getHeader().getBackgroundZone(60, 0).getBackground(), 5.0f, 0.00001f);
		assertEquals(chp.getHeader().getBackgroundZone(120, 120).getBackground(), 10.8f, 0.00001f);
		int n = chp.getHeader().getBackgroundZoneInfo().getNumberZones();
		for (int i = 0; i < n; i++) {
			BackgroundZoneType zone = chp.getHeader().getBackgroundZoneInfo().getZone(i);
			switch (i) {
			case 0:
				assertEquals(zone.getCenterX(), 0.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 0.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 3.2f, 0.0000001f);
				break;
			case 1:
				assertEquals(zone.getCenterX(), 88.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 33.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 8.9f, 0.0000001f);
				break;
			case 2:
				assertEquals(zone.getCenterX(), 0.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 60.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 13.1f, 0.0000001f);
				break;
			case 3:
				assertEquals(zone.getCenterX(), 60.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 0.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 5.0f, 0.0000001f);
				break;
			case 4:
				assertEquals(zone.getCenterX(), 120.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 120.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 10.8f, 0.0000001f);
				break;
			}
		}

		for (int i = 0; i < chp.getHeader().getNumProbeSets(); i++) {
			assertEquals(chp.getExpressionResults(i).getDetectionPValue(), (0.05f - (i / 1000.0)), 0.00001f);
			assertEquals(chp.getExpressionResults(i).getSignal(), (1.1f + i), 0.00001f);
			assertEquals(chp.getExpressionResults(i).getNumPairs().toInt(), 3 + i);
			assertEquals(chp.getExpressionResults(i).getNumUsedPairs().toInt(), 2 + i);
			assertEquals(chp.getExpressionResults(i).getDetection().toShort(), (i % 4));
			switch (chp.getExpressionResults(i).getDetection().toShort()) {
			case (ExpressionProbeSetResults.ABS_PRESENT_CALL):
				assertEquals(chp.getExpressionResults(i).getDetectionString(), "P");
				break;
			case (ExpressionProbeSetResults.ABS_MARGINAL_CALL):
				assertEquals(chp.getExpressionResults(i).getDetectionString(), "M");
				break;
			case (ExpressionProbeSetResults.ABS_ABSENT_CALL):
				assertEquals(chp.getExpressionResults(i).getDetectionString(), "A");
				break;
			case (ExpressionProbeSetResults.ABS_NO_CALL):
				assertEquals(chp.getExpressionResults(i).getDetectionString(), "No Call");
				break;
			default:
				assertEquals(chp.getExpressionResults(i).getDetectionString(), "");
				break;
			}
			assertTrue(chp.getExpressionResults(i).getHasCompResults());
			assertEquals(chp.getExpressionResults(i).getChangePValue(), (0.04f - (i / 1000.0)), 0.0000001f);
			assertEquals(chp.getExpressionResults(i).getSignalLogRatio(), (1.1f + i), 0.0000001f);
			assertEquals(chp.getExpressionResults(i).getSignalLogRatioLow(), (-1.1f + i), 0.0000001f);
			assertEquals(chp.getExpressionResults(i).getSignalLogRatioHigh(), (10.1f + i), 0.0000001f);
			assertEquals(chp.getExpressionResults(i).getNumCommonPairs().toInt(), 2 + i);
			assertEquals(chp.getExpressionResults(i).getChange().toShort(), (i % 6 + 1));
			switch (chp.getExpressionResults(i).getChange().toShort()) {
			case (ExpressionProbeSetResults.COMP_INCREASE_CALL):
				assertEquals(chp.getExpressionResults(i).getChangeString(), "I");
				break;
			case (ExpressionProbeSetResults.COMP_DECREASE_CALL):
				assertEquals(chp.getExpressionResults(i).getChangeString(), "D");
				break;
			case (ExpressionProbeSetResults.COMP_MOD_INCREASE_CALL):
				assertEquals(chp.getExpressionResults(i).getChangeString(), "MI");
				break;
			case (ExpressionProbeSetResults.COMP_MOD_DECREASE_CALL):
				assertEquals(chp.getExpressionResults(i).getChangeString(), "MD");
				break;
			case (ExpressionProbeSetResults.COMP_NO_CHANGE_CALL):
				assertEquals(chp.getExpressionResults(i).getChangeString(), "NC");
				break;
			case (ExpressionProbeSetResults.COMP_NO_CALL):
				assertEquals(chp.getExpressionResults(i).getChangeString(), "No Call");
				break;
			default:
				assertEquals(chp.getExpressionResults(i).getChangeString(), "");
				break;
			}
		}
	}

	public void testmethod_Read_Exp_Abs() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(EXP_ABS_CHP_FILE);
		assertTrue(chp.read());
		assertEquals(chp.getFileName(), EXP_ABS_CHP_FILE);

		// Accessors for header information.
		// assertEquals( chp.getHeader().getVersionNumber(), 1 );
		assertEquals(chp.getHeader().getCols(), 120);
		assertEquals(chp.getHeader().getRows(), 120);
		assertEquals(chp.getHeader().getNumProbeSets(), 8);
		assertEquals(chp.getHeader().getChipType(), "TestExon");
		assertEquals(chp.getHeader().getAlgName(), "Plier");
		assertEquals(chp.getHeader().getAlgorithmParameter("p1"), "1.1.1");
		assertEquals(chp.getHeader().getAlgorithmParameter("p2"), "2.1.1");
		assertEquals(chp.getHeader().getSummaryParameter("cp1"), "1.2.1");
		assertEquals(chp.getHeader().getSummaryParameter("cp2"), "2.2.1");
		assertEquals(chp.getHeader().getAlgVersion(), "1.0");
		assertEquals(chp.getHeader().getParentCellFile(), "TestExon.cel");
		assertEquals(chp.getHeader().getProgID(), "test_id");

		assertEquals(chp.getHeader().getBackgroundZoneInfo().getNumberZones(), 5);
		assertEquals(chp.getHeader().getBackgroundZoneInfo().getSmoothFactor(), 3.8f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(0, 0).getBackground(), 3.2f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(88, 33).getBackground(), 8.9f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(0, 60).getBackground(), 13.1f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(60, 0).getBackground(), 5.0f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(120, 120).getBackground(), 10.8f, 0.000001f);
		int n = chp.getHeader().getBackgroundZoneInfo().getNumberZones();
		for (int i = 0; i < n; i++) {
			BackgroundZoneType zone = chp.getHeader().getBackgroundZoneInfo().getZone(i);
			switch (i) {
			case 0:
				assertEquals(zone.getCenterX(), 0.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 0.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 3.2f, 0.0000001f);
				break;
			case 1:
				assertEquals(zone.getCenterX(), 88.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 33.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 8.9f, 0.0000001f);
				break;
			case 2:
				assertEquals(zone.getCenterX(), 0.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 60.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 13.1f, 0.0000001f);
				break;
			case 3:
				assertEquals(zone.getCenterX(), 60.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 0.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 5.0f, 0.0000001f);
				break;
			case 4:
				assertEquals(zone.getCenterX(), 120.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 120.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 10.8f, 0.0000001f);
				break;
			}
		}

		for (int i = 0; i < chp.getHeader().getNumProbeSets(); i++) {
			assertEquals(chp.getExpressionResults(i).getDetectionPValue(), (float)(0.05 - (i / 1000.0)), 0.000001f);
			assertEquals(chp.getExpressionResults(i).getSignal(), (float)(1.1 + i), 0.000001f);
			assertEquals(chp.getExpressionResults(i).getNumPairs().toInt(), 3 + i);
			assertEquals(chp.getExpressionResults(i).getNumUsedPairs().toInt(), 2 + i);
			assertEquals(chp.getExpressionResults(i).getDetection().toShort(), (i % 4));
			switch (chp.getExpressionResults(i).getDetection().toShort()) {
			case (ExpressionProbeSetResults.ABS_PRESENT_CALL):
				assertEquals(chp.getExpressionResults(i).getDetectionString(), "P");
				break;
			case (ExpressionProbeSetResults.ABS_MARGINAL_CALL):
				assertEquals(chp.getExpressionResults(i).getDetectionString(), "M");
				break;
			case (ExpressionProbeSetResults.ABS_ABSENT_CALL):
				assertEquals(chp.getExpressionResults(i).getDetectionString(), "A");
				break;
			case (ExpressionProbeSetResults.ABS_NO_CALL):
				assertEquals(chp.getExpressionResults(i).getDetectionString(), "No Call");
				break;
			default:
				assertEquals(chp.getExpressionResults(i).getDetectionString(), "");
				break;
			}
			assertFalse(chp.getExpressionResults(i).getHasCompResults());
		}
	}

	public void testmethod_ReadHeader_Exp_Comp() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(EXP_COMP_CHP_FILE);
		assertTrue(chp.readHeader());
		assertEquals(chp.getFileName(), EXP_COMP_CHP_FILE);

		// Accessors for header information.
		// assertEquals( chp.getHeader().getVersionNumber(), 1 );
		assertEquals(chp.getHeader().getCols(), 120);
		assertEquals(chp.getHeader().getRows(), 120);
		assertEquals(chp.getHeader().getNumProbeSets(), 8);
		assertEquals(chp.getHeader().getChipType(), "TestExon");
		assertEquals(chp.getHeader().getAlgName(), "Plier");
		assertEquals(chp.getHeader().getAlgorithmParameter("p1"), "1.1.1");
		assertEquals(chp.getHeader().getAlgorithmParameter("p2"), "2.1.1");
		assertEquals(chp.getHeader().getSummaryParameter("cp1"), "1.2.1");
		assertEquals(chp.getHeader().getSummaryParameter("cp2"), "2.2.1");
		assertEquals(chp.getHeader().getAlgVersion(), "1.0");
		assertEquals(chp.getHeader().getParentCellFile(), "TestExon.cel");
		assertEquals(chp.getHeader().getProgID(), "test_id");

		assertEquals(chp.getHeader().getBackgroundZoneInfo().getNumberZones(), 5);
		assertEquals(chp.getHeader().getBackgroundZoneInfo().getSmoothFactor(), 3.8f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(0, 0).getBackground(), 3.2f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(88, 33).getBackground(), 8.9f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(0, 60).getBackground(), 13.1f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(60, 0).getBackground(), 5.0f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(120, 120).getBackground(), 10.8f, 0.000001f);

		int n = chp.getHeader().getBackgroundZoneInfo().getNumberZones();
		for (int i = 0; i < n; i++) {
			BackgroundZoneType zone = chp.getHeader().getBackgroundZoneInfo().getZone(i);
			switch (i) {
			case 0:
				assertEquals(zone.getCenterX(), 0.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 0.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 3.2f, 0.0000001f);
				break;
			case 1:
				assertEquals(zone.getCenterX(), 88.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 33.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 8.9f, 0.0000001f);
				break;
			case 2:
				assertEquals(zone.getCenterX(), 0.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 60.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 13.1f, 0.0000001f);
				break;
			case 3:
				assertEquals(zone.getCenterX(), 60.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 0.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 5.0f, 0.0000001f);
				break;
			case 4:
				assertEquals(zone.getCenterX(), 120.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 120.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 10.8f, 0.0000001f);
				break;
			}
		}
	}

	public void testmethod_ReadHeader_Exp_Abs() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(EXP_ABS_CHP_FILE);
		assertTrue(chp.readHeader());
		assertEquals(chp.getFileName(), EXP_ABS_CHP_FILE);

		// Accessors for header information.
		// assertEquals( chp.getHeader().getVersionNumber(), 1 );
		assertEquals(chp.getHeader().getCols(), 120);
		assertEquals(chp.getHeader().getRows(), 120);
		assertEquals(chp.getHeader().getNumProbeSets(), 8);
		assertEquals(chp.getHeader().getChipType(), "TestExon");
		assertEquals(chp.getHeader().getAlgName(), "Plier");
		assertEquals(chp.getHeader().getAlgorithmParameter("p1"), "1.1.1");
		assertEquals(chp.getHeader().getAlgorithmParameter("p2"), "2.1.1");
		assertEquals(chp.getHeader().getSummaryParameter("cp1"), "1.2.1");
		assertEquals(chp.getHeader().getSummaryParameter("cp2"), "2.2.1");
		assertEquals(chp.getHeader().getAlgVersion(), "1.0");
		assertEquals(chp.getHeader().getParentCellFile(), "TestExon.cel");
		assertEquals(chp.getHeader().getProgID(), "test_id");

		assertEquals(chp.getHeader().getBackgroundZoneInfo().getNumberZones(), 5);
		assertEquals(chp.getHeader().getBackgroundZoneInfo().getSmoothFactor(), 3.8f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(0, 0).getBackground(), 3.2f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(88, 33).getBackground(), 8.9f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(0, 60).getBackground(), 13.1f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(60, 0).getBackground(), 5.0f, 0.000001f);
		assertEquals(chp.getHeader().getBackgroundZone(120, 120).getBackground(), 10.8f, 0.000001f);

		int n = chp.getHeader().getBackgroundZoneInfo().getNumberZones();
		for (int i = 0; i < n; i++) {
			BackgroundZoneType zone = chp.getHeader().getBackgroundZoneInfo().getZone(i);
			switch (i) {
			case 0:
				assertEquals(zone.getCenterX(), 0.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 0.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 3.2f, 0.0000001f);
				break;
			case 1:
				assertEquals(zone.getCenterX(), 88.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 33.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 8.9f, 0.0000001f);
				break;
			case 2:
				assertEquals(zone.getCenterX(), 0.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 60.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 13.1f, 0.0000001f);
				break;
			case 3:
				assertEquals(zone.getCenterX(), 60.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 0.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 5.0f, 0.0000001f);
				break;
			case 4:
				assertEquals(zone.getCenterX(), 120.0f, 0.0000001f);
				assertEquals(zone.getCenterY(), 120.0f, 0.0000001f);
				assertEquals(zone.getBackground(), 10.8f, 0.0000001f);
				break;
			}
		}
	}

	public void testmethod_ReadHeader_TagV11() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(TAG_V11_FILE);
		assertTrue(chp.readHeader());

		// assertEquals( chp.getHeader().getVersionNumber(), 11 );
		assertEquals(chp.getHeader().getCols(), 105);
		assertEquals(chp.getHeader().getRows(), 105);
		assertEquals(chp.getHeader().getNumProbeSets(), 2050);
		assertEquals(chp.getHeader().getChipType(), "GenFlex");
		assertEquals(chp.getHeader().getAlgName(), "Hybridization");
		assertEquals(chp.getHeader().getAlgVersion(), "4.0");
		assertEquals(chp.getHeader().getParentCellFile(), "c:\\genechip\\testdata\\universal\\1109-A-05a.CEL");
		assertEquals(chp.getHeader().getProgID(), "GeneChipAnalysis.HybBaseCall.1");
		assertEquals(chp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(chp.getHeader().getSummaryParameters(), null);
		assertEquals(chp.getHeader().getAlgorithmParameters(), null);
	}

	public void testmethod_Read_TagV11() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(TAG_V11_FILE);
		assertEquals(chp.read(), true);

		// assertEquals( chp.getHeader().getVersionNumber(), 11 );
		assertEquals(chp.getHeader().getCols(), 105);
		assertEquals(chp.getHeader().getRows(), 105);
		assertEquals(chp.getHeader().getNumProbeSets(), 2050);
		assertEquals(chp.getHeader().getChipType(), "GenFlex");
		assertEquals(chp.getHeader().getAlgName(), "Hybridization");
		assertEquals(chp.getHeader().getAlgVersion(), "4.0");
		assertEquals(chp.getHeader().getParentCellFile(), "c:\\genechip\\testdata\\universal\\1109-A-05a.CEL");
		assertEquals(chp.getHeader().getProgID(), "GeneChipAnalysis.HybBaseCall.1");
		assertEquals(chp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(chp.getHeader().getSummaryParameters(), null);
		assertEquals(chp.getHeader().getAlgorithmParameters(), null);

		assertEquals(chp.getUniversalResults(0).getBackground(), 114.595f, 0.000001f);
		assertEquals(chp.getUniversalResults(1).getBackground(), 118.9f, 0.000001f);
		assertEquals(chp.getUniversalResults(2049).getBackground(), 114.854f, 0.000001f);
	}

	public void testmethod_ReadHeader_TagXDA() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(TAG_XDA_FILE);
		assertTrue(chp.readHeader());

		// assertEquals( chp.getHeader().getVersionNumber(), 1 );
		assertEquals(chp.getHeader().getCols(), 105);
		assertEquals(chp.getHeader().getRows(), 105);
		assertEquals(chp.getHeader().getNumProbeSets(), 2050);
		assertEquals(chp.getHeader().getChipType(), "GenFlex");
		assertEquals(chp.getHeader().getAlgName(), "Hybridization");
		assertEquals(chp.getHeader().getAlgVersion(), "5.0");
		assertEquals(chp.getHeader().getParentCellFile(),
				"C:\\Program Files\\Affymetrix\\GeneChip\\Affy_Data\\Data\\1109-A-05a.CEL");
		assertEquals(chp.getHeader().getProgID(), "GDMTAnalysis.HybBaseCall.1");
		assertEquals(chp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(chp.getHeader().getSummaryParameters().size(), 1);
		TagValuePair param = chp.getHeader().getSummaryParameters().get(0);
		assertEquals(param.getTag(), "WaveLength");
		assertEquals(param.getValue(), null);
		assertEquals(chp.getHeader().getAlgorithmParameters(), null);
	}

	public void testReadTagXDA() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(TAG_XDA_FILE);
		assertTrue(chp.read());

		// assertEquals( chp.getHeader().getVersionNumber(), 1 );
		assertEquals(chp.getHeader().getCols(), 105);
		assertEquals(chp.getHeader().getRows(), 105);
		assertEquals(chp.getHeader().getNumProbeSets(), 2050);
		assertEquals(chp.getHeader().getChipType(), "GenFlex");
		assertEquals(chp.getHeader().getAlgName(), "Hybridization");
		assertEquals(chp.getHeader().getAlgVersion(), "5.0");
		assertEquals(chp.getHeader().getParentCellFile(),
				"C:\\Program Files\\Affymetrix\\GeneChip\\Affy_Data\\Data\\1109-A-05a.CEL");
		assertEquals(chp.getHeader().getProgID(), "GDMTAnalysis.HybBaseCall.1");
		assertEquals(chp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(chp.getHeader().getSummaryParameters().size(), 1);
		TagValuePair param = chp.getHeader().getSummaryParameters().get(0);
		assertEquals(param.getTag(), "WaveLength");
		assertTrue(param.getValue() == null);
		assertEquals(chp.getHeader().getAlgorithmParameters(), null);
		assertEquals(chp.getUniversalResults(1).getBackground(), 114.59585f, 0.000001f);
		assertEquals(chp.getUniversalResults(2).getBackground(), 118.9f, 0.000001f);
		assertEquals(chp.getUniversalResults(2049).getBackground(), 114.85418f, 0.000001f);
	}

	public void testReadReseqXDA_v1() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(RESEQ_V1_FILE);
		assertTrue(chp.read());

		// assertEquals( chp.getHeader().getVersionNumber(), 1 );
		assertEquals(chp.getHeader().getCols(), 488);
		assertEquals(chp.getHeader().getRows(), 639);
		assertEquals(chp.getHeader().getNumProbeSets(), 2);
		assertEquals(chp.getHeader().getChipType(), "DCNtagIQr510989");
		assertEquals(chp.getHeader().getAlgName(), "CustomSeq");
		assertEquals(chp.getHeader().getAlgVersion(), "2");
		assertEquals(chp.getHeader().getParentCellFile(), "C:\\Data\\GCOS\\Data\\5303_DCN_01.CEL");
		assertEquals(chp.getHeader().getProgID(), "GDMTAnalysis.VDABaseCall.1");
		assertEquals(chp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(chp.getHeader().getSummaryParameters(), null);
		assertEquals(chp.getHeader().getAlgorithmParameters().size(), 10);
		TagValuePair param = chp.getHeader().getAlgorithmParameters().get(0);
		assertEquals(param.getTag(), "NoSignal");
		assertEquals(param.getValue(), "1.000000");
		param = chp.getHeader().getAlgorithmParameters().get(1);
		assertEquals(param.getTag(), "WeakSignal");
		assertEquals(param.getValue(), "20.000000");
		param = chp.getHeader().getAlgorithmParameters().get(2);
		assertEquals(param.getTag(), "AberrantSNR2");
		assertEquals(param.getValue(), "20.000000");
		param = chp.getHeader().getAlgorithmParameters().get(3);
		assertEquals(param.getTag(), "StrandLLR");
		assertEquals(param.getValue(), "0.000000");
		param = chp.getHeader().getAlgorithmParameters().get(4);
		assertEquals(param.getTag(), "TotalLLR");
		assertEquals(param.getValue(), "75.000000");
		param = chp.getHeader().getAlgorithmParameters().get(5);
		assertEquals(param.getTag(), "PerfectCallThreshold");
		assertEquals(param.getValue(), "2.000000");
		param = chp.getHeader().getAlgorithmParameters().get(6);
		assertEquals(param.getTag(), "ModelType");
		assertEquals(param.getValue(), "0");
		param = chp.getHeader().getAlgorithmParameters().get(7);
		assertEquals(param.getTag(), "FinalMaxHet");
		assertEquals(param.getValue(), "0.900000");
		param = chp.getHeader().getAlgorithmParameters().get(8);
		assertEquals(param.getTag(), "NeighborhoodRule");
		assertEquals(param.getValue(), "0.500000");
		param = chp.getHeader().getAlgorithmParameters().get(9);
		assertEquals(param.getTag(), "SampleReliability");
		assertEquals(param.getValue(), "0.750000");

		ResequencingResults p = chp.getResequencingResults();

		int n = p.getCalledBasesSize();
		assertEquals(n, 27670);
		assertEquals(p.getCalledBase(0), 'c');
		assertEquals(p.getCalledBase(1), 'c');
		assertEquals(p.getCalledBase(2), 'c');
		assertEquals(p.getCalledBase(3), 'a');
		assertEquals(p.getCalledBase(4), 'g');
		assertEquals(p.getCalledBase(5), 'n');
		assertEquals(p.getCalledBase(6), 'c');

		assertEquals(p.getCalledBase(n - 1), 'n');
		assertEquals(p.getCalledBase(n - 44), 'a');

		final float eps = 1e-5f;
		assertEquals(p.getScore(0), 66.465454, eps);
		assertEquals(p.getScore(1), 120.09094, eps);
		assertEquals(p.getScore(2), 153.29333, eps);
		assertEquals(p.getScore(3), 120.06390, eps);
		assertEquals(p.getScore(n - 1), 0.0, eps);
		assertEquals(p.getScore(n - 6), 75.680908, eps);

		assertEquals(p.getForceCallsSize(), 0);
		assertEquals(p.getOrigCallsSize(), 0);
	}

	public void testReadReseqXDA_v2() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(RESEQ_V2_FILE);
		assertEquals(chp.read(), true);

		// assertEquals( chp.getHeader().getVersionNumber(), 2 );
		assertEquals(chp.getHeader().getCols(), 960);
		assertEquals(chp.getHeader().getRows(), 1008);
		assertEquals(chp.getHeader().getNumProbeSets(), 2);
		assertEquals(chp.getHeader().getChipType(), "Tristezar520098");
		assertEquals(chp.getHeader().getAlgName(), "Reseq2");
		assertEquals(chp.getHeader().getAlgVersion(), "2.0");
		assertEquals(chp.getHeader().getParentCellFile(), "C:\\GeneChip\\Affy_Data\\Data\\Hyb01004 CTV-T36 Expt 1926.CEL");
		assertEquals(chp.getHeader().getProgID(), "GDMTAnalysis.Reseq2BaseCall.1");
		assertEquals(chp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(chp.getHeader().getSummaryParameters(), null);
		assertEquals(chp.getHeader().getAlgorithmParameters().size(), 1);
		TagValuePair param = chp.getHeader().getAlgorithmParameters().get(0);
		assertEquals(param.getTag(), "QualityScore");
		assertEquals(param.getValue(), "3");

		ResequencingResults p = chp.getResequencingResults();

		int n = p.getCalledBasesSize();
		assertEquals(n, 111830);
		assertEquals(p.getCalledBase(0), 'a');
		assertEquals(p.getCalledBase(1), 'c');
		assertEquals(p.getCalledBase(2), 'a');
		assertEquals(p.getCalledBase(3), 'g');
		assertEquals(p.getCalledBase(4), 'c');
		assertEquals(p.getCalledBase(5), 'g');
		assertEquals(p.getCalledBase(6), 'a');

		assertEquals(p.getCalledBase(n - 1), 'g');
		assertEquals(p.getCalledBase(n - 2), 'g');
		assertEquals(p.getCalledBase(n - 3), 't');

		final float eps = 1e-5f;
		assertEquals(p.getScore(0), 18.859673, eps);
		assertEquals(p.getScore(1), 14.588547, eps);
		assertEquals(p.getScore(2), 15.489632, eps);
		assertEquals(p.getScore(3), 25.993202, eps);
		assertEquals(p.getScore(n - 1), 32.353516, eps);
		assertEquals(p.getScore(n - 2), 16.441238, eps);

		assertEquals(p.getOrigCallsSize(), 0);

		n = p.getForceCallsSize();
		assertEquals(n, 61065);

		ForceCallType force;

		force = p.getForceCall(0);
		assertEquals(force.getPosition(), 24);
		assertEquals(force.getCall(), 'n');
		assertEquals(force.getReason(), ForceCallType.QUALITY_SCORE_THR_FORCE_CALL);

		force = p.getForceCall(1);
		assertEquals(force.getPosition(), 31);
		assertEquals(force.getCall(), 'n');
		assertEquals(force.getReason(), ForceCallType.QUALITY_SCORE_THR_FORCE_CALL);

		force = p.getForceCall(2);
		assertEquals(force.getPosition(), 39);
		assertEquals(force.getCall(), 'n');
		assertEquals(force.getReason(), ForceCallType.QUALITY_SCORE_THR_FORCE_CALL);

		force = p.getForceCall(n - 1);
		assertEquals(force.getPosition(), 111795);
		assertEquals(force.getCall(), 'n');
		assertEquals(force.getReason(), ForceCallType.QUALITY_SCORE_THR_FORCE_CALL);
	}

	public void testmethod_Read_Reseq_old_file() {
		CHPFileData chp = new CHPFileData();
		chp.setFileName(RESEQ_OLD_FILE);
		assertTrue(chp.read());

		// assertEquals( chp.getHeader().getVersionNumber(), 13 );
		assertEquals(chp.getHeader().getCols(), 488);
		assertEquals(chp.getHeader().getRows(), 639);
		assertEquals(chp.getHeader().getNumProbeSets(), 0);
		assertEquals(chp.getHeader().getChipType(), "DCNtagIQr510989");
		assertEquals(chp.getHeader().getAlgName(), "CustomSeq");
		assertEquals(chp.getHeader().getAlgVersion(), "1");
		assertEquals(chp.getHeader().getParentCellFile(),
				"S:\\GDAS2\\Test_Files\\Resequence\\DCN Validation Data\\5303_DCN_01.cel");
		assertEquals(chp.getHeader().getProgID(), "GDMTAnalysis.VDABaseCall.1");
		assertEquals(chp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(chp.getHeader().getSummaryParameters(), null);
		assertEquals(chp.getHeader().getAlgorithmParameters().size(), 10);
		for (int i = 0; i < chp.getHeader().getAlgorithmParameters().size(); i++) {
			TagValuePair param = chp.getHeader().getAlgorithmParameters().get(i);
			if (param.getTag().compareTo("SampleReliability") == 0)
				assertEquals(param.getValue(), "0.750000");
			else if (param.getTag().compareTo("NeighborhoodRule") == 0)
				assertEquals(param.getValue(), "0.500000");
			else if (param.getTag().compareTo("FinalMaxHet") == 0)
				assertEquals(param.getValue(), "0.900000");
			else if (param.getTag().compareTo("ModelType") == 0)
				assertEquals(param.getValue(), "0");
			else if (param.getTag().compareTo("PerfectCallThreshold") == 0)
				assertEquals(param.getValue(), "2.000000");
			else if (param.getTag().compareTo("TotalLLR") == 0)
				assertEquals(param.getValue(), "75.000000");
			else if (param.getTag().compareTo("StrandLLR") == 0)
				assertEquals(param.getValue(), "0.000000");
			else if (param.getTag().compareTo("AberrantSNR2") == 0)
				assertEquals(param.getValue(), "20.000000");
			else if (param.getTag().compareTo("WeakSignal") == 0)
				assertEquals(param.getValue(), "20.000000");
			else if (param.getTag().compareTo("NoSignal") == 0)
				assertEquals(param.getValue(), "1.000000");
		}

		ResequencingResults p = chp.getResequencingResults();

		int n = p.getCalledBasesSize();
		assertEquals(n, 29270);
		assertEquals(p.getCalledBase(0), 'c');
		assertEquals(p.getCalledBase(1), 'c');
		assertEquals(p.getCalledBase(2), 'c');
		assertEquals(p.getCalledBase(3), 'a');
		assertEquals(p.getCalledBase(4), 'g');
		assertEquals(p.getCalledBase(5), 't');
		assertEquals(p.getCalledBase(6), 'c');

		assertEquals(p.getCalledBase(n - 1), 'n');
		assertEquals(p.getCalledBase(n - 44), 'a');

		final float eps = 1e-5f;
		assertEquals(p.getScore(0), 90.058533, eps);
		assertEquals(p.getScore(1), 130.95447, eps);
		assertEquals(p.getScore(2), 166.22278, eps);
		assertEquals(p.getScore(3), 143.23492, eps);
		assertEquals(p.getScore(n - 1), 0.0, eps);
		assertEquals(p.getScore(n - 2), 0.0, eps);
		assertEquals(p.getScore(n - 6), 75.936035, eps);

		assertEquals(p.getForceCallsSize(), 0);
		assertEquals(p.getOrigCallsSize(), 0);

	}
}
