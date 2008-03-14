/*
 * FusionCHPLegacyDataTest.java
 * JUnit based test
 *
 * Created on October 13, 2005, 11:00 PM
 */

package affymetrix.fusion.chp;

import java.io.File;
import java.io.IOException;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.fusion.FusionTagValuePair;
import affymetrix.gcos.chp.BackgroundZoneType;
import affymetrix.gcos.chp.ForceCallType;

/**
 * 
 * @author ljevon
 */
public class FusionCHPDataTest extends TestCase {

	public FusionCHPDataTest(String testName) throws Exception {
		super(testName);
		FusionCHPLegacyData.registerReader();
		FusionCHPGenericData.registerReader();
		FusionCHPTilingData.registerReader();
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() throws Exception {
		TestSuite suite = new TestSuite(FusionCHPDataTest.class);
		return suite;
	}

	private static final String ROOT_PATH = "..\\..\\..\\CPPTest\\data\\";

	private static final String EXP_COMP_CHP_FILE = ROOT_PATH + "test.exp.comp.chp";

	private static final String EXP_ABS_CHP_FILE = ROOT_PATH + "test.exp.abs.chp";

	private static final String NO_FILE = ROOT_PATH + "no_file.chp";

	private static final String TAG_V11_FILE = ROOT_PATH + "tag_v11.chp";

	private static final String TAG_XDA_FILE = ROOT_PATH + "tag_XDA.chp";

	private static final String RESEQ_V1_FILE = ROOT_PATH + "reseq_xda_v1.chp";

	private static final String RESEQ_V2_FILE = ROOT_PATH + "reseq_xda_v2.chp";

	private static final String RESEQ_OLD_FILE = ROOT_PATH + "reseq_v13.chp";

	public void testReadCalvinWithParameters() throws Exception {
		String file = new File("..\\..\\AffxFusion\\test\\affymetrix\\fusion\\chp\\Test3-2-121502.calvin.CHP")
				.getCanonicalPath();
		FusionCHPData chp = FusionCHPDataReg.read(file);
		FusionCHPLegacyData data = FusionCHPLegacyData.fromBase(chp);

		assertEquals(data.getHeader().getAlgName(), "ExpressionStat");
		assertEquals(data.getHeader().getAlgVersion(), "5.0");
		assertEquals(data.getHeader().getAlgorithmParameter("Alpha1"), "0.04");
		assertEquals(data.getHeader().getAssayType(), FusionCHPHeader.EXPRESSION_ASSAY);
		assertEquals(data.getHeader().getChipType(), "Test3");
		assertEquals(data.getHeader().getRows(), 126);
		assertEquals(data.getHeader().getCols(), 126);
		assertEquals(data.getHeader().getParentCellFile(), "Test3-2-121502.CEL");
		assertEquals(data.getHeader().getProgID(), "GeneChip.CallGEBaseCall.1");
		assertEquals(data.getHeader().getSummaryParameter("RawQ"), "2.79");

	}

	public void testReadCalvin() throws Exception {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_expression_file");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		FusionCHPLegacyData data = FusionCHPLegacyData.fromBase(chp);

		assertTrue(data.getHeader().getNumProbeSets() == 10);
		assertTrue(data.getHeader().getBackgroundZoneInfo().getNumberZones() == 2);

		FusionExpressionProbeSetResults entry1 = new FusionExpressionProbeSetResults();
		data.getExpressionResults(0, entry1);
		// assertEquals(entry1.getProbeSetName(), "probe set 1");
		assertTrue(entry1.getDetection().toShort() == 10);
		assertTrue(entry1.getDetectionPValue() == 11.0f);
		float signal = entry1.getSignal();
		assertTrue(signal == 17.8f);
		assertTrue(entry1.getNumPairs().toInt() == 6);
		assertTrue(entry1.getNumUsedPairs().toInt() == 5);
		assertTrue(entry1.getChange().toShort() == 2);
		float chgPVal = entry1.getChangePValue();
		assertTrue(chgPVal == 56.9f);
		assertTrue(entry1.getSignalLogRatio() == 45.89f);
		assertTrue(entry1.getSignalLogRatioLow() == 42.0f);
		assertTrue(entry1.getSignalLogRatioHigh() == 47.0f);
		assertTrue(entry1.getNumCommonPairs().toInt() == 2);

		FusionExpressionProbeSetResults entry2 = new FusionExpressionProbeSetResults();
		data.getExpressionResults(1, entry2);
		// assertTrue(entry2.getProbeSetName() == "probe set 2");
		assertTrue(entry2.getDetection().toShort() == 10);
		assertTrue(entry2.getDetectionPValue() == 1.0f);
		assertTrue(entry2.getSignal() == 7.8f);
		assertTrue(entry2.getNumPairs().toInt() == 6);
		assertTrue(entry2.getNumUsedPairs().toInt() == 5);
		assertTrue(entry2.getChange().toShort() == 2);
		assertTrue(entry2.getChangePValue() == 5.9f);
		assertTrue(entry2.getSignalLogRatio() == 4.89f);
		assertTrue(entry2.getSignalLogRatioLow() == 2.0f);
		assertTrue(entry2.getSignalLogRatioHigh() == 7.0f);
		assertTrue(entry2.getNumCommonPairs().toInt() == 2);

		BackgroundZoneType zone1 = data.getHeader().getBackgroundZoneInfo().getZone(0);
		assertTrue(zone1.getCenterX() == 11.0f);
		assertTrue(zone1.getCenterY() == 17.8f);
		assertTrue(zone1.getBackground() == 56.9f);
		// assertTrue(data.getHeader().getBackgroundZoneInfo().getSmoothFactor() == 45.89f);

		BackgroundZoneType zone2 = data.getHeader().getBackgroundZoneInfo().getZone(1);
		assertTrue(zone2.getCenterX() == 1.0f);
		assertTrue(zone2.getCenterY() == 7.8f);
		assertTrue(zone2.getBackground() == 6.9f);
		assertTrue(data.getHeader().getBackgroundZoneInfo().getSmoothFactor() == 5.89f);

	}

	public void testFileId() {
		FusionCHPData chp = FusionCHPDataReg.read(EXP_COMP_CHP_FILE);
		// FusionCHPLegacyData fusionchp = FusionCHPLegacyData.fromBase(chp);
		assertFalse(chp == null);
		assertEquals(chp.getFileId(), null);

	}

	public void testReadHeader_fail() {
		FusionCHPData chp = FusionCHPDataReg.readHeader(NO_FILE);
		assertEquals(chp, null);
	}

	public void testRead_fail() {
		FusionCHPData chp = FusionCHPDataReg.read(NO_FILE);
		assertTrue(chp == null);
	}

	public void testReadExpComp() throws IOException, UnsignedOutOfLimitsException {
		FusionCHPData chp = FusionCHPDataReg.read(new File(EXP_COMP_CHP_FILE).getCanonicalPath());
		FusionCHPLegacyData fusionchp = FusionCHPLegacyData.fromBase(chp);
		assertFalse(chp == null);

		// Accessors for header information.
		assertEquals(fusionchp.getHeader().getCols(), 120);
		assertEquals(fusionchp.getHeader().getRows(), 120);
		assertEquals(fusionchp.getHeader().getNumProbeSets(), 8);
		assertEquals(fusionchp.getHeader().getChipType(), "TestExon");
		assertEquals(fusionchp.getHeader().getAlgName(), "Plier");
		assertEquals(fusionchp.getHeader().getAlgorithmParameter("p1"), "1.1.1");
		assertEquals(fusionchp.getHeader().getAlgorithmParameter("p2"), "2.1.1");
		assertEquals(fusionchp.getHeader().getSummaryParameter("cp1"), "1.2.1");
		assertEquals(fusionchp.getHeader().getSummaryParameter("cp2"), "2.2.1");
		assertEquals(fusionchp.getHeader().getAlgVersion(), "1.0");
		assertEquals(fusionchp.getHeader().getParentCellFile(), "TestExon.cel");
		assertEquals(fusionchp.getHeader().getProgID(), "test_id");

		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo().getNumberZones(), 5);
		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo().getSmoothFactor(), 3.8f, 0.00001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(0, 0).getBackground(), 3.2f, 0.00001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(88, 33).getBackground(), 8.9f, 0.00001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(0, 60).getBackground(), 13.1f, 0.00001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(60, 0).getBackground(), 5.0f, 0.00001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(120, 120).getBackground(), 10.8f, 0.00001f);
		int n = fusionchp.getHeader().getBackgroundZoneInfo().getNumberZones();
		for (int i = 0; i < n; i++) {
			BackgroundZoneType zone = fusionchp.getHeader().getBackgroundZoneInfo().getZone(i);
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

		FusionExpressionProbeSetResults exp = new FusionExpressionProbeSetResults();
		for (int i = 0; i < fusionchp.getHeader().getNumProbeSets(); i++) {
			fusionchp.getExpressionResults(i, exp);
			assertEquals(exp.getDetectionPValue(), (0.05 - (i / 1000.0)), 0.00001f);
			assertEquals(exp.getSignal(), (1.1 + i), 0.00001f);
			assertEquals(exp.getNumPairs().toInt(), 3 + i);
			assertEquals(exp.getNumUsedPairs().toInt(), 2 + i);
			assertEquals(exp.getDetection().toShort(), (i % 4));
			switch (exp.getDetection().toShort()) {
			case (FusionExpressionProbeSetResults.ABS_PRESENT_CALL):
				assertEquals(exp.getDetectionString(), "P");
				break;
			case (FusionExpressionProbeSetResults.ABS_MARGINAL_CALL):
				assertEquals(exp.getDetectionString(), "M");
				break;
			case (FusionExpressionProbeSetResults.ABS_ABSENT_CALL):
				assertEquals(exp.getDetectionString(), "A");
				break;
			case (FusionExpressionProbeSetResults.ABS_NO_CALL):
				assertEquals(exp.getDetectionString(), "No Call");
				break;
			default:
				assertEquals(exp.getDetectionString(), "");
				break;
			}
			assertTrue(exp.hasCompResults());
			assertEquals(exp.getChangePValue(), (0.04f - (i / 1000.0)), 0.0000001f);
			assertEquals(exp.getSignalLogRatio(), (1.1f + i), 0.0000001f);
			assertEquals(exp.getSignalLogRatioLow(), (-1.1f + i), 0.0000001f);
			assertEquals(exp.getSignalLogRatioHigh(), (10.1f + i), 0.0000001f);
			assertEquals(exp.getNumCommonPairs().toInt(), 2 + i);
			assertEquals(exp.getChange().toShort(), (i % 6 + 1));
			switch (exp.getChange().toShort()) {
			case (FusionExpressionProbeSetResults.COMP_INCREASE_CALL):
				assertEquals(exp.getChangeString(), "I");
				break;
			case (FusionExpressionProbeSetResults.COMP_DECREASE_CALL):
				assertEquals(exp.getChangeString(), "D");
				break;
			case (FusionExpressionProbeSetResults.COMP_MOD_INCREASE_CALL):
				assertEquals(exp.getChangeString(), "MI");
				break;
			case (FusionExpressionProbeSetResults.COMP_MOD_DECREASE_CALL):
				assertEquals(exp.getChangeString(), "MD");
				break;
			case (FusionExpressionProbeSetResults.COMP_NO_CHANGE_CALL):
				assertEquals(exp.getChangeString(), "NC");
				break;
			case (FusionExpressionProbeSetResults.COMP_NO_CALL):
				assertEquals(exp.getChangeString(), "No Call");
				break;
			default:
				assertEquals(exp.getChangeString(), "");
				break;
			}
		}
	}

	public void testReadExpAbs() throws Exception {
		FusionCHPData chp = FusionCHPDataReg.read(new File(EXP_ABS_CHP_FILE).getCanonicalPath());
		FusionCHPLegacyData fusionchp = FusionCHPLegacyData.fromBase(chp);
		assertTrue(chp != null);

		// Accessors for header information.
		// assertEquals( fusionchp.getHeader().getVersionNumber(), 1 );
		assertEquals(fusionchp.getHeader().getCols(), 120);
		assertEquals(fusionchp.getHeader().getRows(), 120);
		assertEquals(fusionchp.getHeader().getNumProbeSets(), 8);
		assertEquals(fusionchp.getHeader().getChipType(), "TestExon");
		assertEquals(fusionchp.getHeader().getAlgName(), "Plier");
		assertEquals(fusionchp.getHeader().getAlgorithmParameter("p1"), "1.1.1");
		assertEquals(fusionchp.getHeader().getAlgorithmParameter("p2"), "2.1.1");
		assertEquals(fusionchp.getHeader().getSummaryParameter("cp1"), "1.2.1");
		assertEquals(fusionchp.getHeader().getSummaryParameter("cp2"), "2.2.1");
		assertEquals(fusionchp.getHeader().getAlgVersion(), "1.0");
		assertEquals(fusionchp.getHeader().getParentCellFile(), "TestExon.cel");
		assertEquals(fusionchp.getHeader().getProgID(), "test_id");

		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo().getNumberZones(), 5);
		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo().getSmoothFactor(), 3.8f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(0, 0).getBackground(), 3.2f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(88, 33).getBackground(), 8.9f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(0, 60).getBackground(), 13.1f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(60, 0).getBackground(), 5.0f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(120, 120).getBackground(), 10.8f, 0.000001f);
		int n = fusionchp.getHeader().getBackgroundZoneInfo().getNumberZones();
		for (int i = 0; i < n; i++) {
			BackgroundZoneType zone = fusionchp.getHeader().getBackgroundZoneInfo().getZone(i);
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

		FusionExpressionProbeSetResults exp = new FusionExpressionProbeSetResults();
		for (int i = 0; i < fusionchp.getHeader().getNumProbeSets(); i++) {
			fusionchp.getExpressionResults(i, exp);
			assertEquals(exp.getDetectionPValue(), (0.05 - (i / 1000.0)), 0.000001f);
			assertEquals(exp.getSignal(), (1.1 + i), 0.000001f);
			assertEquals(exp.getNumPairs().toInt(), 3 + i);
			assertEquals(exp.getNumUsedPairs().toInt(), 2 + i);
			assertEquals(exp.getDetection().toShort(), (i % 4));
			switch (exp.getDetection().toShort()) {
			case (FusionExpressionProbeSetResults.ABS_PRESENT_CALL):
				assertEquals(exp.getDetectionString(), "P");
				break;
			case (FusionExpressionProbeSetResults.ABS_MARGINAL_CALL):
				assertEquals(exp.getDetectionString(), "M");
				break;
			case (FusionExpressionProbeSetResults.ABS_ABSENT_CALL):
				assertEquals(exp.getDetectionString(), "A");
				break;
			case (FusionExpressionProbeSetResults.ABS_NO_CALL):
				assertEquals(exp.getDetectionString(), "No Call");
				break;
			default:
				assertEquals(exp.getDetectionString(), "");
				break;
			}
			assertFalse(exp.hasCompResults());
		}
	}

	public void testReadHeaderExpComp() throws Exception {
		FusionCHPData chp = FusionCHPDataReg.readHeader(new File(EXP_COMP_CHP_FILE).getCanonicalPath());
		FusionCHPLegacyData fusionchp = FusionCHPLegacyData.fromBase(chp);
		assertTrue(chp != null);

		// Accessors for header information.
		// assertEquals( fusionchp.getHeader().getVersionNumber(), 1 );
		assertEquals(fusionchp.getHeader().getCols(), 120);
		assertEquals(fusionchp.getHeader().getRows(), 120);
		assertEquals(fusionchp.getHeader().getNumProbeSets(), 8);
		assertEquals(fusionchp.getHeader().getChipType(), "TestExon");
		assertEquals(fusionchp.getHeader().getAlgName(), "Plier");
		assertEquals(fusionchp.getHeader().getAlgorithmParameter("p1"), "1.1.1");
		assertEquals(fusionchp.getHeader().getAlgorithmParameter("p2"), "2.1.1");
		assertEquals(fusionchp.getHeader().getSummaryParameter("cp1"), "1.2.1");
		assertEquals(fusionchp.getHeader().getSummaryParameter("cp2"), "2.2.1");
		assertEquals(fusionchp.getHeader().getAlgVersion(), "1.0");
		assertEquals(fusionchp.getHeader().getParentCellFile(), "TestExon.cel");
		assertEquals(fusionchp.getHeader().getProgID(), "test_id");

		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo().getNumberZones(), 5);
		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo().getSmoothFactor(), 3.8f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(0, 0).getBackground(), 3.2f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(88, 33).getBackground(), 8.9f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(0, 60).getBackground(), 13.1f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(60, 0).getBackground(), 5.0f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(120, 120).getBackground(), 10.8f, 0.000001f);

		int n = fusionchp.getHeader().getBackgroundZoneInfo().getNumberZones();
		for (int i = 0; i < n; i++) {
			BackgroundZoneType zone = fusionchp.getHeader().getBackgroundZoneInfo().getZone(i);
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

	public void testReadHeaderExpAbs() throws Exception {
		FusionCHPData chp = FusionCHPDataReg.readHeader(new File(EXP_ABS_CHP_FILE).getCanonicalPath());
		FusionCHPLegacyData fusionchp = FusionCHPLegacyData.fromBase(chp);
		assertTrue(chp != null);

		// Accessors for header information.
		// assertEquals( fusionchp.getHeader().getVersionNumber(), 1 );
		assertEquals(fusionchp.getHeader().getCols(), 120);
		assertEquals(fusionchp.getHeader().getRows(), 120);
		assertEquals(fusionchp.getHeader().getNumProbeSets(), 8);
		assertEquals(fusionchp.getHeader().getChipType(), "TestExon");
		assertEquals(fusionchp.getHeader().getAlgName(), "Plier");
		assertEquals(fusionchp.getHeader().getAlgorithmParameter("p1"), "1.1.1");
		assertEquals(fusionchp.getHeader().getAlgorithmParameter("p2"), "2.1.1");
		assertEquals(fusionchp.getHeader().getSummaryParameter("cp1"), "1.2.1");
		assertEquals(fusionchp.getHeader().getSummaryParameter("cp2"), "2.2.1");
		assertEquals(fusionchp.getHeader().getAlgVersion(), "1.0");
		assertEquals(fusionchp.getHeader().getParentCellFile(), "TestExon.cel");
		assertEquals(fusionchp.getHeader().getProgID(), "test_id");

		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo().getNumberZones(), 5);
		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo().getSmoothFactor(), 3.8f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(0, 0).getBackground(), 3.2f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(88, 33).getBackground(), 8.9f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(0, 60).getBackground(), 13.1f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(60, 0).getBackground(), 5.0f, 0.000001f);
		assertEquals(fusionchp.getHeader().getBackgroundZone(120, 120).getBackground(), 10.8f, 0.000001f);

		int n = fusionchp.getHeader().getBackgroundZoneInfo().getNumberZones();
		for (int i = 0; i < n; i++) {
			BackgroundZoneType zone = fusionchp.getHeader().getBackgroundZoneInfo().getZone(i);
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

	public void testReadHeaderTagV11() throws Exception {
		FusionCHPData chp = FusionCHPDataReg.readHeader(new File(TAG_V11_FILE).getCanonicalPath());
		FusionCHPLegacyData fusionchp = FusionCHPLegacyData.fromBase(chp);
		assertTrue(chp != null);

		// assertEquals( fusionchp.getHeader().getVersionNumber(), 11 );
		assertEquals(fusionchp.getHeader().getCols(), 105);
		assertEquals(fusionchp.getHeader().getRows(), 105);
		assertEquals(fusionchp.getHeader().getNumProbeSets(), 2050);
		assertEquals(fusionchp.getHeader().getChipType(), "GenFlex");
		assertEquals(fusionchp.getHeader().getAlgName(), "Hybridization");
		assertEquals(fusionchp.getHeader().getAlgVersion(), "4.0");
		assertEquals(fusionchp.getHeader().getParentCellFile(), "c:\\genechip\\testdata\\universal\\1109-A-05a.CEL");
		assertEquals(fusionchp.getHeader().getProgID(), "GeneChipAnalysis.HybBaseCall.1");
		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(fusionchp.getHeader().getSummaryParameters(), null);
		assertEquals(fusionchp.getHeader().getAlgorithmParameters(), null);
	}

	public void testReadTagV11() throws Exception {
		FusionCHPData chp = FusionCHPDataReg.read(new File(TAG_V11_FILE).getCanonicalPath());
		FusionCHPLegacyData fusionchp = FusionCHPLegacyData.fromBase(chp);
		assertTrue(chp != null);

		// assertEquals( fusionchp.getHeader().getVersionNumber(), 11 );
		assertEquals(fusionchp.getHeader().getCols(), 105);
		assertEquals(fusionchp.getHeader().getRows(), 105);
		assertEquals(fusionchp.getHeader().getNumProbeSets(), 2050);
		assertEquals(fusionchp.getHeader().getChipType(), "GenFlex");
		assertEquals(fusionchp.getHeader().getAlgName(), "Hybridization");
		assertEquals(fusionchp.getHeader().getAlgVersion(), "4.0");
		assertEquals(fusionchp.getHeader().getParentCellFile(), "c:\\genechip\\testdata\\universal\\1109-A-05a.CEL");
		assertEquals(fusionchp.getHeader().getProgID(), "GeneChipAnalysis.HybBaseCall.1");
		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(fusionchp.getHeader().getSummaryParameters(), null);
		assertEquals(fusionchp.getHeader().getAlgorithmParameters(), null);

		FusionUniversalProbeSetResults univ = new FusionUniversalProbeSetResults();
		fusionchp.getUniversalResults(0, univ);
		assertEquals(univ.getBackground(), 114.595f, 0.000001f);
		fusionchp.getUniversalResults(1, univ);
		assertEquals(univ.getBackground(), 118.9f, 0.000001f);
		fusionchp.getUniversalResults(2049, univ);
		assertEquals(univ.getBackground(), 114.854f, 0.000001f);
	}

	public void testReadHeaderTagXDA() throws Exception {
		FusionCHPData chp = FusionCHPDataReg.readHeader(new File(TAG_XDA_FILE).getCanonicalPath());
		FusionCHPLegacyData fusionchp = FusionCHPLegacyData.fromBase(chp);
		assertTrue(chp != null);

		// assertEquals( fusionchp.getHeader().getVersionNumber(), 1 );
		assertEquals(fusionchp.getHeader().getCols(), 105);
		assertEquals(fusionchp.getHeader().getRows(), 105);
		assertEquals(fusionchp.getHeader().getNumProbeSets(), 2050);
		assertEquals(fusionchp.getHeader().getChipType(), "GenFlex");
		assertEquals(fusionchp.getHeader().getAlgName(), "Hybridization");
		assertEquals(fusionchp.getHeader().getAlgVersion(), "5.0");
		assertEquals(fusionchp.getHeader().getParentCellFile(),
				"C:\\Program Files\\Affymetrix\\GeneChip\\Affy_Data\\Data\\1109-A-05a.CEL");
		assertEquals(fusionchp.getHeader().getProgID(), "GDMTAnalysis.HybBaseCall.1");
		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(fusionchp.getHeader().getSummaryParameters().size(), 1);
		FusionTagValuePair param = fusionchp.getHeader().getSummaryParameters().get(0);
		assertEquals(param.getTag(), "WaveLength");
		assertEquals(param.getValue(), null);
		assertEquals(fusionchp.getHeader().getAlgorithmParameters(), null);
	}

	public void testReadTagXDA() throws Exception {
		FusionCHPData chp = FusionCHPDataReg.read(new File(TAG_XDA_FILE).getCanonicalPath());
		FusionCHPLegacyData fusionchp = FusionCHPLegacyData.fromBase(chp);
		assertTrue(chp != null);

		// assertEquals( fusionchp.getHeader().getVersionNumber(), 1 );
		assertEquals(fusionchp.getHeader().getCols(), 105);
		assertEquals(fusionchp.getHeader().getRows(), 105);
		assertEquals(fusionchp.getHeader().getNumProbeSets(), 2050);
		assertEquals(fusionchp.getHeader().getChipType(), "GenFlex");
		assertEquals(fusionchp.getHeader().getAlgName(), "Hybridization");
		assertEquals(fusionchp.getHeader().getAlgVersion(), "5.0");
		assertEquals(fusionchp.getHeader().getParentCellFile(),
				"C:\\Program Files\\Affymetrix\\GeneChip\\Affy_Data\\Data\\1109-A-05a.CEL");
		assertEquals(fusionchp.getHeader().getProgID(), "GDMTAnalysis.HybBaseCall.1");
		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(fusionchp.getHeader().getSummaryParameters().size(), 1);
		FusionTagValuePair param = fusionchp.getHeader().getSummaryParameters().get(0);
		assertEquals(param.getTag(), "WaveLength");
		assertEquals(param.getValue(), null);
		assertEquals(fusionchp.getHeader().getAlgorithmParameters(), null);

		FusionUniversalProbeSetResults univ = new FusionUniversalProbeSetResults();
		fusionchp.getUniversalResults(1, univ);
		assertEquals(univ.getBackground(), 114.595f, 0.001f);
		fusionchp.getUniversalResults(2, univ);
		assertEquals(univ.getBackground(), 118.9f, 0.001f);
		fusionchp.getUniversalResults(2049, univ);
		assertEquals(univ.getBackground(), 114.854f, 0.001f);
	}

	public void testReadReseqXDA_v1() throws Exception {
		FusionCHPData chp = FusionCHPDataReg.read(new File(RESEQ_V1_FILE).getCanonicalPath());
		FusionCHPLegacyData fusionchp = FusionCHPLegacyData.fromBase(chp);
		assertTrue(chp != null);

		// assertEquals( fusionchp.getHeader().getVersionNumber(), 1 );
		assertEquals(fusionchp.getHeader().getCols(), 488);
		assertEquals(fusionchp.getHeader().getRows(), 639);
		assertEquals(fusionchp.getHeader().getNumProbeSets(), 2);
		assertEquals(fusionchp.getHeader().getChipType(), "DCNtagIQr510989");
		assertEquals(fusionchp.getHeader().getAlgName(), "CustomSeq");
		assertEquals(fusionchp.getHeader().getAlgVersion(), "2");
		assertEquals(fusionchp.getHeader().getParentCellFile(), "C:\\Data\\GCOS\\Data\\5303_DCN_01.CEL");
		assertEquals(fusionchp.getHeader().getProgID(), "GDMTAnalysis.VDABaseCall.1");
		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(fusionchp.getHeader().getSummaryParameters(), null);
		assertEquals(fusionchp.getHeader().getAlgorithmParameters().size(), 10);
		FusionTagValuePair param = fusionchp.getHeader().getAlgorithmParameters().get(0);
		assertEquals(param.getTag(), "NoSignal");
		assertEquals(param.getValue(), "1.000000");
		param = fusionchp.getHeader().getAlgorithmParameters().get(1);
		assertEquals(param.getTag(), "WeakSignal");
		assertEquals(param.getValue(), "20.000000");
		param = fusionchp.getHeader().getAlgorithmParameters().get(2);
		assertEquals(param.getTag(), "AberrantSNR2");
		assertEquals(param.getValue(), "20.000000");
		param = fusionchp.getHeader().getAlgorithmParameters().get(3);
		assertEquals(param.getTag(), "StrandLLR");
		assertEquals(param.getValue(), "0.000000");
		param = fusionchp.getHeader().getAlgorithmParameters().get(4);
		assertEquals(param.getTag(), "TotalLLR");
		assertEquals(param.getValue(), "75.000000");
		param = fusionchp.getHeader().getAlgorithmParameters().get(5);
		assertEquals(param.getTag(), "PerfectCallThreshold");
		assertEquals(param.getValue(), "2.000000");
		param = fusionchp.getHeader().getAlgorithmParameters().get(6);
		assertEquals(param.getTag(), "ModelType");
		assertEquals(param.getValue(), "0");
		param = fusionchp.getHeader().getAlgorithmParameters().get(7);
		assertEquals(param.getTag(), "FinalMaxHet");
		assertEquals(param.getValue(), "0.900000");
		param = fusionchp.getHeader().getAlgorithmParameters().get(8);
		assertEquals(param.getTag(), "NeighborhoodRule");
		assertEquals(param.getValue(), "0.500000");
		param = fusionchp.getHeader().getAlgorithmParameters().get(9);
		assertEquals(param.getTag(), "SampleReliability");
		assertEquals(param.getValue(), "0.750000");

		FusionResequencingResults p = new FusionResequencingResults();
		fusionchp.getResequencingResults(p);

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

	public void testReadReseqXDA_v2() throws Exception {
		FusionCHPData chp = FusionCHPDataReg.read(new File(RESEQ_V2_FILE).getCanonicalPath());
		FusionCHPLegacyData fusionchp = FusionCHPLegacyData.fromBase(chp);
		assertTrue(chp != null);

		// assertEquals( fusionchp.getHeader().getVersionNumber(), 2 );
		assertEquals(fusionchp.getHeader().getCols(), 960);
		assertEquals(fusionchp.getHeader().getRows(), 1008);
		assertEquals(fusionchp.getHeader().getNumProbeSets(), 2);
		assertEquals(fusionchp.getHeader().getChipType(), "Tristezar520098");
		assertEquals(fusionchp.getHeader().getAlgName(), "Reseq2");
		assertEquals(fusionchp.getHeader().getAlgVersion(), "2.0");
		assertEquals(fusionchp.getHeader().getParentCellFile(),
				"C:\\GeneChip\\Affy_Data\\Data\\Hyb01004 CTV-T36 Expt 1926.CEL");
		assertEquals(fusionchp.getHeader().getProgID(), "GDMTAnalysis.Reseq2BaseCall.1");
		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(fusionchp.getHeader().getSummaryParameters(), null);
		assertEquals(fusionchp.getHeader().getAlgorithmParameters().size(), 1);
		FusionTagValuePair param = fusionchp.getHeader().getAlgorithmParameters().get(0);
		assertEquals(param.getTag(), "QualityScore");
		assertEquals(param.getValue(), "3");

		FusionResequencingResults p = new FusionResequencingResults();
		fusionchp.getResequencingResults(p);

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

		FusionForceCallType force;

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

	public void testReadReseqOldFile() throws Exception {
		FusionCHPData chp = FusionCHPDataReg.read(new File(RESEQ_OLD_FILE).getCanonicalPath());
		FusionCHPLegacyData fusionchp = FusionCHPLegacyData.fromBase(chp);
		assertTrue(chp != null);

		// assertEquals( fusionchp.getHeader().getVersionNumber(), 13 );
		assertEquals(fusionchp.getHeader().getCols(), 488);
		assertEquals(fusionchp.getHeader().getRows(), 639);
		assertEquals(fusionchp.getHeader().getNumProbeSets(), 0);
		assertEquals(fusionchp.getHeader().getChipType(), "DCNtagIQr510989");
		assertEquals(fusionchp.getHeader().getAlgName(), "CustomSeq");
		assertEquals(fusionchp.getHeader().getAlgVersion(), "1");
		assertEquals(fusionchp.getHeader().getParentCellFile(),
				"S:\\GDAS2\\Test_Files\\Resequence\\DCN Validation Data\\5303_DCN_01.cel");
		assertEquals(fusionchp.getHeader().getProgID(), "GDMTAnalysis.VDABaseCall.1");
		assertEquals(fusionchp.getHeader().getBackgroundZoneInfo(), null);
		assertEquals(fusionchp.getHeader().getSummaryParameters(), null);
		assertEquals(fusionchp.getHeader().getAlgorithmParameters().size(), 10);
		for (int i = 0; i < fusionchp.getHeader().getAlgorithmParameters().size(); i++) {
			FusionTagValuePair param = fusionchp.getHeader().getAlgorithmParameters().get(i);
			if (param.getTag().equals("SampleReliability"))
				assertEquals(param.getValue(), "0.750000");
			else if (param.getTag().equals("NeighborhoodRule"))
				assertEquals(param.getValue(), "0.500000");
			else if (param.getTag().equals("FinalMaxHet"))
				assertEquals(param.getValue(), "0.900000");
			else if (param.getTag().equals("ModelType"))
				assertEquals(param.getValue(), "0");
			else if (param.getTag().equals("PerfectCallThreshold"))
				assertEquals(param.getValue(), "2.000000");
			else if (param.getTag().equals("TotalLLR"))
				assertEquals(param.getValue(), "75.000000");
			else if (param.getTag().equals("StrandLLR"))
				assertEquals(param.getValue(), "0.000000");
			else if (param.getTag().equals("AberrantSNR2"))
				assertEquals(param.getValue(), "20.000000");
			else if (param.getTag().equals("WeakSignal"))
				assertEquals(param.getValue(), "20.000000");
			else if (param.getTag().equals("NoSignal"))
				assertEquals(param.getValue(), "1.000000");
		}

		FusionResequencingResults p = new FusionResequencingResults();
		fusionchp.getResequencingResults(p);

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
